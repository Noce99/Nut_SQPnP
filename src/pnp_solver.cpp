#include <iostream>
#include <eigen3/Eigen/Dense>
#include <vector>
#include <iostream>

#define MAX_SINGULAR_VALUE_TO_BEING_NULL 10e-05
#define PERTURBATION_TOLLERANCE 1e-05
#define MAXIMUM_NUMBER_OF_ITERATIONS 15
#define SOLUTION_TOLLERANCE 0.1

class PnPSolver{
	std::vector<Eigen::Matrix<double, 3, 1>> points_3D;
	std::vector<Eigen::Matrix<double, 2, 1>> points_2D;
	Eigen::Matrix<double, 3, 3> K;
	
	int num_of_points;
	std::vector<Eigen::Matrix<double, 2, 1>> normalized_points_2D;
	std::vector<Eigen::Matrix<double, 3, 9>> Ai;
	std::vector<Eigen::Matrix<double, 3, 3>> Qi;
	Eigen::Matrix<double, 3, 3> sum_Qi;
	Eigen::Matrix<double, 3, 9> sum_Qi_Ai;
	Eigen::Matrix<double, 3, 9> P;
	Eigen::Matrix<double, 9, 9> Omega;
	Eigen::Matrix<double, 9, 9> svd_U;
	Eigen::Matrix<double, 9, 1> svd_singular_values;
	int null_space_size;
	std::vector<Eigen::Matrix<double, 9, 1>> eigen_vectors;
	
	std::vector<Eigen::Matrix<double, 9, 1>> solution_r;
	std::vector<double> solution_error;
	
	public:
	PnPSolver(std::vector<Eigen::Matrix<double, 3, 1>> my_points_3D, std::vector<Eigen::Matrix<double, 2, 1>> my_points_2D, Eigen::Matrix<double, 3, 3> my_K){
		points_3D = my_points_3D;
		points_2D = my_points_2D;
		K = my_K;
		num_of_points = my_points_3D.size();
	}
	
	bool Solve(){
		normalize_point_2D();
		computeAi();
		computeQi();
		compute_sum_Qi();
		compute_sum_Qi_Ai();
		computeP();
		computeOmega();
		SVD();
		computeNullSpaceSize();
		computeEigenVector();
		if (null_space_size < 1){
			std::cerr << "Null space size < 1 (" << null_space_size << ")." << std::endl;
			std::cerr << "Singular Values of SVD:" << std::endl;
			std::cerr << svd_singular_values << std::endl;
			return false;
		}
		int k = null_space_size;
		for (int i=0; i<2*k; i++){
			int mu = floor(i/k);
			int ni = 10 - k + i - floor((i+1)/k)*k;
			Eigen::Matrix<double, 9, 1> r;
			if (mu == 0){
				r = nearest_rotation_matrix(sqrt(3)*eigen_vectors[ni]);
			}else{
				r = nearest_rotation_matrix(-sqrt(3)*eigen_vectors[ni]);
			}
			Eigen::Matrix<double, 9, 1> r_hat = SQP(r);
			double epsilon_squared =  r_hat.transpose() * Omega * r_hat;
			solution_r.push_back(r_hat);
			solution_error.push_back(epsilon_squared);
		}
		return checkSolutions();
	}
	
	Eigen::Matrix<double, 4, 4> get_best_solution(){
		int best_solution_index = 0;
		for (int i=1; i<solution_error.size(); i++){
			if (solution_error[i] < solution_error[best_solution_index]){
				best_solution_index = i;
			}
		}
		Eigen::Matrix<double, 9, 1> r = solution_r[best_solution_index];
		Eigen::Matrix<double, 3, 1> t = P * r;
		return (Eigen::Matrix<double, 4, 4>() << r(0), r(1), r(2), t(0),
												 r(3), r(4), r(5), t(1),
												 r(6), r(7), r(8), t(2), 
												 0,    0,    0,    1).finished();
	} 
	
	private:
	void normalize_point_2D(){
		normalized_points_2D.reserve(num_of_points);
		Eigen::Matrix<double, 3, 3> inv_K = K.inverse();
		for (int i=0; i<num_of_points; i++){
			Eigen::Matrix<double, 3, 1> omogenized_point = (Eigen::Matrix<double, 3, 1>() << points_2D[i](0), points_2D[i](1), 1).finished();
			Eigen::Matrix<double, 3, 1> omogenized_normalized_point = inv_K * omogenized_point;
			normalized_points_2D[i] = (Eigen::Matrix<double, 2, 1>() << omogenized_normalized_point(0), omogenized_normalized_point(1)).finished();
		}
	}
	
	void computeAi(){
		Ai.reserve(num_of_points);
		for (int i=0; i<num_of_points; i++){
			Ai[i] = (Eigen::Matrix<double, 3, 9>() <<  points_3D[i](0), points_3D[i](1), points_3D[i](2), 0, 0, 0, 0, 0, 0,
													   0, 0, 0, points_3D[i](0), points_3D[i](1), points_3D[i](2), 0, 0, 0,
													   0, 0, 0, 0, 0, 0, points_3D[i](0), points_3D[i](1), points_3D[i](2)).finished();
		}
	}
	
	void computeQi(){
		Qi.reserve(num_of_points);
		for (int i=0; i<num_of_points; i++){
			Eigen::Matrix<double, 3, 3> tmp = (Eigen::Matrix<double, 3, 3>() << -1, 0, normalized_points_2D[i][0],
																				0, -1, normalized_points_2D[i][1],
																				0,  0,                          0).finished();
			Qi[i] = tmp.transpose() * tmp;
		}
	}
	
	void compute_sum_Qi(){
		sum_Qi = Qi[0];
		for (int i=1; i<num_of_points; i++){
			sum_Qi = sum_Qi + Qi[i];
		}
	}
	
	void compute_sum_Qi_Ai(){
		sum_Qi_Ai = Qi[0] * Ai[0];
		for (int i=1; i<num_of_points; i++){
			sum_Qi_Ai = sum_Qi_Ai + Qi[i] * Ai[i];
		}
	}
	
	void computeP(){
		P = -sum_Qi.inverse() * sum_Qi_Ai;
	}
	
	void computeOmega(){
		Omega = (Ai[0] + P).transpose() * Qi[0] * (Ai[0] + P);
		for (int i=1; i<num_of_points; i++){
			Omega = Omega + (Ai[i] + P).transpose() * Qi[i] * (Ai[i] + P);
		}
	}
	
	void SVD(){
		Eigen::JacobiSVD<Eigen::Matrix<double, 9, 9>> svd(Omega, Eigen::ComputeFullU);
        svd_U = svd.matrixU();
        svd_singular_values = svd.singularValues();
	}
	
	void computeNullSpaceSize(){
		for (int i=num_of_points; i>=0; i--){
			if (svd_singular_values[i] > MAX_SINGULAR_VALUE_TO_BEING_NULL){
				null_space_size = num_of_points - i;
				break;
			}
		}
	}
	
	void computeEigenVector(){
		eigen_vectors.reserve(num_of_points);
		for (int i=1; i<num_of_points; i++){
			eigen_vectors[i] = svd_U.col(i);
		}
	}

	Eigen::Matrix<double, 3, 3> from_raw_to_matrix(Eigen::Matrix<double, 9, 1> x){
		return (Eigen::Matrix<double, 3, 3>() << x(0), x(1), x(2), x(3), x(4), x(5), x(6), x(7), x(8)).finished();
	}
	
	Eigen::Matrix<double, 9, 1> from_matrix_to_raw(Eigen::Matrix<double, 3, 3> X){
		return (Eigen::Matrix<double, 9, 1>() << X(0, 0), X(0, 1), X(0, 2), X(1, 0), X(1, 1), X(1, 2), X(2, 0), X(2, 1), X(2, 2)).finished();
	}
	
	Eigen::Matrix<double, 9, 1> nearest_rotation_matrix(Eigen::Matrix<double, 9, 1> e){
		Eigen::Matrix<double, 3, 3> e_3x3 = from_raw_to_matrix(e);
		Eigen::JacobiSVD<Eigen::Matrix<double, 3, 3>> svd(e_3x3, Eigen::ComputeFullU | Eigen::ComputeFullV);
        Eigen::Matrix<double, 3, 3> svd_U = svd.matrixU();
        Eigen::Matrix<double, 3, 3> svd_V = svd.matrixV();
        double detUV = svd_U.determinant() * svd_V.determinant();
        Eigen::Matrix<double, 3, 3> D = Eigen::Matrix<double, 3, 3>().Zero();
        D(0, 0) = 1;
        D(1, 1) = 1;
        D(2, 2) = detUV;
        Eigen::Matrix<double, 3, 3> r = svd_U * D * svd_V.transpose();
        return from_matrix_to_raw(r);
	}
	
	double detmerninant_of_3x3_raw_matrix(Eigen::Matrix<double, 9, 1> r){
		return r(0)*r(4)*r(8) + r(1)*r(5)*r(6) + r(2)*r(3)*r(7) - r(0)*r(5)*r(7) - r(1)*r(3)*r(8) - r(2)*r(4)*r(6);
	}

	Eigen::Matrix<double, 6, 9> compiute_H(Eigen::Matrix<double, 9, 1> r){
		return (Eigen::Matrix<double, 6, 9>() << 2*r(0),	    2*r(1),		2*r(2),		0, 			0,	 		0,	 		0,	 		0,	 		0,
												 0,			0, 			0, 			2*r(3),		2*r(4), 	    2*r(5), 	0, 			0, 			0,
												 r(3), 		r(4), 		r(5),		r(0), 		r(1), 		r(2), 		0, 			0, 			0,
												 r(6), 		r(7), 		r(8), 		0, 			0, 			0, 			r(0), 		r(1), 		r(2),
												 0, 		0, 			0, 			r(6), 		r(7), 		r(8), 		r(3), 		r(4), 		r(5),
												 r(4)*r(8) - r(5)*r(7), r(5)*r(6) - r(3)*r(8), r(3)*r(7) - r(4)*r(6), r(2)*r(7) - r(1)*r(8), r(0)*r(8) - r(2)*r(6), r(1)*r(6) - r(0)*r(7), r(1)*r(5) - r(2)*r(4), r(2)*r(3) - r(0)*r(5), r(0)*r(4) - r(1)*r(3)).finished();
	}
	
	Eigen::Matrix<double, 6, 1> compiute_h(Eigen::Matrix<double, 9, 1> r){
		return (Eigen::Matrix<double, 6, 1>() << r(0)*r(0) + r(1)*r(1) + r(2)*r(2) - 1,
												 r(3)*r(3) + r(4)*r(4) + r(5)*r(5) - 1,
												 r(0)*r(3) + r(1)*r(4) + r(2)*r(5),
												 r(0)*r(6) + r(1)*r(7) + r(2)*r(8),
												 r(3)*r(6) + r(4)*r(7) + r(5)*r(8),
												 detmerninant_of_3x3_raw_matrix(r) - 1).finished();
	}
	
	Eigen::Matrix<double, 9, 1> SQP(Eigen::Matrix<double, 9, 1> r){
		Eigen::Matrix<double, 9, 1> r_hat = r;
		Eigen::Matrix<double, 9, 1> delta_hat;
		int step = 0;
		do{
			Eigen::Matrix<double, 6, 9> H = compiute_H(r_hat);
			Eigen::Matrix<double, 6, 1> h = compiute_h(r_hat);
			Eigen::Matrix<double, 9, 15> BigMatrix_1 = (Eigen::Matrix<double, 9, 15>() << Omega, H.transpose()).finished();
			Eigen::Matrix<double, 6, 15> BigMatrix_2 = (Eigen::Matrix<double, 6, 15>() << H, Eigen::Matrix<double, 6, 6>().Zero()).finished();
			Eigen::Matrix<double, 15, 15> BigMatrix = (Eigen::Matrix<double, 15, 15>() << BigMatrix_1,
																						  BigMatrix_2).finished();
			Eigen::Matrix<double, 15, 1> BigVector = (Eigen::Matrix<double, 15, 1>() << -Omega*r_hat,
																						-h).finished();
			Eigen::Matrix<double, 15, 1> BigResult = BigMatrix.inverse()*BigVector;
			delta_hat = (Eigen::Matrix<double, 9, 1>() << BigResult.head(9)).finished();
			r_hat = r_hat + delta_hat;
			step = step + 1;
		}while (step <= MAXIMUM_NUMBER_OF_ITERATIONS and delta_hat.norm() > PERTURBATION_TOLLERANCE);
		return r_hat;
	}
	
	bool checkSolutions(){
		for (int i=0; i<solution_error.size(); i++){
			if (solution_error[i] <= SOLUTION_TOLLERANCE){
				return true;
			}
		}
		return false;
	}
};
