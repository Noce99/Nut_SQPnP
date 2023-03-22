#include <iostream>
#include <eigen3/Eigen/Dense>
#include <vector>
#include <iostream>
#include "pnp_solver.cpp"

int main(){
	std::vector<Eigen::Matrix<double, 3, 1>> points_3D{
		(Eigen::Matrix<double, 3, 1>()  << 1.155000000000000000e+02, 1.500000000000000000e+00, 7.000000000000000000e+00).finished(),
		(Eigen::Matrix<double, 3, 1>()  << 1.155000000000000000e+02, 1.500000000000000000e+00, 3.250000000000000000e+01).finished(),
		(Eigen::Matrix<double, 3, 1>()  << 1.155000000000000000e+02, 2.700000000000000000e+01, 3.250000000000000000e+01).finished(),
		(Eigen::Matrix<double, 3, 1>()  << 1.155000000000000000e+02, 2.700000000000000000e+01, 7.000000000000000000e+00).finished(),
		(Eigen::Matrix<double, 3, 1>()  << 1.155000000000000000e+02, 6.000000000000000000e+00, 4.750000000000000000e+01).finished(),
		(Eigen::Matrix<double, 3, 1>()  << 1.155000000000000000e+02, 6.000000000000000000e+00, 6.450000000000000000e+01).finished(),
		(Eigen::Matrix<double, 3, 1>()  << 1.155000000000000000e+02, 2.300000000000000000e+01, 6.450000000000000000e+01).finished(),
		(Eigen::Matrix<double, 3, 1>()  << 1.155000000000000000e+02, 2.300000000000000000e+01, 4.750000000000000000e+01).finished(),
	};
	std::vector<Eigen::Matrix<double, 2, 1>> points_2D{
		(Eigen::Matrix<double, 2, 1>()  << 6.030000000000000000e+02, 2.400000000000000000e+02).finished(),
		(Eigen::Matrix<double, 2, 1>()  << 8.095618286132812500e+02, 2.501897583007812500e+02).finished(),
		(Eigen::Matrix<double, 2, 1>()  << 7.930357666015625000e+02, 4.506818237304687500e+02).finished(),
		(Eigen::Matrix<double, 2, 1>()  << 5.923986206054687500e+02, 4.430352478027343750e+02).finished(),
		(Eigen::Matrix<double, 2, 1>()  << 9.236518554687500000e+02, 2.927422485351562500e+02).finished(),
		(Eigen::Matrix<double, 2, 1>()  << 1.057743896484375000e+03, 2.992045288085937500e+02).finished(),
		(Eigen::Matrix<double, 2, 1>()  << 1.041889892578125000e+03, 4.316576538085937500e+02).finished(),
		(Eigen::Matrix<double, 2, 1>()  << 9.113110961914062500e+02, 4.254056396484375000e+02).finished(),
	};
	Eigen::Matrix<double, 3, 3> K{
		(Eigen::Matrix<double, 3, 3>()  << 9.711433250503055206e+02, 0.000000000000000000e+00, 6.513893362716326010e+02,
										   0.000000000000000000e+00, 9.637184346581425416e+02, 3.631989204633114809e+02,
										   0.000000000000000000e+00, 0.000000000000000000e+00, 1.000000000000000000e+00).finished()
	};
	std::cout << "I have " << points_3D.size() << " 3D points, " << points_2D.size() << " 2D points and a (" << K.rows() << ", " << K.cols() << ") K" << std::endl;
	PnPSolver PnP = PnPSolver(points_3D, points_2D, K);
	if (PnP.Solve()){
		std::cout << "Best solution found:" << std::endl;  
		std::cout << PnP.get_best_solution() << std::endl;
	}
}
