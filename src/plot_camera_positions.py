import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation

def get_data_at_time_t(xs, ys, zs, t):
    new_xs = []
    new_ys = []
    new_zs = []
    for i in range(len(xs)):
        if i <= t:
            new_xs.append(-xs[i])
            new_ys.append(ys[i]*2)
            new_zs.append(-zs[i])
        else:
            new_xs.append(-100000)
            new_ys.append(-100000)
            new_zs.append(-100000)
    return new_xs, new_ys, new_zs

OUTPUT_PATH = "/home/enrico/Progetti/Nut_SQPnP/resources/camera_positions.txt"
MEAN_SIZE = 40
f = open(OUTPUT_PATH, "r")
lines = f.readlines()

row_xs = []
row_ys = []
row_zs = []

for line in lines:
    splitted_line = line.split()
    tmp = []
    row_xs.append(abs(float(splitted_line[0])))
    row_ys.append(abs(float(splitted_line[1])))
    row_zs.append(abs(float(splitted_line[2])))

xs = []
ys = []
zs = []
for i in range(len(row_xs)-MEAN_SIZE):
    xs.append(sum(row_xs[i:i+MEAN_SIZE]) / MEAN_SIZE)
    ys.append(sum(row_ys[i:i+MEAN_SIZE]) / MEAN_SIZE)
    zs.append(sum(row_zs[i:i+MEAN_SIZE]) / MEAN_SIZE)

print(row_xs)
print(xs)
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
scatter = ax.scatter([], [], [], marker='o')

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
plt.xlim([-300, 0])
plt.ylim([0, 600])
ax.set_zlim(-300, 0)

def animate(i):
    new_xs, new_ys, new_zs = get_data_at_time_t(xs, ys, zs, i)
    scatter._offsets3d = (new_xs, new_ys, new_zs)

anim = animation.FuncAnimation(fig, animate, len(xs), interval=20, blit=False)

plt.show()
