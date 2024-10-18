import numpy as np
import matplotlib.pyplot as plt

import position

pos1 = position.Position(l_deg=1, b_deg=15, solid_angle=1e-1)

N = 1000
st_dir = np.random.uniform(0, 2 * np.pi, size=N)  # Phi in the paper
st_rad = np.sqrt(np.random.uniform(0, pos1.cone_angle ** 2, size=N))  # Theta in the paper
delta_l_rad = st_rad * np.sin(st_dir)
delta_b_rad = st_rad * np.cos(st_dir)

import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as axes3d

fig = plt.figure()
ax = fig.add_subplot(1, 2, 1, projection='3d')

for l in range(0, 360, 5):
    for b in range(-90, 90, 5):
        pos = position.Position(l_deg=l, b_deg=b, solid_angle=1e-6)
        _, _, _, _, starl_deg, starb_deg = pos.draw_random_point_in_slice(0, 1, 1)
        starl_rad = starl_deg / 180 * np.pi
        starb_rad = starb_deg / 180 * np.pi
        ax.scatter3D(np.cos(starb_rad) * np.cos(starl_rad), np.cos(starb_rad) * np.sin(starl_rad),
            np.sin(starb_rad), c='k', alpha=0.2)

ax.set_xlim(-1, 1)
ax.set_ylim(-1, 1)
ax.set_zlim(-1, 1)

plt.title("this should be look like a sphere")

ax = fig.add_subplot(1, 2, 2)

l = np.random.uniform(-180, 180)
pos = position.Position(l_deg=l, b_deg=0, solid_angle=1e-3)
_, _, _, _, starl_deg, starb_deg = pos.draw_random_point_in_slice(0, 1, 10000)
ax.scatter(starl_deg, starb_deg)
pos = position.Position(l_deg=l, b_deg=0, solid_angle=1e-4)
_, _, _, _, starl_deg, starb_deg = pos.draw_random_point_in_slice(0, 1, 10000)
ax.scatter(starl_deg, starb_deg)

plt.title("this should be look like a big  and a small circle ")
ax.axis('equal')

plt.show()
