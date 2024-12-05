"""
Validate the Position module, which generates random positions within the conical field of view.
"""

#External imports
import os
DIRNAME=os.path.dirname(os.path.abspath(__file__))
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as axes3d

#SynthPop imports
from synthpop import position
from synthpop.synthpop_utils.coordinates_transformation import CoordTrans

#first position object
pos1 = position.Position(coord_trans=CoordTrans())
pos1.update_location(l_deg=1, b_deg=15, solid_angle=1e-1)

N = 1000
st_dir = np.random.uniform(0, 2 * np.pi, size=N)  # Phi in the paper
st_rad = np.sqrt(np.random.uniform(0, pos1.cone_angle ** 2, size=N))  # Theta in the paper
delta_l_rad = st_rad * np.sin(st_dir)
delta_b_rad = st_rad * np.cos(st_dir)

fig1 = plt.figure()
ax = fig1.add_subplot(1, 1, 1, projection='3d')

pos = position.Position(coord_trans=CoordTrans())
for l in range(0, 360, 5):
    for b in range(-90, 90, 5):
        pos.update_location(l_deg=l, b_deg=b, solid_angle=1e-6)
        _, _, _, _, starl_deg, starb_deg = pos.draw_random_point_in_slice(0, 1, 1)
        starl_rad = starl_deg / 180 * np.pi
        starb_rad = starb_deg / 180 * np.pi
        ax.scatter3D(np.cos(starb_rad) * np.cos(starl_rad), np.cos(starb_rad) * np.sin(starl_rad),
            np.sin(starb_rad), c='k', alpha=0.2)

ax.set_xlim(-1, 1)
ax.set_ylim(-1, 1)
ax.set_zlim(-1, 1)

plt.title("This should show circles on a sphere.")
plt.tight_layout()

fig2 = plt.figure()
ax = fig2.add_subplot(1, 1, 1)

l = np.random.uniform(-180, 180)
pos.update_location(l_deg=l, b_deg=0, solid_angle=1e-3)
_, _, _, _, starl_deg, starb_deg = pos.draw_random_point_in_slice(0, 1, 10000)
ax.scatter(starl_deg, starb_deg)
pos.update_location(l_deg=l, b_deg=0, solid_angle=1e-4)
_, _, _, _, starl_deg, starb_deg = pos.draw_random_point_in_slice(0, 1, 10000)
ax.scatter(starl_deg, starb_deg,marker=',')

plt.title("This should show a large circle\nand a small circle\ncentered on the same location.")
plt.tight_layout()
ax.axis('equal')

if not os.path.isdir(DIRNAME+'/validation_figures'):
    os.mkdir(DIRNAME+'/validation_figures')
fig1.savefig(DIRNAME+'/validation_figures/'+'position_sphere.png')
fig2.savefig(DIRNAME+'/validation_figures/'+'position_circle.png')
