import sys
import os
import matplotlib.pyplot as plt
import numpy as np
sys.path.append(f'{os.path.dirname(__file__)}/..')
from modules.kinematics.velocity_gradient import VelocityGradient


kin = VelocityGradient(0,0,0,60)
x = np.linspace(-10,0)
y = 0
z = 0
u,v,w = kin.draw_random_velocity(x,y,z)

plt.plot(x,v)
plt.show()