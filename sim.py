import numpy as np 
import math

# Define some natural constants

eps_0 = 8.854 * math.pow(10, -12)
mu_0 = 1.256 * math.pow(10, -6)
c = 3 * math.pow(10, 8)

# Define space and time increments and coordinate mesh

dt = .001
dx, dy = .01, .01

x_coord = np.arange(0,1,dx)
y_coord = np.arange(0,1,dy)
(meshX, meshY) = np.meshgrid(x_coord, y_coord)

# Define material and particle properties

eps_rel = 1.33
mu_rel = 1.33
index_ref = np.sqrt(eps_rel*mu_rel)

v_light = c/index_ref
v_part = v_light*1.2

print(v_part)