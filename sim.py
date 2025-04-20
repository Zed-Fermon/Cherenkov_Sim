import numpy as np 
import math

import matplotlib.pyplot as plt 
import matplotlib.animation as animation

# Define some natural constants

eps_0 = 8.854 * math.pow(10, -12)
mu_0 = 1.256 * math.pow(10, -6)
c = 3 * math.pow(10, 8)

# Define space and time increments and coordinate mesh

dt = .001
dx, dy, dz = .01, .01, .01

x_coord = np.arange(0,1,dx)
y_coord = np.arange(0,1,dy)
z_coord = np.arange(0,1,dz)
(meshX, meshY, meshZ) = np.meshgrid(x_coord, y_coord, z_coord, indexing='xy')

# Define plot objects
fig, ax = plt.subplots()

# Define material and particle properties
eps_rel = 1.33
mu_rel = 1.33
index_ref = np.sqrt(eps_rel*mu_rel)

v_light = c/index_ref
v_part = v_light*1.2

part_charge = 1.602 * math.pow(10,-19)


# Define scaled parameters such that grid is 1 light-second (in vacuum) squared

v_part_scaled = v_part/c

#	Helper functions for animation
def update(frame):
	redraw_figure_E(fig, ax, frame)

def redraw_figure_E(fig, ax, Efield):
	(Ex, Ey) = Efield 
	ax.clear()
	ax.quiver(x_coord, y_coord, Ex, Ey, angles='xy', scale=(part_charge/(4*np.pi*eps_0*eps_rel))*100000)
	fig.canvas.draw()

def redraw_figure_B(fig, ax, B):
	ax.clear()
	ax.matshow(B)
	fig.canvas.draw()


def get_particle_Efield(particle_location):
	#	Returns the X and Y components of the electric field of a point source at the given X, Y

	#Unpack particle location
	(partX, partY, partZ) = particle_location

	#	Start by defining grids for electric field X and Y components
	Ex, Ey, Ez = np.zeros_like(meshX), np.zeros_like(meshY), np.zeros_like(meshZ)

	#	Find the squared distance of each grid point from the source
	d2 = (meshX - partX)**2 + (meshY - partY)**2 + (meshZ - partZ)**2

	#	Find the angle of rhat in xy plane at each point in the grid
	rhat_angle = np.atan2((partY-meshY), (partX-meshX))+(np.pi)
	
	#	Calculate electric field magnitude using that distance
	E = part_charge/(4*np.pi*eps_0*eps_rel*d2)

	#	Calculate field components and return
	Ex = E*np.cos(rhat_angle)
	Ey = E*np.sin(rhat_angle)
	return Ex, Ey, Ez

def get_Bz(Ex, Ey):
	#	Returns z-component of magnetic field
	return (np.gradient(Ex, axis=1)-np.gradient(Ey, axis=0))*dt*(-1)


def main():
	part_location = (.1, .5, .5)

	Eframes, Bframes = list(), list()

	Ex, Ey, Ez= get_particle_Efield(part_location)
	Bz = get_Bz(Ex, Ey)
	Eframes.append((Ex, Ey))
	Bframes.append(Bz)
	
	while part_location[0]<1:
		#	loop until particle leaves frame

		#	start by moving the particle
		part_location = (part_location[0] + dt*v_part_scaled, part_location[1], part_location[2])

		#	Compute new E field
		Ex, Ey, Ez = get_particle_Efield(part_location)

		#	Compute new B field
		Bz = get_Bz(Ex, Ey)

		Eframes.append((Ex, Ey))
		Bframes.append(Bz)

	#	Show animation
	ani = animation.FuncAnimation(fig, update, frames=Eframes, interval=10, repeat_delay=1000)
	#ani.save("Cherenkov_sim_animation.mp4")
	plt.show()

if __name__ == '__main__':
	main()