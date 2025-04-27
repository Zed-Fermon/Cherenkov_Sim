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
dx, dy = .005, .005
dz = .01

x_coord = np.arange(0,1,dx)
y_coord = np.arange(0,1,dy)
#z_coord = np.arange(0,1,dz)
(meshX, meshY) = np.meshgrid(x_coord, y_coord, indexing='xy')

# Define plot objects
fig, ax = plt.subplots()

# Define material and particle properties
eps_rel = 1.77
mu_rel = 1
index_ref = np.sqrt(eps_rel*mu_rel)

v_light = c/index_ref
global v_part
global v_part_scaled
#v_part = v_light*1.2
#v_part_scaled = v_part/c

part_charge = 1.602 * math.pow(10,-19)
#part_charge = 1


# Define scaled parameters such that grid is 1 light-second (in vacuum) squared

v_light_scaled = v_light/c

#	Helper functions for animation
def update(frame):
	(Ex, Ey) = frame
	redraw_figure_E(fig, ax, (Ex, Ey))
	#redraw_figure_B(fig, ax[1], Bz)
	fig.canvas.draw()

def update2(frame):
	ax.clear()
	ax.imshow(frame)
	fig.canvas.draw()

#def update_B(frame):
#	redraw_figure_B(fig, ax[1], frame)

def redraw_figure_E(fig, ax, Efield):
	(Ex, Ey) = Efield 
	ax.clear()
	ax.quiver(x_coord, y_coord, Ex, Ey, angles='xy', scale=(part_charge/(4*np.pi*eps_0*eps_rel))*150000)
	#fig.canvas.draw()

def redraw_figure_B(fig, ax, B):
	ax.clear()
	ax.imshow(B)
	#fig.canvas.draw()


def get_particle_Efield(path):
	#	Returns the X and Y components of the electric field of a point source at the given X, Y

	fields = list()

	for loc in path:
		#Unpack particle location
		(partX, partY, partZ) = loc

		#	Start by defining grids for electric field X and Y components
		Ex, Ey = np.zeros_like(meshX), np.zeros_like(meshY)

		#	Find the squared distance of each grid point from the source
		d2 = (meshX - partX)**2 + (meshY - partY)**2

		#	Find the angle of rhat in xy plane at each point in the grid
		rhat_angle = np.atan2((partY-meshY), (partX-meshX))+(np.pi)
		
		#	Calculate electric field magnitude using that distance
		E = part_charge/(4*np.pi*eps_0*eps_rel*d2)

		#	Calculate field components and return
		Ex = E*np.cos(rhat_angle)
		Ey = E*np.sin(rhat_angle)
		fields.append((Ex, Ey))
	return fields

def get_Bz(Ex, Ey):
	#	Returns z-component of magnetic field
	return (np.gradient(Ex, axis=0)-np.gradient(Ey, axis=1))*dt*(-1)

def get_particle_path(starting_location):
	#	Returns a list of the particle's location at each timestep

	#Unpack particle location
	(partX, partY, partZ) = starting_location

	path = list()
	path.append(starting_location)

	while(partX < .9):
		partX = partX + v_part_scaled*dt 
		path.append((partX, partY, partZ))

	return path

def get_masked_field(field, location, num_timesteps_passed):
	#	Takes the input field, produced by particle at input location, and returns only the field that would have propogated after however many timesteps
	(partX, partY, partZ) = location
	d = np.sqrt((meshX - partX)**2 + (meshY - partY)**2)

	#	How many timesteps it would take to get to each point
	ret_d = d/(v_light_scaled*dt)

	#	Mask any points that the radiation hasn't reached yet
	field = np.where(ret_d<num_timesteps_passed, field, 0)

	#	Mask any points that the radiation passed too long ago
	prop_t = (dx/v_part_scaled)/dt
	field = np.where(ret_d>(num_timesteps_passed-prop_t), field, 0)
	return field

def do_the_thing(animation_file_name):
	#	First, make an array of all the locations the particle will be for its entire path
	path = get_particle_path((.1, .5, .5))
	#	Then calculate the field of a point source at each of those locations, split into X and Y components
	fields = get_particle_Efield(path)

	frames = list()
	#	Now, loop over all the locations in that path array
	for t_curr, location in enumerate(path):
		#	Set up grids for this timestep
		Ex = np.zeros_like(fields[0][0])
		Ey = np.zeros_like(fields[0][1])
		for t_past, location_past in enumerate(path[0:t_curr]):
			#	This is where all the work is done: for each timestep, iterate through all the locations the particle has already been
			#	and sum up the contributions of all the past fields, but only where the influence of that past time would be 
			#	at the current time
			#	You can imagine concentric rings, where the ring surrounding the particle's current position at timestep t_curr has
			#	the electric field produced by the particle at location(t_curr), one ring out has the electric field produced by the 
			#	particle at location(t_curr-1), and so on
			timesteps_passed = t_curr-t_past
			#test = test + get_masked_field(np.ones_like(fields[0][0]), location_past, timesteps_passed)
			Ex = Ex + get_masked_field(fields[t_past][0], location_past, timesteps_passed)
			Ey = Ey + get_masked_field(fields[t_past][1], location_past, timesteps_passed)
		print(f"Working on frame: {t_curr} out of {len(path)}")
		frames.append((Ex, Ey))
		#frames.append(test)
	#print(frames)
	ani = animation.FuncAnimation(fig, update, frames=frames, interval=50, repeat_delay=1000)
	ani.save(animation_file_name)


def main():
	global v_part
	global v_part_scaled
	v_part = v_light*.8
	v_part_scaled = v_part/c
	do_the_thing('Cherenkov_sim_v.8.mp4')
	v_part = v_light*1.2
	v_part_scaled = v_part/c
	do_the_thing('Cherenkov_sim_v1.2.mp4')


if __name__ == '__main__':
	main()