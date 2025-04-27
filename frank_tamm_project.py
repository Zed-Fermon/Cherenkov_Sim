import numpy as np
import matplotlib.pyplot as plt
import colorpy.plots

#	Constants:
c = 3*(10**8)
eps0 = 8.854*(10**-12)
mu0 = 1.256*(10**-6)
q = 1.602*(10**-19)

def main():
	base_spec = get_ft_spectrum()
	colorpy.plots.spectrum_plot(base_spec, 'Emitted Radiation', 'Cherenkov_ft_plot', 'Wavelength (nm)', 'Intensity per unit length')
	colorpy.plots.spectrum_plot(do_raleigh_scattering(base_spec), 'Raleigh Scattered Emitted Radiation', 'Scattered_Cherenkov_plot', 'Wavelength (nm)', 'Intensity per unit length')

def plot_basic_spectrum():
	fig, ax = plt.subplots()
	w_min = 400*(10**12)*(2*np.pi)
	w_max = 790*(10**12)*(2*np.pi)
	dw = round((w_max-w_min)/1000)
	#print(w_min, w_max, dw)
	w = np.arange(w_min, w_max, dw)
	#print(w)
	n = get_n(w)
	mu = get_mu(w)
	v = get_v()

	spec = frank_tamm(n, mu, v, q, w)
	ax.plot(w, spec)
	ax.set_xlabel('Angular Frequency')
	ax.set_ylabel('Energy per particle unit length travelled')
	plt.show()

def frank_tamm(n, mu, v, q, w):
	return ((q**2)/(4*np.pi))*mu*w*(1-((c**2)/((v**2)*(n**2))))

def get_mu(w):
	return mu0*.999992

def get_n(w):
	return 1.33

def get_v():
	return (c/1.33)*1.3

def get_ft_spectrum():
	spec = empty_spectrum()
	(rows, cols) = spec.shape
	for i in range(rows):
		#	convert wavelength to angular frequency
		w = (2*np.pi*get_v())/spec[i][0]
		spec[i][1] = frank_tamm(get_n(w), get_mu(w), get_v(), q, w)*(10**40)
	return spec

def do_raleigh_scattering(spec):
	(rows, cols) = spec.shape
	for i in range(rows):
		spec[i][1] = spec[i][1]/((spec[i][0]*(10**-9))**4)
	return spec

def empty_spectrum():
	spec = np.zeros((831-360, 2))
	for i, wl in enumerate(range(360, 831)):
		spec[i][0] = wl
	return spec

if __name__ == '__main__':
	main()