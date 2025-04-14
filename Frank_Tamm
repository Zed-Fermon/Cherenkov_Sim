import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import c, eV, h, pi


# Constants
Mass_of_e = 9.31e-31 #kg
Epsilon0 = 8.85e-12 #SI
e = 1.6e-19


class Chernekov_Suite:
    def __init__(self, Omega, speed, Charge):
        ### Add Inputs here
        self.Mass_of_Charge = 9.31e-31 #kg
        self.Epsilon0 = 8.85e-12 #SI
        self.Base_Charge = 1.6e-19
        self.Electron_density = 1 # For water
        self.Omega = Omega
        self.Speed = speed
        self.Charge = Charge
        self.beta = self.Speed/c

        
        
    def Drude(self, Omega):
        
        '''
        How the Permittivity of a material depends on Omega. Assumes only free electrons influence permittivity.
        Inputs: (As stated)
        Outputs: Epsilon as a function of frequency.
        Implementation: Just Read Jackson Bro.... This implementation assumes 1 electron per atom and 1 mole of electrons per volume
        Scope = Plasma Frequency (if needed), Deduce Perimittivity.
        
        '''
        
        Plasma_freq = (self.Electron_density*(self.Base_Charge)**2)/(self.Epsilon0 * self.Mass_of_Charge)
        return self.Epsilon0*(1 - (Plasma_freq/Omega)**2)
    
    def Refractive_Index(self):
        
        Epsilon = self.Drude(self.Omega)
        '''
        Inputs: Value of epsilon from Drude Model
        Outputs: Refractive Index.
        Scope: Add to Frank_Tamm.
        Caveats: Assumes Non-magnetic materials ( We don't have mu dependence on Omega')
        '''
        return np.sqrt(self.Epsilon0/Epsilon)
    
    def Frank_Tamm(self, Omega):
        
        '''
        Inputs: Particle Speed, Refractive Index, Omega
        Output: Power v/s Omega Spectrum
        '''
        Eta = self.Refractive_Index()
        beta = self.beta
        
        return (((self.Charge*self.Base_Charge)/c)**2)*self.speed*(1 - 1/(beta*Eta)**2)*Omega
        
    
    

### Testing Area

Omega = 2*pi*1e14
Speed = 3e8
Charge = 1


Code = Chernekov_Suite(Omega, Speed , Charge)
print(Code.Drude(Omega))
