import os,sys
import numpy as np
from ase.visualize import view
from ase.io import read, write
from ase.units import mol, kJ, kcal, Pascal, m, Bohr
from ase.vibrations import Vibrations
from ase.parallel import world
from ase.thermochemistry import HarmonicThermo
from ase import units

home=os.getcwd()
path_info=home.split('/')
atoms=read('neb_band.traj@3')
vib=Vibrations(atoms,indices=[95,96,97,98])
vib.summary()
vib.write_mode(0)
out=open('vibrations.out','w')
out.write(' '.join([str(i) for i in vib.get_frequencies()]))
out.close()

vibenergies = vib.get_energies()

thermo = HarmonicThermo(vib_energies=vibenergies[2:],
                        potentialenergy=0)
G = thermo.get_helmholtz_energy(temperature=298.15)


sys.exit()
