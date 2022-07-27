import sys,os
import time
import numpy as np
#from sco2_tools.sco2_neb_tools import neb_calculator
from ase.visualize import view
from gpaw.utilities import h2gpts
from gpaw import GPAW, FermiDirac
from ase.io import write, read
from ase.units import Pascal, m, Bohr
# from ase.neb import NEB
from ase.vibrations import Vibrations

# Import solvation modules
from ase.data.vdw import vdw_radii
from gpaw.solvation import (LinearDielectric, GradientSurface,
                            SurfaceInteraction, EffectivePotentialCavity)

from gpaw.solvation.sjm import SJM, SJMPower12Potential

# Solvent parameters from JCP 141, 174108 (2014):
u0 = 0.180  # eV
epsinf = 78.36  # dielectric constant of water at 298 K
gamma = 18.4*1e-3 * Pascal * m  # Surface tension
T = 298.15  # Temperature
atomic_radii = lambda atoms: [vdw_radii[n] for n in atoms.numbers]

home=os.getcwd()
path_info = home.split('/')
images=read('neb_band.traj@:')

highest_point = np.argsort([atoms.get_potential_energy() for atoms in images])[-1]
#print(highest_point)
atoms=images[highest_point]

potential = atoms.calc.results['electrode_potential']
ne = atoms.calc.results['ne']

cell_len = np.linalg.norm(atoms.get_cell(),axis=1)
calc = SJM(poissonsolver={'dipolelayer': 'xy'},
                   doublelayer={'upper_limit': 26},
                   symmetry={'point_group': False},
                   gpts=h2gpts(0.18,atoms.get_cell(),idiv=8),
                   kpts=(np.around(35/cell_len[0],0),np.around(35/cell_len[1],0),1),
                   xc='BEEF-vdW',
                   spinpol=False,
                   potential=potential,
                   dpot=0.005,
                   ne=ne,
                   txt='vibs.txt',
                   occupations=FermiDirac(0.1),
                   maxiter=400,
                   cavity=EffectivePotentialCavity(
                       effective_potential=SJMPower12Potential(
                           atomic_radii, u0, H2O_layer=False),
                       temperature=T,
                       surface_calculator=GradientSurface()),
                   dielectric=LinearDielectric(epsinf=epsinf),
                   interactions=[SurfaceInteraction(surface_tension=gamma)],
                   )

atoms.set_calculator(calc)

vibind=[95,96,97,98]
vib = Vibrations(atoms,indices=vibind)
vib.clean(empty_files=True)

vib.run()
vib.summary()

