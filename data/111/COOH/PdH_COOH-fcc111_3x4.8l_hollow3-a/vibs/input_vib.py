#!/usr/bin/env python

import os, sys
import numpy as np
asa = np.asarray
from math import sqrt, pi
from ase.io import write, read
from gpaw import GPAW, FermiDirac
from gpaw import PW
from ase.optimize import BFGS, FIRE
from ase import Atom,Atoms
from ase.build import fcc111
from ase.constraints import FixAtoms, FixInternals
from ase.units import Bohr
from ase.vibrations import Vibrations
from gpaw.poisson import PoissonSolver
from gpaw import GPAW, FermiDirac, Mixer, MixerSum, MixerDif, Davidson
from gpaw.solvation.sjm import SJM, SJMPower12Potential
from gpaw.utilities import h2gpts
import ase.units as units

# Import solvation modules
from ase.data.vdw import vdw_radii
from gpaw.solvation import (
    EffectivePotentialCavity,
    LinearDielectric,
    GradientSurface,
    SurfaceInteraction)

filename='CO2'
run=0
pot=4.1
ne = 0.0
dpot = 0.005

# Solvent parameters
u0 = 0.180  # eV
epsinf = 78.36  # Dielectric constant of water at 298 K
gamma = 0.0 #0114843767916  # 18.4*1e-3 * Pascal* m
T = 298.15   # K

def atomic_radii(atoms):
    return [vdw_radii[n] for n in atoms.numbers]

# File naming
# File naming
namedict={'filename':filename,'run':run,'pot':pot}
restartfile="{filename}_run{run}_{pot}.pckl"
txtfile="{filename}_run{run}_{pot}.txt"
logfile="{filename}_run{run}_{pot}.log"
trajfile="{filename}_run{run}_{pot}.traj"
gpwfile="{filename}_run{run}_{pot}.gpw"
dens="density_run{run}_{pot}.cube"

atoms = read('../'+trajfile.format(**namedict))
c1 = atoms.constraints
cell = np.diag(atoms.cell)
#ne=atoms.calc.results['ne']
name = f'{filename}_{pot}_vib'

#Find adsorbate indices
iatm=0
vibind=[]
while iatm < len(atoms):
    if atoms.get_chemical_symbols()[iatm] == 'Au':
        iatm+=1
        continue
    if atoms.get_chemical_symbols()[iatm] == 'O' and iatm + 1 < len(atoms):
        if (atoms.get_chemical_symbols()[iatm+1] == 'H' and
            atoms.get_chemical_symbols()[iatm+2] == 'H'):
            if (np.linalg.norm(atoms.positions[iatm]-atoms.positions[iatm+1]) < 1.2 and
                np.linalg.norm(atoms.positions[iatm]-atoms.positions[iatm+2]) < 1.2):
                iatm+=3
                continue

    vibind.append(iatm)
#    print('1',iatm,atoms.get_chemical_symbols()[iatm])
    iatm+=1

vibind=[]
for i in range(len(atoms)-96):
    vibind.append(i+96)
#print(vibind)
# The calculator
calc = SJM(doublelayer={'upper_limit': 26},
           potential=pot,
           dpot=dpot,
           ne=ne,
           verbose=True,
           #write_grandcanonical_energy=False,
           symmetry={'point_group': False},
           gpts=h2gpts(0.18,atoms.get_cell(),idiv=8),
           poissonsolver={'dipolelayer': 'xy'},
           kpts=(np.around(35/cell[0],0),np.around(35/cell[1],0),1),
           xc='BEEF-vdW',
           txt=txtfile.format(**namedict),
           occupations=FermiDirac(0.1),
           cavity=EffectivePotentialCavity(
               effective_potential=SJMPower12Potential(atomic_radii, u0, H2O_layer=False),
               temperature=T,
               surface_calculator=GradientSurface()),
           dielectric=LinearDielectric(epsinf=epsinf),
           interactions=[SurfaceInteraction(surface_tension=gamma)])

cell = np.diag(atoms.cell)
#nel_m = 10.
#atoms, spinpol = add_magnetic_moment(atoms)

#calc = GPAW(
#           symmetry={'point_group': False},
#           gpts=h2gpts(0.18,atoms.get_cell(),idiv=8),
#           poissonsolver={'dipolelayer': 'xy'},
#           kpts=(np.around(35/cell[0],0),np.around(35/cell[1],0),1),
#           xc='BEEF-vdW',
#           txt='out.txt',
#           occupations=FermiDirac(0.1))
# Attach to atoms and get energy/forces
atoms.set_calculator(calc)
#dyn=BFGS(atoms,restart=restartfile.format(**namedict),trajectory=trajfile.format(**namedict),logfile=logfile.format(**namedict),maxstep=0.1)
#dyn.run(fmax=0.03)
potentialenergy = atoms.get_potential_energy()
print(potentialenergy)

atoms_to_vibrate = vibind#[27,28,29]
vib = Vibrations(atoms, indices = atoms_to_vibrate, nfree=2, delta=0.01)
vib.clean(empty_files=True)

vib.run()
vib_energies = vib.get_energies()
vib.summary()
f = open('vib.list','w')
for i in vib_energies:
    if np.real(i)==0:
        f.write(str(np.imag(i / units.invcm)) + 'i' + '\n')
    else:
        f.write(str(np.real(i / units.invcm)) + '\n')
f.close()



vib.summary(method = 'frederiksen', log = open(name+'.txt', mode = 'w'))

# Make trajectory files to visualize normal modes
for i,mode in enumerate(vib.modes):
	vib.write_mode(i)


