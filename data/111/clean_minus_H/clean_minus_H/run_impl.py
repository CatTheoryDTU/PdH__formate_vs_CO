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
from gpaw.poisson import PoissonSolver
from gpaw import GPAW, FermiDirac, Mixer, MixerSum, MixerDif, Davidson
from gpaw.solvation.sjm import SJM, SJMPower12Potential
from gpaw.utilities import h2gpts

# Import solvation modules
from ase.data.vdw import vdw_radii
from gpaw.solvation import (
    EffectivePotentialCavity,
    LinearDielectric,
    GradientSurface,
    SurfaceInteraction)
#==========================================================

filename='CO2'
run=0

# Solvent parameters
u0 = 0.180  # eV
epsinf = 78.36  # Dielectric constant of water at 298 K
gamma = 0.0 #0114843767916  # 18.4*1e-3 * Pascal* m
T = 298.15   # K

def atomic_radii(atoms):
    return [vdw_radii[n] for n in atoms.numbers]

# Structure is created
atoms = read('out.traj')
c1 = atoms.constraints
cell = np.diag(atoms.cell)
if cell[2]<20.0:
    atoms.cell[2,2]+=5.0
#write('test.traj',atoms)
atoms.pbc = (True, True, False)

# SJM parameters
pot = 3.4
ne = 0.0
dpot = 0.01

# File naming
namedict={'filename':filename,'run':run,'pot':pot}
restartfile="{filename}_run{run}_{pot}.pckl"
txtfile="{filename}_run{run}_{pot}.txt"
logfile="{filename}_run{run}_{pot}.log"
trajfile="{filename}_run{run}_{pot}.traj"
gpwfile="{filename}_run{run}_{pot}.gpw"
dens="density_run{run}_{pot}.cube"

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

# Run the calculation
atoms.set_calculator(calc)
# atoms.get_potential_energy()

for i in range(11):
    calc.set(potential=pot)
    calc.set(txt=txtfile.format(**namedict))
    dyn=BFGS(atoms,restart=restartfile.format(**namedict),trajectory=trajfile.format(**namedict),logfile=logfile.format(**namedict),maxstep=0.1)
    dyn.run(fmax=0.05)
    #calc.write(gpwfile)
    density = calc.get_all_electron_density() * Bohr**3
    write(dens.format(**namedict), atoms, data=density)
    pot=round(pot+0.1,1)
    namedict.update(pot=pot)

