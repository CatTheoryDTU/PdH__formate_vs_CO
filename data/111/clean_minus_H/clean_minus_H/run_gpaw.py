import os,sys
import numpy as np
from ase.visualize import view
from ase.io import read, write
sys.path.append('/home/cat/geokast/scripts/')
from adsorb_functions import *
from general_tools import *
asa = np.asarray
from math import sqrt, pi
from gpaw import GPAW, FermiDirac
from gpaw import PW
from ase.optimize import BFGS, FIRE
from ase import Atom,Atoms
from ase.build import fcc111
from ase.constraints import FixAtoms, FixInternals
from ase.units import Bohr
from gpaw.poisson import PoissonSolver
from gpaw import GPAW, FermiDirac, Mixer, MixerSum, MixerDif, Davidson
from gpaw.utilities import h2gpts

from ase.data.vdw import vdw_radii
#==========================================================

home = os.getcwd()
path_info = home.split('/')

atoms = read('init.traj')
atoms.center(axis=2,vacuum=8.)
cell = np.diag(atoms.cell)
#nel_m = 10.
atoms, spinpol = add_magnetic_moment(atoms)

calc = GPAW(
           symmetry={'point_group': False},
           gpts=h2gpts(0.18,atoms.get_cell(),idiv=8),
           poissonsolver={'dipolelayer': 'xy'},
           kpts=(np.around(35/cell[0],0),np.around(35/cell[1],0),1),
           xc='BEEF-vdW',
           txt='out.txt',
           occupations=FermiDirac(0.1))


atoms.set_calculator(calc)
calc.initialize(atoms)
#nelect0 = calc.default_nelect_from_ppp()

#outdir = 'output/'
#create_dir_and_save_old_output(outdir)

#os.system('cp run_gpaw.py %s'%outdir)
#os.chdir(outdir)

dyn=BFGS(atoms,trajectory='out.traj',logfile='out.log',maxstep=0.1)
dyn.run(fmax=0.05)
calc.write('out.gpw')

#os.chdir(home)
write('Relaxed.traj',atoms)
