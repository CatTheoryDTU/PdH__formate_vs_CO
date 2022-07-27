import sys,os
import time
import numpy as np
sys.path.append("/home/cat/geokast/scripts")
from glob import glob
from sjm_neb_tools import save_old_output,set_up_images,run_endstates,restart_NEB,look_for_previous_run,adapt_springs
#from sco2_tools.sco2_neb_tools import neb_calculator
from sco2_tools.requeue import ReQueue
from ase.parallel import world, parprint, paropen
from ase.optimize import BFGS,FIRE
from ase.visualize import view
from gpaw.utilities import h2gpts
from gpaw import GPAW, FermiDirac
from ase.io import write, read
from ase.units import Pascal, m, Bohr
# from ase.neb import NEB
from ase.neb import NEB

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

# SJM parameters
dire = os.getcwd().split('/')[-1]
pot = float(dire.split('_')[-1])
ne = 0.0
dpot = 0.01


# The calculator
method = 'sjm'
def calculator(txtfile,image,potential,xc=None,kpts=None,gpts=None,H2O_layer=None,epsinf=None,dl=None,method='sjm'):
    cell_len = np.linalg.norm(image.get_cell(),axis=1)
    if method != None:
        method = method.lower()
    # Obviously this calculator should be adapted
    if method == 'sjm':
        calc = SJM(poissonsolver={'dipolelayer': 'xy'},
                   doublelayer={'upper_limit': 26},
                   symmetry={'point_group': False},
                   gpts=h2gpts(0.18,image.get_cell(),idiv=8),
                   kpts=(np.around(35/cell_len[0],0),np.around(35/cell_len[1],0),1),
                   xc='BEEF-vdW',
                   spinpol=False,
                   potential=potential,
                   dpot=dpot,
                   #ne=ne,
                   txt=txtfile,
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
    elif method == 'gpaw':
        calc =GPAW(
                   poissonsolver={'dipolelayer': 'xy'},
                   symmetry={'point_group': False},
                   gpts=h2gpts(0.18,image.get_cell(),idiv=8),
                   kpts=(np.around(35/cell_len[0],0),np.around(35/cell_len[1],0),1),
                   xc='BEEF-vdW',
                   spinpol=False,
                   occupations=FermiDirac(0.1),
                   maxiter=400,
                   txt=txtfile,
                   )
    return calc

home=os.getcwd()
path_info = home.split('/')
requeue=ReQueue(maxtime=49.89,checktime=0.1,log='requeue.log')

# NEB default parameters
nimg = 7            # Number of neb images, only applies for a fresh start
potential = pot#float(path_info[-1].split('_')[1])     # Desired potential
fmax = 0.05    #Force threshold
climb = False   #Climb
#Spring constants if int all springs are equal, if list they are different else they will be set to 0.1, if set to a string the code will try to set them automatically
springs=[0.2] + [0.3] + [0.4]*3 + [0.3] + [0.2]*2#[0.2]+[0.3]*3+[0.2]+[0.1]*4
starting_potential=3.9
automatic_restart=True


initial_file = 'IS_start.traj'#'IS_start.traj'#relaxed_pot_4.40_ne_-2.06702.traj'
final_file = 'FS_start.traj'#'FS_relaxed_pot_2.90_ne_1.4074.traj'
try:
    ne_IS=read(initial_file).calc.results['ne']
except:
    ne_IS=ne#0.61177   # First guess of charge on IS
try:
    ne_FS=read(final_file).calc.results['ne']
except:
    ne_FS=ne#0.98067#1.46962  # First guess of charge on FS

#Copy away the last run
save_old_output('previous_run',substitute_init_file=True)
if world.size > 1: time.sleep(10)

#Get the important parameters from the last run (if present)
fresh_start,nimg_last,relax_end_states,climb=\
        look_for_previous_run(['neb_start.traj'],home,fmax=fmax,
                              endfiles=[initial_file,final_file],
                              potential=potential,check_cell_size=False,check_latconst=False)
if nimg_last: nimg=nimg_last
if glob(home+'/IS_relaxed*'):
    IS_paths=glob(home+'/IS_relaxed*')
    IS_paths.sort(key=os.path.getmtime)
    initial_file=IS_paths[-1]
    if glob(home+'/FS_relaxed*'):
        FS_paths=glob(home+'/FS_relaxed*')
        FS_paths.sort(key=os.path.getmtime)
        final_file=FS_paths[-1]
        relax_end_states = False
    else:
        relax_end_states='FS_only'
elif glob(home+'/FS_relaxed*'):
    FS_paths=glob(home+'/FS_relaxed*')
    FS_paths.sort(key=os.path.getmtime)
    final_file=FS_paths[-1]
    relax_end_states='IS_only'
parprint(relax_end_states)
#relax_end_states = False
#relax_end_states='FS_only'

#If the NEB ist started from a converged potential
if abs(potential - starting_potential) > 0.1:
    climb=True
    automatic_restart=False


#Making sure the job doesn't crash because the springs list is too short
if isinstance(springs,list):
    if nimg + 1 > len(springs):
        springs.extend([0.1]*(nimg+1-len(springs)))



if fresh_start:    restart_file = None
else:    restart_file = 'neb_start.traj'
#restart_file='new_band.traj'

#Read out initial charge guesses or set default, align the endstate geometries, check the lattice constant
images,ne_img = set_up_images(nimg,relax_end_states,restart_file,initial_file=initial_file,final_file=final_file,ne_IS=ne_IS,ne_FS=ne_FS,check_cell_size=False,check_latconst=False)

#Automatically set the springs based on charge transfer, manual setting works better normally
if not isinstance(springs,(float,list)):
    springs=adapt_springs(images,springs)

#Some debugging
if world.size == 1:
    print('fresh_start:',fresh_start,
          'nimg: ', nimg,
          'relax_endstates: ', relax_end_states,
          'climb: ', climb)
#Relax endstates, will be passes if relax_end_states is False
else:
    images = run_endstates(images,calculator,ne_img,relax_end_states,potential=potential,fmax=fmax)

#Restart the job if endstates were relaxed, because requeue can not handle more  than one optimization
    if relax_end_states and automatic_restart:
        if not fresh_start:
            write('neb_start.traj',images)
            restart_NEB(sys.argv,debug=False)
        else:
            neb = NEB(images)
            neb.interpolate(mic=True)
            write('neb_start.traj',images)
            restart_NEB(sys.argv,debug=False)
        sys.exit()

# Combine NEB images with their respective calculators
for i in range(1, nimg + 1):
        images[i].set_calculator(calculator(f'image_{i}.txt',images[i],potential=potential,epsinf=epsinf,method=method))
        images[i].calc.ne = ne_img[i]

# Run the NEB
neb = NEB(images, dynamic_relaxation=True,climb=climb,fmax=fmax*2,scale_fmax=0,k=springs)

if not restart_file:
    neb.interpolate(method='idpp',mic=True)

if world.size == 1:
        write('neb_for_view.traj',images)
        sys.exit()

if not climb:
    qn = BFGS(neb, logfile='neb.log', trajectory='neb.traj')
    if automatic_restart:
        status=requeue.run(qn.run,fmax=fmax*2)
        if status != 'time_elapsed':
            time.sleep(60)
            write('neb_GC_final_not_climbed.traj', read(f'neb.traj@-{nimg+2}:'))
        restart_NEB(sys.argv,debug=False)
        sys.exit()
    else:
        qn.run(fmax=fmax*2)
        write('neb_GC_final_not_climbed.traj', images)

neb = NEB(images, dynamic_relaxation=True,climb=True,fmax=fmax,scale_fmax=2,k=springs)
qn = FIRE(neb, logfile='neb.log', trajectory='neb.traj')
requeue=ReQueue(maxtime=49.89,checktime=0.1,log='requeue_climb.log')

if automatic_restart:
        status=requeue.run(qn.run,fmax=fmax)
        if status == 'time_elapsed':
            restart_NEB(sys.argv,debug=False)
        else:
            #write('neb_GC_final_climbed.traj', images)
            time.sleep(60)
            write('neb_GC_final_climbed.traj', read(f'neb.traj@-{nimg+2}:'))
        sys.exit()
else:
        qn.run(fmax=fmax)
        write('neb_GC_final_climbed.traj', images)

