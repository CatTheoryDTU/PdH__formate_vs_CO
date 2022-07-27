from catmap import ReactionModel
from matplotlib import pyplot as plt
import matplotlib
import sys
import os



#Perform a run for solving the microkinetic model
mkm_file = 'model.mkm'
model = ReactionModel(setup_file=mkm_file)
model.output_variables += ['rate', 'production_rate','coverage','free_energy']
model.run()

exit()

output = 'output/'
os.system('mkdir -p ' + output )
#Analysis
from catmap import analyze
vm = analyze.VectorMap(model)
vm.plot_variable = 'production_rate' #tell the model which output to plot
vm.log_scale = True #rates should be plotted on a log-scale
vm.min = 1e-5 #minimum rate to plot
vm.max = 1e15 #maximum rate to plot
vm.threshold = 1e-25 #anything below this is considered to be 0
vm.descriptor_labels = ['Potential vs SHE', 'TOF [site$^{-1}$s$^{-1}$]']
#vm.subplots_adjust_kwargs = {'left':0.2,'right':0.8,'bottom':0.15}
vm.plot(save=output+'production_rate.pdf')

vm.plot_variable = 'rate' #tell the model which output to plot
vm.log_scale = True #rates should be plotted on a log-scale
vm.min = 1e-10 #minimum rate to plot
vm.max = 1e5 #maximum rate to plot
vm.plot(save=output+'rate.pdf')

vm.plot_variable = 'coverage'
vm.log_scale = True
vm.min = 1e-20
vm.max = 1

vm.descriptor_labels = ['$\Delta$E$_{b}$(CO) [eV]', '$\Delta$E$_{b}$(CO2) [eV]']
vm.plot(save=output+'coverage.pdf')

ma = analyze.MechanismAnalysis(model)
ma.energy_type = 'free_energy' #can also be free_energy/potential_energy
ma.pressure_correction = False #assume all pressures are 1 bar (so that energies are the same as from DFT)
ma.include_labels = True
for i in [-0.4]:
    fig = ma.plot(plot_variants=[i], save=output+f'FED_{i+0.4:.1f}VRHE.png')

#sa = analyze.ScalingAnalysis(model)
#sa.plot(save='output/scaling.pdf')
