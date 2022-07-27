import pickle as pkl
from general_tools import get_reference_energies,lin_fun
from matplotlib import pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

a=pkl.load(open('/Users/geokast/Papers/PdH_CO2/results/parsed_all_data.pkl','rb'))
ETS=a['H-CO2']['111']['barrier']['E']
Eslab=a['clean_slab']['111']['slab']['E_v_pot']
pltdat=[]
for i in a['H-CO2']['111']['barrier']['E']:
  gasref=get_reference_energies('COO',code='GPAW',references={'C':'CO2'})
  slabE=Eslab[0]*i[0]**2+Eslab[1]*i[0]+Eslab[2]
  pltdat.append([i[0],i[1]-slabE-gasref])
pltdat=np.array(pltdat)
plt.plot(pltdat[:,0],pltdat[:,1],'ok')

coeff,d=curve_fit(lin_fun,pltdat[:,0],pltdat[:,1])
plt.plot(pltdat[:,0],coeff[0]*pltdat[:,0]+coeff[1],label='Slope: %1.2feV/V'%coeff[0])
plt.legend()

plt.show()
