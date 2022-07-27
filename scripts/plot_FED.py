import pickle as pkl
import sys,os
sys.path.append('scripts')
from FED_tools import plot_FED_with_barrier
import numpy as np

SHE_potentials=[-0.4,-0.9]
pH=6.8
pka_formic=3.75
facet='111'
SHE_abs=4.4

ads_and_electron=pkl.load(open('results/parsed_most_stable_data.pkl','rb'))
ads_and_electron['HCOO$^-_{(aq)}$']['G_vs_pot_%s'%facet]=ads_and_electron['HCOOH$_{(aq)}$']['G_vs_pot_%s'%facet]
ads_and_electron['HCOO$^-_{(aq)}$']['G_vs_pot_%s'%facet][1]-=0.059*(pH-pka_formic)

plt=plot_FED_with_barrier(ads_and_electron,facets=[facet],
            included_steps=[
            ['CO$_{2(g)}$','CO2','COOH','CO','CO$_{(g)}$'],
            ['CO$_{2(g)}$','','HCOO$^-_{(aq)}$+clean_minus_H','HCOO$^-_{(aq)}$'],
            ['CO$_{2(g)}$'], #This line is just for making CO2_g black
            ],
            potentials=[i+SHE_abs for i in SHE_potentials],
            pH=[pH],
            ylim=[-1.1,1.5],
            view=True,
            annotate_intermediates=False,
            annotate_reaction_conditions=False,
            energy_type='G',
            figsize=(10,7),
            check_barriers_with_line=False,
            outdir='results',outformat='png',labelsize=None,colors=['k','r','b','r','g','y'],
            linestyles=['-','--'],
            return_plt=True,
            title='')

plt.yticks(fontsize=28)

plt.plot(np.nan,np.nan,'-k',label='0V$_{\mathrm{RHE}}$',linewidth=3)
plt.plot(np.nan,np.nan,'--k',label='-0.5V$_{\mathrm{RHE}}$',linewidth=3)
plt.annotate('CO$_{2(g)}$',(1,0.05),ha='center',fontsize=20).draggable()

#Write out the intermediates of the CO path ontop
COpath_y=1.47
plt.annotate('CO pathway:',(0.6,COpath_y),ha='left',va='top',fontsize=20,color='b').draggable()
plt.annotate('*CO$_{2}$',(2,COpath_y),ha='center',va='top',fontsize=20,color='b').draggable()
plt.annotate('*COOH',(3,COpath_y),ha='center',va='top',fontsize=20,color='b').draggable()
plt.annotate('*CO',(4,COpath_y),ha='center',va='top',fontsize=20,color='b').draggable()
plt.annotate('CO$_{(g)}$',(5,COpath_y),ha='center',va='top',fontsize=20,color='b').draggable()

#Write out the intermediates of the formate path ontop
formpath_y=-1.1
plt.annotate('Formate\npathway:',(0.6,formpath_y),ha='left',va='bottom',fontsize=20,color='r').draggable()
plt.annotate('*H-CO$_{2}$+\nPd$_x$H$_{(x-1)}$',(2,formpath_y),ha='center',va='bottom',fontsize=20,color='r').draggable()
plt.annotate('HCOO$^-_{(\mathrm{aq})}$+\nPd$_x$H$_{(x-1)}$',(3,formpath_y),ha='center',va='bottom',fontsize=20,color='r').draggable()
plt.annotate('HCOO$^-_{(\mathrm{aq})}$+\nPd$_x$H$_{(x)}$',(4,formpath_y),ha='center',va='bottom',fontsize=20,color='r').draggable()

# Write the values of the slope explicitly
TSslope=ads_and_electron['CO$_{2(g)}$']['G_ddag_vs_pot_%s'%facet]['HCOO$^-_{(aq)}$']['base'][0]
plt.arrow(1.9,0.92,0,0.1,length_includes_head=True,color='r',head_width=0.05,width=0.02,head_length=0.03)
plt.annotate(r'%1.2f $\frac{eV}{V}$'%(TSslope),(2.15,1.02),va='top',ha='center',color='r')

CO2slope=ads_and_electron['CO2']['G_vs_pot_111'][0]
plt.arrow(2.2,0.42,0,0.14,length_includes_head=True,color='b',head_width=0.05,width=0.02,head_length=0.03)
plt.annotate(r'%1.2f $\frac{eV}{V}$'%CO2slope,(1.95,0.55),va='top',ha='center',color='b').draggable()

COOHslope=ads_and_electron['COOH']['G_vs_pot_111'][0]+1 #The +1 is the CHE contribution
plt.arrow(3.0,0.66,0,0.5,length_includes_head=True,color='b',head_width=0.05,width=0.02,head_length=0.03)
plt.annotate(r'%1.2f $\frac{eV}{V}$'%COOHslope,(3.0,0.88),va='center',ha='center',color='b',bbox={'fc':'w','ec':'w'}).draggable()

plt.legend(loc='upper right',facecolor='w',edgecolor='w',bbox_to_anchor=(0.22,0.94,0,0),fontsize=16)
plt.show()
