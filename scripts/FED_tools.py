import numpy as np
import math
import sys,os
import pickle as pckl
from matplotlib import pyplot as plt
from general_tools import get_reference_vibrational_contribution
from ase  import units
#from mkm_creator import *

plt.rcParams["figure.figsize"] = (10,5)
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 16
plt.rc('axes', labelsize=28)    # fontsize of the x and y labels
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
plt.rcParams['figure.figsize'] = (10,5)
markersize=10


def  add_CHE_and_energy_vs_RHE(ads_and_electron,facet,show_reaction_steps=True,
        add_field=False,PZC=None,Capacitance=None,pH=None,T=298.,V0SHE=4.4,add_CO_to_C1s=True,add_CO_to_H=True):

    #if 'G_vs_pot_%s'%facet not in ads_and_electron['CO'].keys():
    add_vibrational_contribution(ads_and_electron,facet,barrier=True)

    for entype in ['E','G']:
        enstring=entype+'_vs_pot_%s'%facet
        for ads in ads_and_electron.keys():
            if (enstring not in ads_and_electron[ads].keys() and
                ads[-2:] != '_g'): continue
            if 'nHe' not in ads_and_electron[ads].keys(): continue
            nhe = ads_and_electron[ads]['nHe']

            #Gas  phase species:
            if ads[-2:] == '_g':
                print(ads, ads_and_electron[ads].keys())
                if entype in ads_and_electron[ads]:
                    ads_and_electron[ads][entype+'_vs_RHE']=np.array([nhe,ads_and_electron[ads][entype]])
                continue

            e_0V = np.array(ads_and_electron[ads][enstring])@[-np.log(10)*units.kB*T*pH+V0SHE,1]
            ads_and_electron[ads][entype+'_vs_RHE_%s'%facet]=\
                np.array([ads_and_electron[ads][enstring][0]+nhe,e_0V])

            #For C1 path an extra CO is added in the energy for consistency with C2 pathway
            #TODO: Should maybe be somewhere else
            if add_CO_to_C1s:
             if  ads in ['CHO','COH','CHOH','CH','CH2','CH3','C','CHOH','CH2OH']:
                #TODO: DO not remember what the following two lines did
                #if show_reaction_steps:
                #    e+=RHE_potential+ads_and_electron['CO'][energy_type+'_vs_RHE']
                #else:
                ads_and_electron[ads][entype+'_vs_RHE_%s'%facet][1]+=ads_and_electron['CO'][entype+'_vs_RHE_%s'%facet][1]
            if add_CO_to_H and ads == 'H':
                ads_and_electron[ads][entype+'_vs_RHE_%s'%facet][1]+=ads_and_electron['CO'][entype+'_vs_RHE_%s'%facet][1]
            #For O and OH ethylene is added for consistency with C2 pathway
            elif ads  in ['O','OH']:
                ads_and_electron[ads][entype+'_vs_RHE_%s'%facet][1]+=ads_and_electron['C2H4_g'][entype]

def plot_FED_with_barrier(alldata,facets,included_steps,potentials=[],pH=[13],ylim=[-1.2,0.6],view=False, annotate_intermediates=True,energy_type='G',proton_donor='base',V0_SHE=4.4,normalize_to_IS=False,figsize=(15,5),title=None,annotate_reaction_conditions=True,
        check_barriers_with_line=False,outdir='results',outformat='png',labelsize=None,colors=['k','b','r','g','y'],
        linestyles=['-','--','-.',':'],return_plt=False):
   # plt=plt
#    plt.plot([-0.25,0.25],[0,0],'-k')
    plt.rcParams["figure.figsize"] = figsize

    if not potentials:
        print('No potentials where given for FED with barriers')
        return

    if isinstance(included_steps[0],str):
        included_steps=[included_steps]

    if isinstance(facets,str):
        facets=[facets]

    if len(facets) > 1:
        # If more than one facet is included the facets have different colors
        # and the mechanism different linestyles
         colors=list(reversed(colors[:len(facets)]))
         linestyles=list(reversed(linestyles[:len(included_steps)]))
         print('Several facets have been given, they will have varying '
               'colors and the mechanism will change linestyle')

    else:
        # If only one facet is included the mechanism have different colors
        # and the potentials different linestyles
         colors=list(reversed(colors[:len(included_steps)]))
         #linestyles=list(linestyles[0])*len(included_steps)
         linestyles=list(linestyles)*len(included_steps)
         if len(list(linestyles)) >= len(potentials):
             linestyles=list(linestyles)*len(included_steps)
         print('Only one facet has been given the mechanisms will have varying '
               'colors and the potentials will change linestyle')


    all_ens=[]
    for ifac,facet in enumerate(facets):
     #print(len(facets),ifac,colors)
     enstring='%s_vs_pot_%s'%(energy_type,facet)
     barenstring='%s_ddag_vs_pot_%s'%(energy_type,facet)
     for imech,single_mech_steps_in in enumerate(included_steps):

        #If facets are compared color by facet otherwise by mechanism
        if len(facets) > 1:
            color = colors[ifac]
            print(color)

        else:
            color = colors[imech]

        single_mech_steps=single_mech_steps_in
        added_ads=[None]*len(single_mech_steps_in)
        #If C2  is plotted in the same plot as C1 intermediates the name is given
        # with a plus and will be split here
        if np.any([['+' in i for i in single_mech_steps_in]]):
            single_mech_steps = [i.split('+')[0] for i in single_mech_steps_in]
            for istep,step in enumerate(single_mech_steps_in):
                if step.split('+')[-1] != single_mech_steps_in[istep]:
                    added_ads[istep] = step.split('+')[-1]
                else:
                    added_ads.append(None)

        #Recognize if adsorbate is given with a 2 infront e.g. 2CO
        if np.any([[i[0] == '2' for i in single_mech_steps_in if len(i)]]):
            single_mech_steps = [i.lstrip('2') for i in single_mech_steps]
            for istep,step in enumerate(single_mech_steps_in):
                if not len(step): continue
                if step[0] == '2':#single_mech_steps[istep]:
                    added_ads[istep] = step.lstrip('2')#split('+')
                else:
                    added_ads.append(None)

        if enstring not in alldata[single_mech_steps[0]].keys():
            add_CHE_and_energy_vs_RHE(alldata,facet)
            #add_vibrational_contribution(alldata,facet,barrier=[])


        for ph in pH:
         if isinstance(proton_donor,dict):
            pdonor=proton_donor[ph]
         else:
             pdonor=proton_donor

         for ipot,potential in enumerate(potentials):
            rhe_pot=potential-(V0_SHE-0.059*ph)
            IS_normalization=0
            for istep,current_intermediate in enumerate(single_mech_steps):
#                current_intermediate=single_mech_steps[istep]
                if istep < len(single_mech_steps)-1:
                    next_intermediate,inext_step=single_mech_steps[istep+1],istep+1
                    if next_intermediate == '':
                        next_intermediate=single_mech_steps[istep+2]
                        inext_step=istep+2
                if current_intermediate not in alldata.keys():
                    print('Could not find the intermediate ',current_intermediate)
                    continue
                if enstring not in alldata[current_intermediate]:
                    print('Intermdiate %s does not seem to have an energy with name %s'%(current_intermediate,enstring))
                    continue
                En_IS=np.poly1d(alldata[current_intermediate][enstring])(potential)
                #Add CHE to endstate
                En_IS+=alldata[single_mech_steps[istep]]['nHe']*rhe_pot

                if added_ads[istep] is not None:
                    En_IS+=np.poly1d(alldata[added_ads[istep]][enstring])(potential)+\
                            alldata[added_ads[istep]]['nHe']*rhe_pot

                if normalize_to_IS and istep==0:
                    IS_normalization=-En_IS

                if current_intermediate == 'CO' and 'OCCO' in single_mech_steps:
                    if barenstring not in alldata['CO']:
                        alldata['CO'][barenstring] = {}

                    if barenstring in alldata['CO']:
                        if 'OCCO' not in alldata['CO'][barenstring].keys():
                            alldata['CO'][barenstring]['OCCO']= alldata['COCO'][barenstring]['OCCO']



                if annotate_intermediates:
                    stepout=current_intermediate
                    if added_ads[istep] is not None:
                        stepout+='+'+added_ads[istep]
                        if current_intermediate == added_ads[istep]:
                            stepout = '2'+current_intermediate

                    # Subscript numbers - TODO
#                    inum=[]
#                    for i in stepout:
#                        if i.isnumeric(): inum.append(i)
#                    for i in reversed(setpout)
#                    for i in reversed(inum):
#                        stepout

                    if labelsize is None:
                        plt.annotate(stepout,xy=(inext_step,ylim[1]-0.1*(imech+1)),color=color).draggable()
                    else:
                        plt.annotate(stepout,xy=(inext_step,ylim[1]-0.1*(imech+1)),color=color,fontsize=labelsize).draggable()

                #if istep == len(single_mech_steps)-1:
                #    plt.plot([istep+0.75,istep+1.25],[En_IS+IS_normalization,En_IS+IS_normalization],'-'+colors[imech])
                #    continue

                if istep < len(single_mech_steps)-1:
                    if next_intermediate == '':
                        print('Skipping the step %i'%inext_step)
                    elif next_intermediate not in alldata.keys():
                        print('Could not find the intermediate for FS:',next_intermediate)
                        continue
                    else:
                        if enstring not in alldata[next_intermediate]:
                            print('Intermdiate %s_%s does not seem to have an energy'%(next_intermediate,facet))
                            continue

                    print(next_intermediate,current_intermediate)
                    En_FS=np.poly1d(alldata[next_intermediate][enstring])(potential)+\
                            alldata[next_intermediate]['nHe']*rhe_pot

                    if added_ads[inext_step] is not None:
                        En_FS+=np.poly1d(alldata[added_ads[inext_step]][enstring])(potential)+\
                                        alldata[added_ads[inext_step]]['nHe']*rhe_pot


                    #If no  barrier at all has  been calculated from the current intermediate
                    #Draw a straight line to the next intermediate
                    Eddag=None
                    #if step == 'CO2':
                    #    print(single_mech_steps)
                    #    print(alldata['CO2'][barenstring],single_mech_steps[istep+1])
                    #    print(barenstring)
                    #    das
                    if barenstring not in alldata[current_intermediate].keys():
                        plt.arrow(istep+1.25,En_IS,0.5,En_FS-En_IS,head_width=0.0,length_includes_head=True, linewidth=0.1,color='r')
                        pass
                        #continue

                    #If the barrier between IS and FS has been calculated
                    elif next_intermediate in alldata[current_intermediate][barenstring]:
                     toads=next_intermediate
                     if pdonor in alldata[current_intermediate][barenstring][toads]:
                        Eddag=np.poly1d(alldata[current_intermediate][barenstring][toads][pdonor])(potential)
                    #For a chemical step
                    elif next_intermediate == 'H':
                        #Eddag=np.poly1d(alldata['clean'][barenstring]['H'][pdonor])(potential)
                        Eddag=np.poly1d(alldata['clean'][barenstring]['H']['base'])(potential)

                    print(facet,Eddag,single_mech_steps[istep])

                     #Add CHE to barriers (alkaline has CHE like IS, acid has CHE like FS)
                    if Eddag is not None:
                        Eddag+=alldata[current_intermediate]['nHe']*rhe_pot
                        if added_ads[istep] is not None and added_ads[inext_step] is not None:
                                Eddag+=np.poly1d(alldata[added_ads[istep]][enstring])(potential)


                        if pdonor == 'acid':
                                Eddag+=rhe_pot
                        #        plt.plot([istep+1.35,istep+1.65],[Eddag,Eddag],'--b')
                        #else:
                        #        plt.plot([istep+1.35,istep+1.65],[Eddag,Eddag],'--'+colors[imech])

                        #Connect the states by lines (TODO: Maybe polynoms?)

                #plot everything
                #Plot IS
                if len(facets) > 1:
                    linestyle=linestyles[imech]
                else:
                    linestyle=linestyles[ipot]
                print(current_intermediate,potential,En_IS+IS_normalization)
                plt.plot([istep+0.75,istep+1.25],
                        [En_IS+IS_normalization,En_IS+IS_normalization],linestyle=linestyle,
                        color=color,linewidth=4)
                all_ens.append(En_IS+IS_normalization)
                if istep == len(single_mech_steps)-1:
                    continue

                #If the barrer to the FS has not been calculated
                #Draw a straight line to the next intermediate
                #print(alldata[step].keys())
                #das
                all_ens.append(En_FS+IS_normalization)
                arrow_xlen=0.5+(inext_step-istep-1)
                straight_linestyle=linestyle
                if linestyle=='--': straight_linestyle=':' #Hack because -- looks like a solid line
                if barenstring not in alldata[current_intermediate]:# and single_mech_steps[istep+1] != 'H':
                        print(alldata[current_intermediate].keys())
                        print(f'Barrier for step {current_intermediate} not found')
                        plt.arrow(istep+1.25,En_IS+IS_normalization,arrow_xlen,
                                En_FS-En_IS,head_width=0.0,length_includes_head=True,
                                linewidth=1.5,color=color,linestyle=straight_linestyle)

                elif (next_intermediate not in alldata[current_intermediate][barenstring] and
                        next_intermediate != 'H'):
                        plt.arrow(istep+1.25,En_IS+IS_normalization,arrow_xlen,
                                En_FS-En_IS,head_width=0.0,length_includes_head=True,
                                linewidth=1.5,color=color,linestyle=straight_linestyle)

                #Draw parabolic barrier
                elif Eddag is not None:
                    parpts=np.array([[istep+1.25,En_IS],
                        [istep+1+(inext_step-istep)/2+0.25*((En_FS-En_IS)/3.),Eddag],
                        [inext_step+0.75,En_FS]])
                    from general_tools import quad_fun
                    from scipy.optimize import curve_fit
                    coeff,dummy=curve_fit(quad_fun,parpts[:,0],parpts[:,1])
                    fitpts=np.linspace(istep+1.25,inext_step+0.75,30)#,include_endpoints=True)
                    fit=np.poly1d(coeff)(fitpts)+IS_normalization
                    #print(fit)
                    if pdonor == 'acid':
                            plt.plot(fitpts,fit,'--b')
#                            plt.plot([istep+1.35,istep+1.65],[Eddag+IS_normalization,Eddag+IS_normalization],'--b')
                    else:
                            #plt.plot(fitpts,fit,linestyles[imech],color=color)
                            plt.plot(fitpts,fit,linestyle,color=color)
                            if check_barriers_with_line:
                               plt.plot([inext_step+0.35,inext_step+0.65],[Eddag+IS_normalization,Eddag+IS_normalization],'--'+colors[imech])

                    all_ens.append(Eddag+IS_normalization)
#                    plt.arrow(istep+1.25,En_IS+IS_normalization,0.1,Eddag-En_IS,head_width=0.0,
#                              length_includes_head=True, linewidth=0.1,color=colors[imech])
#                    plt.arrow(istep+1.65,Eddag+IS_normalization,0.1,En_FS-Eddag,head_width=0.0,
#                            length_includes_head=True, linewidth=0.1,color= colors[imech])



    plt.ylim(ylim)
    plt.xticks([])
    if proton_donor=='base':
        donor='H$_2$O'
    elif proton_donor == 'acid':
        donor='H$_3$O$^+$'
    elif proton_donor == 'chemical':
        donor='chemical'

    elif isinstance(proton_donor,dict):
        donor=list(proton_donor.items())
#    print(','.join(['-'.join(i) for i in included_steps]))
    if title is None:
        plt.title(#','.join(['-'.join(i) for i in included_steps])+
            ', facet: '+facet+
            ', WF='+'-'.join([str(np.around(i,2)) for i in potentials])+
            ', pH='+'-'.join([str(np.around(i,1)) for i in pH])+
            ', proton donor: %s '%donor,fontsize=15)
        #'-'.join(included_steps))+
    else:
        plt.title(title)

    plt.ylabel('$\Delta$G$^\phi$ / eV')

    if annotate_reaction_conditions:
        sheout=','.join([str(np.around(i-V0_SHE,2)) for i in potentials])+'V$_{\mathrm{SHE}}$'
        phout = ','.join([str(np.around(i,1)) for i in pH])
        plt.annotate(f"Cu({facet})\n{sheout}\npH={phout}",(0.8,min(all_ens)),ha='left',va='bottom',fontsize=23).draggable()

    plt.tight_layout()
    if view:
        if return_plt:
            return plt
        else:
            plt.show()
            return
    if potentials:
        print(included_steps)
        if isinstance(included_steps[0],list):
            stepsout=['-'.join(steps) for steps in included_steps]
        else:
            stepsout=included_steps
        print(stepsout)
        plt.savefig(outdir+'/FED_w_barrier_'+
                '_'.join(stepsout)+
                '_pot_'+'-'.join([str(np.around(i,2)) for i in potentials])+
                '_pH_'+'-'.join([str(np.around(i,2)) for i in pH])+
                '.'+outformat,transparent=True)
    else:
        plt.savefig(outdir+'/FED_w_barrier_'+'-'.join(included_steps)+'.'+outformat,transparent=True)
    plt.close()

#main()
def read_beta_from_catmap_input(line):
            beta=None
            if 'beta' in line:
                beta_string=line.split('beta=(')[1]
                beta_string=beta_string.split(')')[0]
                beta=[float(i.replace(',','')) for i in beta_string.split()]
            return beta

