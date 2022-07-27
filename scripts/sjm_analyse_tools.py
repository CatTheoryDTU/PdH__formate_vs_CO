from ase.io import read,write
from scipy.optimize import curve_fit,OptimizeWarning

import matplotlib
import pickle
from matplotlib import pyplot as plt
from matplotlib.pyplot import cm
import os,sys
import numpy as np
from ase import units
import warnings
from general_tools import get_reference_energies,lin_fun,quad_fun,get_reference_vibrational_contribution

#matplotlib.rc('text', usetex=True)
#matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
plt.rcParams["font.family"] = "Times New Roman"
plt.rc('axes', labelsize=28)    # fontsize of the x and y labels
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
plt.rcParams['figure.figsize'] = (7,5)
markersize=10

from matplotlib.ticker import FormatStrFormatter

def get_E_vs_pot(alldata,ads,potentials,enkey,deeperkey=None):
        E_v_pot=[]
        for potential in np.around(potentials,2):
            if deeperkey is not None:
              if isinstance(deeperkey,str):
                if str(potential) in alldata[ads][enkey][deeperkey].keys():
                    E_v_pot.append([potential,alldata[ads][enkey][deeperkey][str(potential)]])
              #if the dict goes deeper (e.g. for barriers)
              elif isinstance(deeperkey,list):
                    fulldepth = len(deeperkey)
                    en=alldata[ads][enkey]
                    for depth in range(fulldepth):
                        en=en[deeperkey[depth]]
                    if str(potential) in en.keys():
                        E_v_pot.append([potential,en[str(potential)]])
                    elif '%1.2f'%potential in en.keys():
                        E_v_pot.append([potential,en['%1.2f'%potential]])

            else:
                if potential in alldata[ads][enkey].keys():
                    E_v_pot.append([potential,alldata[ads][enkey][potential]])
        return np.array(E_v_pot)

def fit_potential_response(alldata,facet,potentials,include_barriers=True,plot=True,plotoutname='E_',specific_ads=None,quad_fit=False,not_enough_potentials=None):
    #TURN THIS THE FOLLOWING OFF FOR DEBUGGING!
    warnings.simplefilter('ignore',category=OptimizeWarning)
    if not_enough_potentials is None:
        not_enough_potentials=[]
    for iads,ads in enumerate(alldata.keys()):
     #If only a single adsorbate should be fitted
     if specific_ads:
         if isinstance(specific_ads,str):
             specific_ads=[specific_ads]
         if ads not in specific_ads: continue
     if ads[-2:] == '_g': continue
     if 'E_%s'%facet in alldata[ads].keys():
        E_v_pot = get_E_vs_pot(alldata,ads,potentials,'E_%s'%facet)
        if quad_fit:
            if len(E_v_pot) > 2:
                coeff,dummy = curve_fit(quad_fun,E_v_pot[:,0],E_v_pot[:,1])
                alldata[ads]['E_vs_pot_%s'%facet]=coeff
                if len(E_v_pot) == 3:
                    not_enough_potentials.append(ads)
        else:
            if len(E_v_pot) > 1:
                coeff,dummy = curve_fit(lin_fun,E_v_pot[:,0],E_v_pot[:,1])
                alldata[ads]['E_vs_pot_%s'%facet]=coeff
            if len(E_v_pot) == 2:
                not_enough_potentials.append(ads)

     if 'E_C_%s'%facet in alldata[ads].keys():
        E_v_pot = get_E_vs_pot(alldata,ads,potentials,'E_%s'%facet)
        if len(E_v_pot) > 1:
            coeff,onset = curve_fit(lin_fun,E_v_pot[:,0],E_v_pot[:,1])
            alldata[ads]['E_C_vs_pot_%s'%facet]=coeff
        if len(E_v_pot) == 2:
                not_enough_potentials.append(ads)

    if include_barriers:
        alldata = fit_barriers_vs_pot(alldata,facet,potentials,plot=plot)

    if plot:
        plot_E_vs_pot(alldata,'E',facet,potentials,quad_fit=quad_fit)
        plot_betas(alldata,facet)

    #warnings.simplefilter('default')#,category='OptimizeWarning')
    return alldata,not_enough_potentials


def fit_barriers_vs_pot(alldata,facet,potentials,plot=False):
    for iads,ads in enumerate(alldata.keys()):
     if 'E_ddag_%s'%facet not in alldata[ads].keys():
         continue
     for toads in alldata[ads]['E_ddag_%s'%facet].keys():

      for pH in ['base','acid','chemical']:
        for quad_fit in [False,True]:

          if pH not in alldata[ads]['E_ddag_%s'%facet][toads].keys():
              continue
          en = alldata[ads]['E_ddag_%s'%facet][toads][pH]
          if len(en.keys()) > 2:
             E_v_pot = get_E_vs_pot(alldata,ads,potentials,'E_ddag_%s'%facet,deeperkey=[toads,pH])
             #print(ads,pH,en,E_v_pot)
             if quad_fit:
                 if len(E_v_pot) < 3:
                     print("Quadratic fit for %s to %s failed due to lack of  potentials"%(ads,toads))
                     continue
                 coeff,dummy = curve_fit(quad_fun,E_v_pot[:,0],E_v_pot[:,1])
             else:
                 coeff,dummy = curve_fit(lin_fun,E_v_pot[:,0],E_v_pot[:,1])
          elif quad_fit: continue
          elif len(en.keys()) > 1:
             E_v_pot = get_E_vs_pot(alldata,ads,potentials,'E_ddag_%s'%facet,deeperkey=[toads,pH])
             if len(E_v_pot) > 1:
                 coeff,dummy = curve_fit(lin_fun,E_v_pot[:,0],E_v_pot[:,1])
             else:
                 continue
          else:
              continue

          fitname='E_ddag_vs_pot'
          if quad_fit:
              fitname+='_quad'#E_ddag_vs_pot_quad'
          fitname+='_%s'%facet
          if fitname not in alldata[ads].keys():
                 alldata[ads][fitname]={}
          if toads not in alldata[ads][fitname].keys():
                 alldata[ads][fitname][toads]={}
          alldata[ads][fitname][toads][pH]=coeff

          #Add barrier referenced to initial state to dictionary
          if 'E_vs_pot_%s'%facet not in alldata[ads]:
              print('Thermodynamics of IS for %s to %s on facet %s are missing'
                      %(ads,toads,facet))
              continue

          E_v_pot2=np.array(E_v_pot)
          E_v_pot2[:,1]=1
          E_v_pot[:,1]-=E_v_pot2@np.array(alldata[ads]['E_vs_pot_%s'%facet])

          if 'E_ddag_rel_%s'%facet not in alldata[ads].keys():
                 alldata[ads]['E_ddag_rel_%s'%facet]={}
          if toads not in alldata[ads]['E_ddag_rel_%s'%facet].keys():
                 alldata[ads]['E_ddag_rel_%s'%facet][toads]={}
          alldata[ads]['E_ddag_rel_%s'%facet][toads][pH]={}
          for pot_E in E_v_pot:
              alldata[ads]['E_ddag_rel_%s'%facet][toads][pH]['%s'%pot_E[0]]=pot_E[1]

          fitname='E_ddag_rel_vs_pot'
          if quad_fit:
              fitname+='_quad'#E_ddag_vs_pot_quad'
          fitname+='_%s'%facet
          if quad_fit:
            coeff,dummy = curve_fit(quad_fun,E_v_pot[:,0],E_v_pot[:,1])
          else:
              coeff,dummy = curve_fit(lin_fun,E_v_pot[:,0],E_v_pot[:,1])
          if fitname not in alldata[ads].keys():
                 alldata[ads][fitname]={}
          if toads not in alldata[ads][fitname].keys():
                 alldata[ads][fitname][toads]={}
          alldata[ads][fitname][toads][pH]=coeff
    if plot:
        quad_fit=False
        plot_E_vs_pot(alldata,'E_ddag',facet,potentials,True,quad_fit)
        plot_E_vs_pot(alldata,'E_ddag_rel',facet,potentials,True,quad_fit,ylabel='E$_a$ [eV]')
        plot_Eddag_from_specific_IS(alldata,'E_ddag',facet,potentials,True,quad_fit)

    return alldata

    #Plot BEP
def plot_E_vs_pot(alldata,enkey,facet,potentials,deeperkey=None,quad_fit=False,plot_separate=True,ylabel='Formation energy [eV]'):
    counter=0
    for iads,ads in enumerate(alldata.keys()):
     if enkey+'_%s'%facet in alldata[ads].keys():
         counter+=1

    colors=cm.gist_ncar(np.linspace(0,1,counter+1))
    counter=0
    symbols=['+','s','d']
    for iads,ads in enumerate(alldata.keys()):
     if ads[-2:] == '_g': continue
     if enkey+'_vs_pot_%s'%facet in alldata[ads].keys() and ads not in ['clean']:
        if deeperkey:
          counter2=0
          #print(ads,enkey)
          for toads in alldata[ads][enkey+'_vs_pot_%s'%facet].keys():
           for pH in alldata[ads][enkey+'_vs_pot_%s'%facet][toads].keys():
            coeff=alldata[ads][enkey+'_vs_pot_%s'%facet][toads][pH]
            fit = []
            if quad_fit:
                for pot in np.linspace(potentials[0],potentials[-1],10):
                    fit.append([pot,np.array([pot**2,pot,1])@coeff])
            else:
                for pot in np.linspace(potentials[0],potentials[-1],2):
                    fit.append([pot,np.array([pot,1])@coeff])
            fit=np.array(fit)
            plt.plot(fit[:,0],fit[:,1],color=colors[counter%len(colors)])
            E_v_pot=get_E_vs_pot(alldata,ads,potentials,enkey+'_%s'%facet,[toads,pH])
            if len(E_v_pot):
                plt.plot(E_v_pot[:,0],E_v_pot[:,1],symbols[counter2%len(symbols)],label=ads+'-'+toads+', '+pH+r',dE/d$\Phi$=%1.2f'%coeff[0],color=colors[counter%len(colors)])
                #plt.plot(E_v_pot[:,0],E_v_pot[:,1],symbols[counter%len(symbols)],label=ads+'-'+toads+r',dE/d$\Phi$=%1.2f'%coeff[0],color=colors[counter%len(colors)])
            counter2+=1
        else:
            coeff=alldata[ads][enkey+'_vs_pot_%s'%facet]
            fit = np.array([[potentials[0],potentials[0]*coeff[0]+coeff[1]],
                           [potentials[-1],potentials[-1]*coeff[0]+coeff[1]]])
            plt.plot(fit[:,0],fit[:,1],color=colors[counter%len(colors)])
            E_v_pot=get_E_vs_pot(alldata,ads,potentials,enkey+'_%s'%facet)
            plt.plot(E_v_pot[:,0],E_v_pot[:,1],symbols[counter%len(symbols)],label=ads+r',dE/d$\Phi$=%1.2f'%coeff[0],color=colors[counter%len(colors)])
        counter+=1

    plt.xlabel('Work function [eV]')
    plt.ylabel(ylabel)
    #print(enkey,facet,counter)
    plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=4,fontsize=3)
    plt.tight_layout()
    plt.savefig('results/%s_vs_pot_%s.pdf'%(enkey,facet))
    if plot_separate:
        plt.close()


def write_catinput(alldata,facets,potential,outfile,include_barriers=True,
                    catbasein='/Users/geokast/SelectCO2/endstates/tools_for_analysis/catinput_base.txt',
                    catbasefreein='/Users/geokast/SelectCO2/endstates/tools_for_analysis/catinput_base_freeen.txt'):
    out=open(outfile,'w')
    out_G=open('.'.join(outfile.split('.')[:-1])+'_freEn.txt','w')
    basein=open(catbasein,'r').read()
    basefreein=open(catbasefreein,'r').read()
    out.write(basein)
    out_G.write(basefreein)

    if isinstance(facets,str):facets=[facets]
    ##Get the free energy contribution of one water on the clean slab for referencing barriers
    #from ase.thermochemistry import HarmonicThermo
    #from ase.units import invcm

    for facet in facets:
        written_barriers=[]
        for ads in alldata.keys():
            ads_short = ads.lstrip('md-').lstrip('bdo-')
            vib_string=''
            if 'vibs_%s'%facet in alldata[ads].keys():
                vib_string_ads=alldata[ads]['vibs_%s'%facet]
                for vib in alldata[ads]['vibs_%s'%facet]:
                    vib_string+=str(np.around(vib,6))+', '
                vib_string='['+vib_string[:-1]+']'

            if not len(vib_string):
                vib_string = '[]'

            #testprint(ads,facet)
            #try:
            #    testprint(alldata[ads]['E_vs_pot_%s'%facet])
            #except:
            #    testprint(alldata[ads].keys())
            #    sys.exit()
            if 'E_vs_pot_%s'%facet in alldata[ads].keys():

                testprint(ads,facet)
                coeff = tuple(alldata[ads]['E_vs_pot_%s'%facet])

                out.write('Cu\t%s\t%s\t%f\t%s\tbeta=%s\n'%
                          (facet,ads_short,np.array([potential,1])@coeff,vib_string,coeff))
                if 'free_en_corr_%s'%facet in alldata[ads]:
                    out_G.write('Cu\t%s\t%s\t%f\t%s\tbeta=%s\n'%
                          (facet,ads_short,(np.array([potential,1])@coeff)+alldata[ads]['free_en_corr_%s'%facet]
                              ,'[]',coeff))
                else:
                    print(f'Could not write the G of {ads_short} on {facet} in catmap energy file. missing vibs?')


            if 'E_ddag_vs_pot_%s'%facet in alldata[ads].keys():
                for toads in alldata[ads]['E_ddag_vs_pot_%s'%facet].keys():
                  ads_short = ads.lstrip('md-').lstrip('bdo-')
                  toads_short = toads.lstrip('md-').lstrip('bdo-')

                  for pH in ['base','acid','chemical']:
                    if pH not in alldata[ads]['E_ddag_vs_pot_%s'%facet][toads].keys():
                        continue

                    #if pH == 'acid':
                    #    print('catmap in',ads,toads,
                    #            alldata[ads]['E_ddag_vs_pot_%s'%facet][toads][pH],
                    #           alldata[ads]['E_ddag_vs_pot_%s'%facet][toads][pH][0]*4.4+alldata[ads]['E_ddag_vs_pot_%s'%facet][toads][pH][1])
                    outname,written_barriers,TSnamediff=get_unique_TS_name(ads,toads,pH,written_barriers)
                    coeff=tuple(alldata[ads]['E_ddag_vs_pot_%s'%facet][toads][pH])
                    out_string='Cu\t%s\t%s\t%f\t%s\tbeta=%s %s_to_%s_%s %s\n'

                    #Add vibrations to transition states
                    ## If no vibrations are found use the IS vibrations for
                    ## and FS vibrations otherwise
                    no_vibs=True
                    if 'vibs_ddag_%s'%facet in alldata[ads]:
                        if toads in alldata[ads]['vibs_ddag_%s'%facet]:
                          if pH in alldata[ads]['vibs_ddag_%s'%facet][toads]:
                            vib_string=alldata[ads]['vibs_ddag_%s'%facet][toads][pH][1:]
                            img_freq=alldata[ads]['vibs_ddag_%s'%facet][toads][pH][0]
                            no_vibs=False

                    if no_vibs:
                            #Add FS vibs as default for acidic barriers
                            if pH in ['acid','chemical']:
                                substitute_endstate = toads
                            #Add IS vibs as default for acidic barriers
                            else:
                                substitute_endstate = ads

                            if 'vibs_%s'%facet in alldata[substitute_endstate]:
                                vib_string=alldata[substitute_endstate]['vibs_%s'%facet]
                            else:
                                vib_string='[]'
#                            else:
#                                vib_string=alldata[ads]['vibs_%s'%facet]
                            img_freq=0

                    if len(coeff) == 3:
                        out.write(out_string%
                              (facet,outname,np.array([potential**2,potential,1])@coeff,
                              vib_string,coeff,ads_short,toads_short,pH,
                              np.around(img_freq,0)))
                    else:
                        out.write(out_string%
                              (facet,outname,np.array([potential,1])@coeff,
                               vib_string,coeff,ads_short,toads_short,pH,
                               np.around(img_freq,0)))
                        print(ads,out_string)
                    if 'free_en_corr_ddag_%s'%facet in alldata[ads]:
                     print(ads,toads)
                     if toads in alldata[ads]['free_en_corr_ddag_%s'%facet]:
                        if len(coeff) == 3:
                            out_G.write(out_string%
                                  (facet,outname,np.array([potential**2,potential,1])@coeff+alldata[ads]['free_en_corr_ddag_%s'%facet][toads][pH],
                                  '[]',coeff,ads_short,toads_short,pH,
                                  np.around(img_freq,0)))
                        else:
                            out_G.write(out_string%
                                  (facet,outname,np.array([potential,1])@coeff+alldata[ads]['free_en_corr_ddag_%s'%facet][toads][pH],
                                  '[]',coeff,ads_short,toads_short,pH,
                                  np.around(img_freq,0)))


def add_vibrational_free_energy_corrections(alldata,facet,no_water_layer=False,references={'C':'CO'}):
    from ase.thermochemistry import HarmonicThermo
    from ase.units import invcm

    #Get the free energy contribution of one water on the clean slab for referencing barriers
    if not no_water_layer:
        if 'free_en_corr' not in alldata['clean'].keys():
            try:
                vibens=[i*invcm for i in alldata['clean']['vibs_%s'%facet]]
            except:
                vibens=[0]*9#[i*invcm for i in alldata['clean']['vibs_%s'%facet]]
            alldata['clean']['free_en_corr_%s'%facet]=clean_G=HarmonicThermo(vibens).get_helmholtz_energy(298,verbose=False)
        else:
            print('It seems the add_vibrational_free_energy_corrections class is called more than  once')
            return

    #Add vibration correction to adsorbates
    no_freen=[]
    for ads in alldata.keys():
      if ads[-2:] == '_g':
          continue


      if 'E_vs_pot_%s'%facet not in alldata[ads]: continue

      ads_short = ads.lstrip('md-').lstrip('bdo-')
      if ads == 'clean':
          ads_short = ''
      vib_string=''
      vibens=[]
      if 'vibs_%s'%facet in alldata[ads].keys():
          adsvibs=alldata[ads]['vibs_%s'%facet]
          for vib in alldata[ads]['vibs_%s'%facet]:
              vib_string+=str(np.around(vib,4))+', '
              vibens.append(vib*invcm)

          alldata[ads]['free_en_corr_%s'%facet]=HarmonicThermo(vibens).get_helmholtz_energy(298,verbose=False)
          #Subtract gas phase reference free energies
          ads_short_for_vibs=ads_short.replace('2','H').replace('3','HH')
          alldata[ads]['free_en_corr_%s'%facet]-=get_reference_vibrational_contribution(ads_short_for_vibs,references=references)
      else:
          no_freen.append(ads)
          alldata[ads]['free_en_corr_%s'%facet]=0

      alldata[ads]['G_vs_pot_%s'%facet] = alldata[ads]['E_vs_pot_%s'%facet].copy()
      alldata[ads]['G_vs_pot_%s'%facet][1] += alldata[ads]['free_en_corr_%s'%facet]

      #Transition states
      if 'E_ddag_vs_pot_%s'%facet in alldata[ads].keys():
          for toads in alldata[ads]['E_ddag_vs_pot_%s'%facet].keys():
            toads_short = toads.lstrip('md-').lstrip('bdo-')
            for pH in alldata[ads]['E_ddag_vs_pot_%s'%facet][toads].keys():
      #        ads_short = ads.lstrip('md-').lstrip('bdo-')
              ## If no vibrations are found use the IS vibrations
              no_vibs=True
              if 'vibs_ddag_%s'%facet in alldata[ads]:
                  if toads in alldata[ads]['vibs_ddag_%s'%facet]:
                    #if pH == 'chemical':
                    #      print(alldata[ads]['vibs_ddag_%s'%facet][toads])
                    if pH in alldata[ads]['vibs_ddag_%s'%facet][toads].keys():
                      vibens=[vib*invcm for vib in alldata[ads]['vibs_ddag_%s'%facet][toads][pH][1:]]
                      no_vibs=False

              #TODO: Add acidic defaulting to FS not IS!
              if no_vibs:
                  vibens=[vib*invcm for vib in adsvibs]
                  if f'vibs_from_IS_{facet}' not in alldata[ads]:
                      alldata[ads][f'vibs_from_IS_{facet}']={}
                  alldata[ads][f'vibs_from_IS_{facet}'][toads]=True
                  no_freen.append(f'{ads}-{toads}{facet}{pH}')

              if 'free_en_corr_ddag_%s'%facet not in alldata[ads]:
                  alldata[ads]['free_en_corr_ddag_%s'%facet]={}
              if toads not in alldata[ads]['free_en_corr_ddag_%s'%facet]:
                  alldata[ads]['free_en_corr_ddag_%s'%facet][toads]={}
              eng_corr=HarmonicThermo(vibens).get_helmholtz_energy(298,verbose=False)

              if pH == 'base':
                  name_for_reference=ads_short.replace('2','H').replace('3','HH')
              elif pH in ['acid','chemical']:
                  name_for_reference=toads_short.replace('2','H').replace('3','HH')

              outname,dummy,TSnamediff=get_unique_TS_name(ads,toads,pH,[])
              #The free energy contribution of one water molecule is subtracted, because endstates do not contain
              #the water, while barriers need to in order to get the free energy
              if (any(i in outname for i in ['-H2O-','-H-']) or outname == 'H2O-ele') and pH != 'chemical':
                  if not no_vibs: #and (toads not in ['H2CCO'] and not facet == '100'):
                   eng_corr-=clean_G
                  #HCCO to H2CCO is excluded because the desorption not protonation are the barrier
                  if toads in ['H2CCO'] and facet == '100' and not no_water_layer: #not no_vibs:
                      eng_corr+=clean_G

              print(name_for_reference)
              eng_corr-=get_reference_vibrational_contribution(name_for_reference,references=references)
              alldata[ads]['free_en_corr_ddag_%s'%facet][toads][pH]=eng_corr

              if 'G_ddag_vs_pot_%s'%facet not in alldata[ads].keys():
                  alldata[ads]['G_ddag_vs_pot_%s'%facet]={}
              if toads not in alldata[ads]['G_ddag_vs_pot_%s'%facet].keys():
                  alldata[ads]['G_ddag_vs_pot_%s'%facet][toads]={}

              alldata[ads]['G_ddag_vs_pot_%s'%facet][toads][pH] = alldata[ads]['E_ddag_vs_pot_%s'%facet][toads][pH].copy()

              alldata[ads]['G_ddag_vs_pot_%s'%facet][toads][pH][1] += alldata[ads]['free_en_corr_ddag_%s'%facet][toads][pH]
    print('Vibrations could not be found for: ', no_freen)

def get_data(potentials,facets,include_barriers=True,basepath=os.getcwd(),cavity_file='cavity_pot_2.out',
             barrier_path='/Users/geokast/SelectCO2/barriers', outfile = 'results/parsed_data.pckl',plot=False,alldata={}):

    #If a single facet is given as a string

    if isinstance(facets,str):
        facets = [facets]
    for facet in facets:
        print('-'*15)
        print(f'Parsing {facet}')
        facpath=os.path.join(basepath,facet)
        not_enough_potentials=[]
        for ads in os.listdir(facpath):
            final_paths=[]
            adspath=os.path.join(facpath,ads)
            if os.path.isdir(adspath):
                for inads in os.listdir(adspath):
                    if inads.split('_')[0] == ads:
                        final_paths.append(os.path.join(adspath,inads))

            if not len(final_paths):  continue

            if ads not in alldata.keys():
                alldata[ads]={}
            adsdata={}
            for final_path in  final_paths:
                sitename=final_path.split('/')[-1]
                adsdata[sitename] = {}
                for adsfile in os.listdir(final_path):
                    for potential in np.around(potentials,2):
                        if len(adsfile.split('_')) > 1:
                            if '_'.join(adsfile.split('_')[:2]) == 'Pot_%1.2f'%potential:
                                if 'E_%s'%facet not in adsdata[sitename].keys():
                                    adsdata[sitename]['E_%s'%facet]={}
                                    adsdata[sitename]['E_C_%s'%facet]={}
                                if 'ne_%s'%facet not in adsdata[sitename].keys():
                                    adsdata[sitename]['ne_%s'%facet]={}
                                #print(sitename,facet,potential)
                                adsdata[sitename]['E_%s'%facet][potential] = read(final_path+'/'+adsfile).get_potential_energy()
                                adsdata[sitename]['ne_%s'%facet][potential] = read(final_path+'/'+adsfile).calc.results['ne']
                                elpot=read(final_path+'/'+adsfile).calc.results['electrode_potential']
                                adsdata[sitename]['E_C_%s'%facet][adsdata[sitename]['ne_%s'%facet][potential]] =\
                                        read(final_path+'/'+adsfile).get_potential_energy() -\
                                        adsdata[sitename]['ne_%s'%facet][potential]*elpot

                                if ads == 'clean':
                                    adsdata[sitename]['cell'] = read(final_path+'/'+adsfile).cell
                    if adsfile == cavity_file:
                        adsdata[sitename]['cavity_%s'%facet] = np.loadtxt(final_path+'/'+adsfile)

            # Find the most stable site from selection
            #if facet == '211' and ads == 'CO':
                #testprint(alldata[ads])
            #    alldata[ads].update(find_most_stable_site(adsdata,facet,potentials).copy())
            #    testprint(ads,alldata[ads].keys())
            #    das
            #else:
            alldata[ads].update(find_most_stable_site(adsdata,facet,
                potentials,plot=plot,ads=ads,not_enough_potentials=not_enough_potentials).copy())

        E_C_clean_fit_data=[]
        for ne in alldata['clean']['E_C_%s'%facet].keys():
            E_C_clean_fit_data.append([ne,alldata['clean']['E_C_%s'%facet][ne]])
            #print(ne,alldata['clean']['E_C_%s'%facet][ne])
        E_C_clean_fit_data=np.array(E_C_clean_fit_data)
        coeff,d=curve_fit(quad_fun,E_C_clean_fit_data[:,0],E_C_clean_fit_data[:,1])
        alldata['clean']['E_C_vs_ne_%s'%facet]=coeff
        alldata,not_enough_potentials=fit_potential_response(alldata,facet,potentials,include_barriers=False,
                plot=False,specific_ads='clean',quad_fit=True,not_enough_potentials=not_enough_potentials)

        for ads in alldata.keys():
            #Loop over lots of potentials for checking inhomogeneous potential grids
            for potential in np.around(potentials,2):
                if 'E_%s'%facet in alldata[ads].keys() and ads != 'clean':
                    #Reference energies to the clean slab
                    try:
                        #alldata[ads]['E_%s'%facet][potential] -= \#alldata['clean']['E_%s'%facet][potential]
                        alldata[ads]['E_%s'%facet][potential] -= \
                                alldata['clean']['E_vs_pot_%s'%facet][0]*potential**2+\
                                alldata['clean']['E_vs_pot_%s'%facet][1]*potential+\
                                alldata['clean']['E_vs_pot_%s'%facet][2]
                    except KeyError:
                        #print('Potential %1.2f seems  to be missing for adsorbate %s'%(potential,ads))
                        pass
                    else:
                        ads_short=ads.lstrip('md-').lstrip('bdo-').replace('2','H')
                        ads_short=ads_short.replace('3','HH').replace('-uw','')
                        ads_short=ads_short.replace('-hollow','').rstrip('old').lstrip('x-')
                        alldata[ads]['E_%s'%facet][potential] -= get_reference_energies(ads_short,code='GPAW')
                    #Add canonical energies
                    try:
                        alldata[ads]['E_C_%s'%facet][alldata[ads]['ne_%s'%facet][potential]] -=\
                                alldata['clean']['E_C_vs_ne_%s'%facet][0]*alldata[ads]['ne_%s'%facet][potential]**2+\
                                alldata['clean']['E_C_vs_ne_%s'%facet][1]*alldata[ads]['ne_%s'%facet][potential]+\
                                alldata['clean']['E_C_vs_ne_%s'%facet][2]
                    except KeyError:
                        pass
                    else:
                        alldata[ads]['E_C_%s'%facet][alldata[ads]['ne_%s'%facet][potential]] -=\
                                get_reference_energies(ads_short,code='GPAW')

        if include_barriers:
            alldata = get_barriers(alldata,barrier_path,facet,plot_charge_transfer=plot)
        add_gas_phase_data_to_dict(alldata,facet)

        alldata,not_enough_potentials=fit_potential_response(alldata,facet,potentials,include_barriers=True,plot=plot,not_enough_potentials=not_enough_potentials)
        if len(not_enough_potentials):
            print('More potential should be calculated for: ',not_enough_potentials)
        read_vibrational_frequencies(alldata,None,'Cu_surface_sampling_constq.txt',facet)
        add_vibrational_free_energy_corrections(alldata,facet)
        if plot:
            plot_Eddag_from_specific_IS(alldata,'G_ddag',facet,potentials,True,False)

        #testcoeff=alldata['CO']['E_ddag_vs_pot_%s'%facet]['COH']['acid']
        #print('getdata after addvibs CO COH',testcoeff,testcoeff[0]*4.4+testcoeff[1])

        if 0:
            for ads in alldata.keys():
                if 'E_vs_pot_%s'%facet in alldata[ads]:
                    testprint(ads,facet,alldata[ads]['E_vs_pot_%s'%facet])
    import pickle
    out=open(outfile,'wb')
    pickle.dump(alldata,out)

    return alldata

def add_gas_phase_data_to_dict(alldata,facet=None,temperature=298.15,references={'C':'CO'}):
    #Reference gas phase molecules
    from ase.thermochemistry import  IdealGasThermo
    from ase.units import invcm
    for gas in ['H2_g','H2O_g','CO_g','CH4_g','C2H4_g','CH3CH2OH_g']:
        alldata[gas]={}

    alldata['H2_g']['E']=0
    alldata['H2O_g']['E']=0
    alldata['CO_g']['E']=0

    alldata['H2_g']['vibs']=[0, 123.3, 182.2, 304.6, 427.9, 4470.7]
    alldata['CO_g']['vibs']=[0, 0, 78.6, 263.9, 284.3, 2115.0]
    alldata['H2O_g']['vibs']=[0, 0, 121.9, 259.2, 290.5, 358.9, 1623.4, 3753.9,  3874.9]

    alldata['CH4_g']['E']=-2.486809
    alldata['CH4_g']['vibs']=[0, 0, 0,  0,  83.0, 143.7, 1288.3, 1302.9, 1304.3, 1504.7, 1515.7, 3022.0, 3108.2, 3121.0, 3126.1]
    alldata['C2H4_g']['E']=-2.76975805
    alldata['C2H4_g']['vibs']=[0, 0, 0, 0, 0, 0, 826.36096307,  959.45228605,  960.34934469, 1053.97242644, 1229.69837936, 1348.85228724, 1472.24547713, 1668.96347041, 3105.59055667, 3125.65012816, 3173.82607546, 3196.77938303]

    alldata['CH3CH2OH_g']['E']=-3.255991999999992
    alldata['CH3CH2OH_g']['vibs'] = [0,  0,   0,  134.7,  149.6,  181.7,  334.3,  384.9,  495.9,  832.0, 907.2,  992.8, 1081.6, 1157.1, 1271.7, 1295.2, 1403.1, 1426.8, 1465.6, 1502.8, 1519.1, 2959.8, 2971.0, 3032.1, 3073.8, 3090.5, 3799.5]

    alldata['CO_g'].update({'pressure':101325,'geometry':'linear','symmetry': 1})
    alldata['H2_g'].update({'pressure':101325,'geometry':'linear','symmetry': 1})
    alldata['H2O_g'].update({'pressure':0.035*101325,'geometry':'nonlinear','symmetry':2})

    alldata['C2H4_g'].update({'pressure': 1,'geometry':'nonlinear','symmetry': 4})
    alldata['CH3CH2OH_g'].update({'pressure': 1, 'geometry':'nonlinear','symmetry':1})
    alldata['CH4_g'].update({'pressure': 1,'geometry':'nonlinear','symmetry':12})

    for ads in ['H2_g','H2O_g','CO_g','CH4_g','C2H4_g','CH3CH2OH_g']:
        gibbs = IdealGasThermo(vib_energies = np.array(alldata[ads]['vibs'])*invcm,
                                    geometry=alldata[ads]['geometry'],
                                    spin=0,
                                    symmetrynumber=alldata[ads]['symmetry'],
                                    atoms=read('/Users/geokast/SelectCO2/endstates/gas_geometries/'+ads.rstrip('_g')+'.traj'))

        alldata[ads]['free_en_corr'] = gibbs.get_gibbs_energy(
                                                pressure=alldata[ads]['pressure'],
                                                temperature=temperature,
                                                verbose=False)

        alldata[ads]['G'] = alldata[ads]['E']+alldata[ads]['free_en_corr']


    vibnames={'CH4_g':'CHHH','C2H4_g':'CCHHHH','CH3CH2OH_g':'CHHHCHHOH'}
    for ads in ['CH4_g','C2H4_g','CH3CH2OH_g']:
        alldata[ads]['G'] -= get_reference_vibrational_contribution(vibnames[ads],references=references)

    #print('Equilibrium potential of CH4:',-alldata['CH4_g']['G']/6)#*0.255)
    #print('Equilibrium potential of C2H4_g:',-alldata['C2H4_g']['G']/8)#+8*0.14)
    #print('Equilibrium potential of CH3CH2OH_g:',-alldata['CH3CH2OH_g']['G']/8)#+8*0.14)
#    das



def find_most_stable_site(data,facet,potentials,potential_to_check=2.9,plot=False,
        plotdir='results/E_on_varying_sites/',ads=None,quad_fit=True,not_enough_potentials=None):
    E_at_pots=[]
    if not_enough_potentials is None:
        not_enough_potentials=[]
    data,not_enough_potentials=fit_potential_response(data,facet,potentials,include_barriers=False,plot=False,quad_fit=quad_fit,not_enough_potentials=not_enough_potentials)
    coeffs=[]
    for site in data:
        if 'E_vs_pot_%s'%facet not in data[site].keys():continue
        coeff=data[site]['E_vs_pot_%s'%facet].copy()
        #the "fit_potential_response" function defaults to "E_vs_pot" for the fit
        #It's manually renamed here
        data[site]['E_abs_vs_pot_%s'%facet]=data[site]['E_vs_pot_%s'%facet].copy()
        del data[site]['E_vs_pot_%s'%facet]

        if quad_fit:
            E_at_pots.append([site,potential_to_check**2*coeff[0]+potential_to_check*coeff[1]+coeff[2]])
        else:
            E_at_pots.append([site,potential_to_check*coeff[0]+coeff[1]])
        coeffs.append([site,coeff])

    if not len(E_at_pots): return data[site]
    i_most_stable_site=np.argsort(np.array(E_at_pots)[:,1])[-1]

    if plot:
        colors=cm.gist_ncar(np.linspace(0,0.5,len(coeffs)))
        for ic,coeff in enumerate(coeffs):
            #print(coeff,data[coeff[0]]['E_%s'%facet])
            E_at_site=[]
            for pot in np.linspace(2.4,3.4,5):
                E_at_site.append([pot,coeff[1][0]*pot**2+coeff[1][1]*pot+coeff[1][2]])
            E_at_site=np.array(E_at_site)
            if len(coeff[0].split('_')) > 3:
                label = ' '.join(coeff[0].split('_')[3:])
            else:
                label = ' '.join(coeff[0].split('_')[1:])
            plt.plot(E_at_site[:,0],E_at_site[:,1],'-',
                         label=label,
                         color=colors[ic])
            points=[]
            for pot in data[coeff[0]]['E_%s'%facet]:
                points.append([float(pot),data[coeff[0]]['E_%s'%facet][pot]])
            points=np.array(points)
            plt.plot(points[:,0],points[:,1],'o',color=colors[ic],markeredgecolor='k')
        plt.ylabel('$\Omega$ [eV]')
        plt.xlabel('Work function [eV]')
        plt.legend()
        plt.tight_layout()
        plt.savefig(plotdir+'E_vs_site_%s_%s.pdf'%(ads,facet))
        plt.close()

    return data[E_at_pots[i_most_stable_site][0]]


def get_barriers(alldata,barrier_path,facet,plot_charge_transfer=False):
    barrier_dirs=[]
    #Check which barriers might be there
    for bardir in os.listdir(barrier_path):
        if (len(bardir.split('_')) > 2 and bardir != 'CO_to_OCCO'):
            if bardir.split('_')[1] == 'to':
                if facet in os.listdir(barrier_path+'/'+bardir):
                    barrier_dirs.append(barrier_path+'/'+bardir+'/'+facet)

    #Check at which potential the barriers have been calculated
    for bardir_iter in barrier_dirs:
      #Identify IS and FS
      ads=bardir_iter.split('/')[-2].split('_')[0]
      toads=bardir_iter.split('/')[-2].split('_')[2]

      if ads not in alldata.keys():
            if 'bdo-'+ads in alldata.keys():
                ads='bdo-'+ads
            elif 'md-'+ads in alldata.keys():
                ads='md-'+ads
            elif ads == 'COCOH':
                ads='COH'
            elif ads == 'COCHO':
                ads='CHO'
            else:
                print('Could not find the initial state for barrier %s'%bardir_iter.split('/')[-2])
                continue


      if 'E_ddag_%s'%facet not in alldata[ads].keys():
              alldata[ads]['E_ddag_%s'%facet]={}

      toads_short=toads.lstrip('md-').lstrip('bdo-')
      alldata[ads]['E_ddag_%s'%(facet)][toads_short]={}
      ads_short=ads.lstrip('md-').lstrip('bdo-').replace('2','H').replace('3','HH')
      toads_short2=toads.lstrip('md-').lstrip('bdo-').replace('2','H').replace('3','HH')
      letter_diff=toads_short2
      for letter in ads_short:
          letter_diff=letter_diff.replace(letter,'',1)

      for pH in ['base','acid','chemical']:
       charges,energies={},{}
       potentials=[]
       if pH == 'acid':
           if 'acidic' not in os.listdir(bardir_iter):  continue
           bardir=bardir_iter+'/acidic'
       elif pH == 'chemical':
           if 'chemical' not in os.listdir(bardir_iter):  continue
           bardir=bardir_iter+'/chemical'
       else:
           bardir=bardir_iter

       # Check which potentials have been calculated
       for potdir in os.listdir(bardir):
            if len(potdir.split('_')) > 1:
                if potdir.split('_')[0] == 'pot':
                    if 'neb_GC_final_climbed.traj' in os.listdir(bardir+'/'+potdir):
                        potentials.append(potdir.split('_')[1])

       #Determine the difference of IS and FS for formation energies
       for potential in sorted(potentials):
          #print(bardir)
          atoms_list=read(bardir+'/pot_%s/neb_GC_final_climbed.traj@:'%potential)

          Eddag=max([atoms.get_potential_energy() for atoms in atoms_list])
#          if pH == 'acid':
#              print('1',ads,potential,Eddag,alldata['clean']['E_%s'%facet][float(potential)])
          if float(potential) in alldata['clean']['E_%s'%facet]:
              Eddag-=alldata['clean']['E_%s'%facet][float(potential)]
          elif 'E_vs_pot_%s'%facet in alldata['clean']:
              Eddag-= alldata['clean']['E_vs_pot_%s'%facet][0]*float(potential)**2+\
                      alldata['clean']['E_vs_pot_%s'%facet][1]*float(potential)+\
                      alldata['clean']['E_vs_pot_%s'%facet][2]
          else:
              print('Potential %s of the clean slab could not be found for barrier from %s to %s.'%(potential,ads,toads))
              continue

          #It's always the IS we reference to in alkaline
          #TODO: Check whether this works

          #Protonation
          if letter_diff == 'H':
              if pH == 'base':
               reference_string=ads_short#+'HHO'
              elif pH in ['acid','chemical']:
               reference_string=toads_short#+'HHO'
               #reference_string=ads_short#+'HHO'

          else:
              if pH == 'base':
                reference_string=ads_short
              elif pH in ['acid']:
                reference_string=toads_short+'HHO'
              elif pH in ['chemical']:
                  raise NotImplementedError('Chemical OH desorption is not '
                          'implemented')

          if toads_short in ['OCCO']:
              Eddag -= get_reference_energies(toads_short,code='GPAW')
          elif (toads_short in ['OCCOH'] and ads == 'COH') or\
                  (toads_short in ['OCCHO'] and ads == 'CHO'):
              Eddag -= get_reference_energies(toads_short,code='GPAW')
          elif reference_string == 'clean':
              pass
          else: # len(ads_short) > len(toads_short):
              Eddag -= get_reference_energies(reference_string,code='GPAW')

          #if pH == 'acid':
              #print('2',ads,potential,Eddag,alldata[toads]['E_vs_pot_%s'%facet][0]*float(potential)+alldata[toads]['E_vs_pot_%s'%facet][1]-alldata['clean']['E_%s'%facet][float(potential)]-get_reference_energies(toads_short,code='GPAW'))
              #print('2',ads,potential,Eddag,Eddag+(float(potential)-4.4),
              #        alldata[toads]['E_vs_pot_%s'%facet][0]*float(potential)
              #        +alldata[toads]['E_vs_pot_%s'%facet][1]+((float(potential)-4.4)))

          try:
              charges[potential]= [i.calc.results['ne'] for i in atoms_list]
              energies[potential] = [i.calc.results['energy'] for i in atoms_list]
          except:
              pass

          if pH not in alldata[ads]['E_ddag_%s'%facet][toads_short].keys():
              alldata[ads]['E_ddag_%s'%facet][toads_short][pH] = {}
          alldata[ads]['E_ddag_%s'%facet][toads_short][pH][potential]=Eddag
          #if ads=='COCO':
          #  if 'Eddag' not in alldata['CO']:
          #      alldata['CO']['E_ddag_%s'%facet]={toads: {pH: {potential: Eddag}}}
          #  alldata['CO']['E_ddag_%s'%facet][toads][pH][potential]=Eddag
#       if facet == '211':
#           print(facet,ads,toads,potential,alldata[ads]['E_ddag_%s'%facet][toads])
       if len(charges.keys()) and plot_charge_transfer:
          plot_ne_over_band(charges,energies,ads,toads_short,bardir_iter,pH)
      #print(bardir,potentials)
    #das
    #print(alldata['clean'])
    return alldata

def plot_ne_over_band(charges,energies,ads,toads,bardir,pH):
    fig, ax1 = plt.subplots()
    ax1.set_xlabel('#image')
    ax1.set_ylabel(r'$\Delta$E [eV]')
    ax2=ax1.twinx()
    ax2.set_ylabel('Charge transfer [e]')
    for pot in charges.keys():
        ax2.plot(np.array(charges[pot])-charges[pot][0],'--')
        ax1.plot(np.array(energies[pot])-energies[pot][0],label='WF=%1.2feV'%float(pot))
    ax1.set_title(bardir.split('/')[-1])
    fig.tight_layout()
    fig.legend()
    fig.savefig('results/charge_transfer_along_band/charge_transfer_'+bardir.split('/')[-1]+'_'+pH+'.pdf')
    plt.close()

def read_vibrational_frequencies(alldata,line,backup_vibfile=None,facet='100',vibfile='vibrations.pckl',use_HHH_vibs=False,substrates=['Cu'],no_water_layer=False):
    # "no_water_layer" has been added for the calculation without water layer in the barrier calculations

    allvibs=pickle.load(open(vibfile,'rb'))
    #print(allvibs.keys())
    #print(alldata.keys())

    missing_vibs=[]
    used_vibs=[]
    for ads_long in alldata.keys():
     if no_water_layer and ads_long == 'clean_slab': continue
     for pH in ['base']:
        no_vibs=True
        if ads_long[-2:] == '_g': continue
        ads=ads_long.lstrip('md-').lstrip('bdo-')
        # Read vibrational frequencies of thermodynamics
        if ads_long in allvibs.keys() and not use_HHH_vibs:
         if facet in allvibs[ads_long].keys():
          if 'vibs' in allvibs[ads_long][facet].keys():
           if pH in allvibs[ads_long][facet]['vibs'].keys():
            if 'vibs_%s'%facet not in alldata[ads_long].keys():
                alldata[ads_long]['vibs_%s'%facet] = {}
            alldata[ads_long]['vibs_%s'%facet] = [i/units.invcm for i in allvibs[ads_long][facet]['vibs'][pH]]
            no_vibs=False
        #    continue

        # If the vibs are found everything is good and we go to the next
#        print(ads_long)
        if not no_vibs:
            used_vibs.append(ads_long)
            continue

        if backup_vibfile:
            missing_vibs.append(ads_long)
            #print(ads_long+' vibs have not been found in my vibs')
            viblines=open(backup_vibfile,'r').readlines()[1:]
            for ivibline in viblines:
                if ivibline.split()[0] in substrates+['None']:
                    #if ivibline.split()[1] in [facet,'gas']:
                    if ivibline.split()[1] in [facet]:
                        if ivibline.split()[2] == ads:
                           vibline=ivibline
                           break
            else:
      #          print(ads + ' is not in any vibfile')
                continue
                #ads_and_electron[ads]['vibs_%s'%facet] = []
#                return []
        else:
            vibline=line

        if vibline is None:
            continue
        freq_inline=[None,None]
     #   print(ads_long)
        for isplit,splitline in  enumerate(vibline.split()):
            if splitline[0] == '[':
                freq_inline[0]=isplit
            elif splitline[-1] == ']':
                freq_inline[1]=isplit+1
                break

        if None not in freq_inline:
            frequencies = [float(vib.replace(',','').replace('[','').replace(']',''))
                    for vib in vibline.split()[freq_inline[0]:freq_inline[1]]]

        else:
            print('No frequencies given for '+ads)
            frequencies=[]

        alldata[ads_long]['vibs_%s'%facet] = frequencies

    print(f'Missing vibrations for {facet}:', missing_vibs)

    if no_water_layer:
        bars=[barads
              for barads in allvibs.keys()
              if '-' in barads]
    else:
        bars=[barads
              for barads in allvibs.keys()
              if '_to_'  in barads]

    for bar in bars:
        for pH in ['base','acid','chemical']:
            if facet not in allvibs[bar].keys(): continue
            if 'vibs' not in allvibs[bar][facet].keys(): continue
            if pH not in allvibs[bar][facet]['vibs'].keys(): continue
            if not len(allvibs[bar][facet]['vibs'][pH]): continue

            # The following "if no_water_layer" has been added for the calculation without water layer
            if no_water_layer:
              IS_FS={'OC-CO': ['COCO','OCCO'],
                    'OCCO-H': ['OCCO','OCCOH'],
                    'H-HCCO': ['HCCO','H2CCO'],
                    'HCCO-H': ['HCCO','HCCOH'],
                    'HOCCO-H':['OCCOH','HOCCOH'],
                    'OCC-OH':['OCCOH','CCO'],
                    'CO$_{2(g)}$_to_HCOO$^-_{(aq)}':['CO$_{2(g)}$','HCOO$^-_{(aq)}$']}
              barads=IS_FS[bar][0]
              bartoads=IS_FS[bar][1]
            else:
                barads=bar.split('_')[0]
                bartoads=bar.split('_')[-1]
                bartoads=bartoads.lstrip('md-').lstrip('bdo-')
                if barads == 'COCOH':
                    barads='COH'
                elif barads == 'COCHO':
                    barads='CHO'
            #print(alldata.keys())
            if 'vibs_ddag_%s'%facet not in alldata[barads].keys():
                alldata[barads]['vibs_ddag_%s'%facet]={}

            if bartoads not in alldata[barads]['vibs_ddag_%s'%facet].keys():
                alldata[barads]['vibs_ddag_%s'%facet][bartoads]={}
            alldata[barads]['vibs_ddag_%s'%facet][bartoads][pH]=\
                    [i/units.invcm for i in allvibs[bar][facet]['vibs'][pH]]
            used_vibs.append(bar)

    unused_vibs=[]
    for vib in allvibs:
        if vib not in used_vibs:
            unused_vibs.append(vib)
    print(f'Unused vibrations for {facet}:', unused_vibs)

def plot_imaginary_frequencies(alldata, facet,outfilebase='results/img_frequencies.pdf'):
    for pH in ['base','acid','chemical']:
      outfile=outfilebase.rstrip('.pdf')+'_%s.pdf'%pH
      xlabels=[]
      for iads,ads in enumerate(alldata.keys()):
        if 'vibs_ddag_%s'%facet not in alldata[ads].keys(): continue
        for itoads,toads in enumerate(alldata[ads]['vibs_ddag_%s'%facet].keys()):
          #for pH in alldata[ads]['vibs_ddag_%s'%facet][toads].keys():
            if pH not in alldata[ads]['vibs_ddag_%s'%facet][toads]: continue
            if not np.imag(alldata[ads]['vibs_ddag_%s'%facet][toads][pH][0]): continue

            plt.plot(len(xlabels),np.imag(alldata[ads]['vibs_ddag_%s'%facet][toads][pH][0])
                    ,'o',color=get_intcolors(ads,toads),markeredgecolor='k',markersize=8)

            ads_short = ads.lstrip('md-').lstrip('bdo-')
            toads_short = toads.lstrip('md-').lstrip('bdo-')
            xlabels.append(get_intcolors(ads,toads,return_name=True))

      plt.plot(np.nan,np.nan,'o',color=get_intcolors('CO','CHO'),label='C protonation')
      plt.plot(np.nan,np.nan,'o'+get_intcolors('HOCCOH','CCOH'),label='OH desorption')
      plt.plot(np.nan,np.nan,'o'+get_intcolors('OCCO','OCCOH'),label='O protonation')
      plt.plot(np.nan,np.nan,'o'+get_intcolors('clean','H'),label='Other')
      #plt.legend(bbox_to_anchor=(1,1))
      plt.legend()
      plt.xticks(np.arange(len(xlabels)),xlabels,rotation='vertical')
      plt.ylim=[0,1800]
      plt.ylabel('Imaginary frequency [cm$^{-1}$]',fontsize=14)
      plt.tight_layout()
      #plt.show()
      plt.savefig(outfile)
      plt.close()

def apply_Wigner_correction(alldata,facets):
    from ase.units import invcm
    kb=8.617333262145e-5
    T=300.
    hbar=6.582119569e-16
    fig,ax=plt.subplots(1,2,sharex=True,figsize=(12,6))
    fig2,ax2=plt.subplots(1,2,sharex=True,figsize=(12,6))
    second_order_wigfacs,first_order_wigfacs,imgvibs=[],[],[]
    for facet in facets:
     for iads,ads in enumerate(alldata.keys()):
        if 'vibs_ddag_%s'%facet not in alldata[ads].keys(): continue
        for itoads,toads in enumerate(alldata[ads]['vibs_ddag_%s'%facet].keys()):
            if not np.imag(alldata[ads]['vibs_ddag_%s'%facet][toads]['base'][0]): continue
            imgvib=np.imag(alldata[ads]['vibs_ddag_%s'%facet][toads]['base'][0])
            print(ads,toads,imgvib)
            beta_vib=imgvib*invcm/(kb*T)
            print(imgvib,imgvib*invcm,beta_vib,imgvib*invcm/(2*np.pi*kb))

            #Tunneling factor from
            #https://pubs.rsc.org/en/content/articlelanding/2014/CP/C4CP03235G#!divAbstract
            #qi0=1/
            #vibs=np.array(alldata[ads]['vibs_ddag_%s'%facet][toads])


            if 0:
               for N in range(10,100):
                    A_N=1
                    for j in range(3,N-1):

                        eta_0_j = np.sqrt((2*j*np.pi*kb*T/hbar)**2+(imgvib*invcm/(hbar*2*np.pi))**2)
                        print(eta_0_j,kb*T/(hbar*2*np.pi*eta_0_j),hbar*2*np.pi*eta_0_j)
                        A_N*=kb*T/(hbar*2*np.pi*eta_0_j)
               print(A_N)
            #Rate from New Journal ofPhysics 12 (2010) 055002
            om_0=3000
            om_b=1000
            E_b=.7
            T1=300

            for om_b in np.linspace(1,2000,30):
            #for T1 in [100,300,1000]:
                k=om_b*invcm/(4*np.pi*hbar)*\
                  np.sinh(om_0*invcm/(2*kb*T1))/np.sin(om_b*invcm/(2*kb*T1))*np.exp(-E_b/(kb*T1))
                kTST=kb*T1/(2*np.pi*hbar)*np.exp(-E_b/(kb*T1))
                print(T1,k,kTST,k/kTST)
                plt.plot(om_b,k/kTST,'o')
            #plt.ylim([-10,10])
            plt.show()
            asd

            #for k in range(10):
            #    eta_0_k = np.sqrt((2*k*np.pi*kb*T/hbar)**2+(imgvib*invcm/(hbar*2*np.pi))**2)

            #wigfac = beta_vib/2/(np.sin(beta_vib/2))
            crossover_T=imgvib*invcm/(2*np.pi*kb)
            print('Crossover temperature: %s-%s%s, '%(ads,toads,facet), crossover_T)
            second_order_wigfac=1+1/24*(beta_vib)**2
            second_order_wigfacs.append(second_order_wigfac)
            first_order_wigfac=(beta_vib)/2* 1/(np.sin(beta_vib/2.))
            first_order_wigfacs.append(first_order_wigfac)
            imgvibs.append(imgvib)
            ax[1].plot(imgvib,second_order_wigfac,'o',color=get_intcolors(ads,toads))#,label='Second order Wigner')
            ax2[1].plot(imgvib,second_order_wigfac,'o',color=get_intcolors(ads,toads))#,label='Second order Wigner')
            ax[1].plot(imgvib,first_order_wigfac,'d',color=get_intcolors(ads,toads))#,label='Standard Wigner')
            ax2[0].plot(imgvib,first_order_wigfac,'d',color=get_intcolors(ads,toads))#,label='Standard Wigner')
            ax[0].plot(imgvib,crossover_T,'o',color=get_intcolors(ads,toads))#,label='Second order Wigner')
    imgvibs=np.array(imgvibs)
    sort_img=np.argsort(imgvibs)
    imgvibs=imgvibs[sort_img]
    first_order_wigfacs=np.array(first_order_wigfacs)[sort_img]
    second_order_wigfacs=np.array(second_order_wigfacs)[sort_img]

    #Imaginary vibration of crossover at given T
    crossover_vib=(2*np.pi*kb*T)/invcm
    ax[1].axvline(x=crossover_vib,linestyle='--',color='k')
    ax2[0].axvline(x=crossover_vib,linestyle='--',color='k')
    ax2[1].axvline(x=crossover_vib,linestyle='--',color='k')
    ax2[0].annotate(r'$\omega_c$',(1310,1),fontsize=20)
    ax2[1].annotate(r'$\omega_c$',(1310,1),fontsize=20)
#    ax[1].plot(imgvibs,first_order_wigfacs,'--k')
#    ax[1].plot(imgvibs,second_order_wigfacs,'-k')
    ax[0].axhline(y=T,linestyle='--',color='k')
    ax[1].annotate('T=%iK'%int(T),(1300,0.1),fontsize=20)
    ax[0].plot(np.nan,np.nan,'s',color=get_intcolors('CO','CHO'),label='C protonation')
    ax[0].plot(np.nan,np.nan,'s',color=get_intcolors('CO','COH'),label='O protonation')
    ax[0].plot(np.nan,np.nan,'s',color=get_intcolors('COH','C'),label='OH desorption')
    ax[0].plot(np.nan,np.nan,'s',color=get_intcolors('clean','H'),label='Other')
    ax2[0].plot(np.nan,np.nan,'s',color=get_intcolors('CO','CHO'),label='C protonation')
    ax2[0].plot(np.nan,np.nan,'s',color=get_intcolors('CO','COH'),label='O protonation')
    ax2[0].plot(np.nan,np.nan,'s',color=get_intcolors('COH','C'),label='OH desorption')
    ax2[0].plot(np.nan,np.nan,'s',color=get_intcolors('clean','H'),label='Other')
    ax[1].plot(np.nan,np.nan,'ok',label='Second order Wigner')
    ax[1].plot(np.nan,np.nan,'dk',label='Standard Wigner')
    ax2[1].set_title('Second order Wigner')
    ax2[0].set_title('Standard Wigner')
    ax[0].set_ylabel('T$_{cross} [K]$')
    ax[1].set_xlabel('Imaginary frequency [cm$^{-1}$]')
    ax[0].set_xlabel('Imaginary frequency [cm$^{-1}$]')
    ax2[1].set_xlabel('Imaginary frequency [cm$^{-1}$]')
    ax2[0].set_xlabel('Imaginary frequency [cm$^{-1}$]')
    ax[1].set_ylabel('$\kappa_t$')
    ax2[0].set_ylabel('$\kappa_t$')
    ax[1].set_ylim([0,5])

    ax[0].legend()
    ax2[0].legend()
    ax[1].legend()
    fig.tight_layout()
    fig2.tight_layout()
    fig.savefig('results/Wigner_CrossoverT_and_tunneling.pdf')
    fig2.savefig('results/Wigner_tunneling_first_and_second_order.pdf')
    #plt.show()

def plot_Eddag_from_specific_IS(alldata,enkey,facet,potentials,deeperkey=None,quad_fit=False,plot_separate=True,ylabel='Formation energy [eV]'):

    symbols=['+','s','d']
    enkey2=enkey+'_vs_pot'
    for iads,ads in enumerate(alldata.keys()):
     plotted=False
     for quad_fit in [False,True]:
         enkey2=enkey+'_vs_pot'
         if quad_fit: enkey2+='_quad'
         if ads[-2:] == '_g': continue
         if (enkey2+'_%s'%facet not in alldata[ads].keys() or
             ads in ['clean']): continue
         #if enkey+'_vs_pot_%s'%facet in alldata[ads].keys() and ads not in ['clean']:
         if deeperkey:
         #  print(ads,enkey)
           counter=0
           for itoads,toads in enumerate(alldata[ads][enkey2+'_%s'%facet]):
              if enkey2+'_%s'%facet in alldata[ads].keys():
                 counter+=1
           colors=cm.jet(np.linspace(0,1,counter+1))

           counter2=0
           for toads in alldata[ads][enkey2+'_%s'%facet].keys():
            for pH in alldata[ads][enkey2+'_%s'%facet][toads].keys():

             coeff=alldata[ads][enkey2+'_%s'%facet][toads][pH]

             fit = []
             #print(ads,toads,enkey2,coeff,quad_fit)
             if quad_fit:
                 for pot in np.linspace(potentials[0],potentials[-1],10):
                     fit.append([pot,np.array([pot**2,pot,1])@coeff])
             else:
                 for pot in np.linspace(potentials[0],potentials[-1],2):
                     fit.append([pot,np.array([pot,1])@coeff])
             fit=np.array(fit)

             if quad_fit:
                label=ads+'-'+toads+', '+pH+r',d%s/d$\Phi$=(%1.2f,%1.2f)'%(enkey[0],2*coeff[0],coeff[1])
                plt.plot(fit[:,0],fit[:,1],'--',color=colors[counter2%len(colors)],
                     label=label)
             else:
                label=ads+'-'+toads+', '+pH+r',d%s/d$\Phi$=%1.2f'%(enkey[0],coeff[0])
                plt.plot(fit[:,0],fit[:,1],'-',color=colors[counter2%len(colors)],
                     label=label)
             try:
                E_v_pot=get_E_vs_pot(alldata,ads,potentials,enkey+'_%s'%facet,[toads,pH])
             except:
                 pass
             else:
                if len(E_v_pot):
                 plt.plot(E_v_pot[:,0],E_v_pot[:,1],
                         symbols[counter2%len(symbols)],
                         #label=ads+'-'+toads+', '+pH+r',dE/d$\Phi$=%1.2f'%coeff[0],
                         color=colors[counter2%len(colors)])
             plotted=True
             counter2+=1
     if not plotted: continue
     plt.xlabel('Work function [eV]')
     plt.ylabel(ylabel)
     plt.title('Barriers from %s on facet %s'%(ads,facet))
     #plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
     #           mode="expand", borderaxespad=0, ncol=4,fontsize=3)
     plt.legend()
     plt.tight_layout()
     plt.savefig('results/E_ddag_vs_pot_from_specific_ads/%s_vs_pot_%s_%s.pdf'%(enkey,ads,facet))
     plt.close()

def testprint(*string):
    print(*string)
