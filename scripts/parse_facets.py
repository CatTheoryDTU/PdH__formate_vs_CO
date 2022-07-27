import sys

from general_imports import *
from sjm_analyse_tools import read_vibrational_frequencies,add_vibrational_free_energy_corrections,testprint,add_gas_phase_data_to_dict


# alldata tree: adsorbates, facet, sites ('most_stable' is most stable site), scanning potential, quantities
# ads_and_electrons tree: adsorbates, quantities_with_facetname_in the end (only most stable site)

def main():
    alldata={}
    parse_data(alldata,ads_and_electron=ads_and_electron)

def parse_data(alldata,specific_ads=None,ads_and_electron=None,scanning_potential=3.4):

        if specific_ads is not None:
            if isinstance(specific_ads,str): adsorbates=[specific_ads]
            elif isinstance(specific_ads,list): adsorbates=specific_ads
            else:
                print('Do not understand specific ads given to parse_data in scripts/parse_facets.py')
                exit()
        else:
            adsorbates=[]
            for fac in os.listdir(rundir):
                facpath=rundir+'/'+fac
                if not os.path.isdir(facpath): continue
                if len(fac) != 3: continue
                if not fac.isnumeric(): continue

                for ads in os.listdir(facpath):
                    if ads[0] in ['m','b']: continue

                    if ads not in adsorbates and ads not in ['clean_slab']:
                        adsorbates.append(ads)
        parse_dirs(home,alldata=alldata)
        for ads in adsorbates:
            if '-' in ads:
                parse_dirs(home,ads,alldata=alldata,barrier=True)
            else:
#                continue
                parse_dirs(home,ads,alldata=alldata,barrier=False)

        get_most_stable_binding(alldata)

        if ads_and_electron is not None:
            transform_alldata_for_old_routines(alldata,ads_and_electron,scanning_potential)

        # Calculating the free energies
        for ifac,facet in enumerate(alldata['CO'].keys()):
            read_vibrational_frequencies(ads_and_electron,None,None,facet=facet,vibfile='results/vibrations.pckl',no_water_layer=True)
            add_vibrational_free_energy_corrections(ads_and_electron,facet,no_water_layer=True,references={'C':'CO2'})

        # All empirical corrections are in the following
        surfH_free_en_corr=0.14 #Value gotten from H desoprtion study, used for subtracting at the TS, where H is now part of the adsorbate
        # Manually calculating the free energy correction of the TS to formate
        from ase.thermochemistry import HarmonicThermo
        from ase.units import invcm
        from general_tools import get_reference_vibrational_contribution
        ads,toads,pH='CO$_{2(g)}$','HCOO$^-_{(aq)}$','base'
        adsvibs=ads_and_electron[ads]['vibs_ddag_%s'%facet][toads][pH]
        vibens=[vib*invcm for vib in adsvibs]
        ads_and_electron[ads]['free_en_corr_ddag_%s'%facet]={toads:{}}

        eng_corr=HarmonicThermo(vibens).get_helmholtz_energy(298,verbose=False)
        eng_corr-=get_reference_vibrational_contribution('HCOO',references={'C':'CO2'})
        #Subtracting the vibrational contribution of one H on the surface as it is now part of the adsorbate and not the slab anymore
        eng_corr-=surfH_free_en_corr
        ads_and_electron[ads]['free_en_corr_ddag_%s'%facet][toads][pH]=eng_corr

        ads_and_electron[ads]['G_ddag_vs_pot_%s'%facet]={toads:{pH:ads_and_electron[ads]['E_ddag_vs_pot_%s'%facet][toads][pH].copy()}}
        ads_and_electron[ads]['G_ddag_vs_pot_%s'%facet][toads][pH][1] += ads_and_electron[ads]['free_en_corr_ddag_%s'%facet][toads][pH]


        ads_and_electron['CO$_{2(g)}$']['E_111']=[0,0.0]
        ads_and_electron['CO$_{2(g)}$']['E_vs_pot_111']=[0,0.0]
        ads_and_electron['CO$_{2(g)}$']['G_111']=[0,0]
        ads_and_electron['CO$_{2(g)}$']['G_vs_pot_111']=[0,0]
        ads_and_electron['CO$_{(g)}$']={'nHe':2}

        # Actual calculated Formation free energy of CO leading to a determined correction factor of 0.33
        ads_and_electron['CO$_{(g)}$']['G_vs_pot_%s'%facet]=[0,0.7354919708262795]
        #Applying the CO2 gasphase correction to *CO and CO_g
        ads_and_electron['CO$_{(g)}$']['G_vs_pot_%s'%facet][1]-=0.33
        ads_and_electron['CO']['G_vs_pot_%s'%facet][1]-=0.33

        ads_and_electron['HCOOH$_{(aq)}$']={'nHe':2}
        ads_and_electron['HCOOH$_{(aq)}$']['G_vs_pot_%s'%facet]=[0,0.4440000000000026]

        pkl.dump(alldata,open(pkldir,'wb'))
        if ads_and_electron is not None:
            pkl.dump(ads_and_electron,open(pkldir2,'wb'))
        return

def transform_alldata_for_old_routines(alldata,ads_and_electron,scanning_potential):
    print('Creating oldstyle dict of results for using with scripts made for 100 study')
#ads_and_electron['CO'].keys()
#dict_keys(['E_100', 'E_C_100', 'ne_100', 'E_C_vs_pot_100', 'E_abs_vs_pot_100', 'E_ddag_100', 'E_vs_pot_100', 'E_ddag_vs_pot_100', 'E_ddag_rel_100', 'E_ddag_rel_vs_pot_100', 'E_ddag_vs_pot_quad_100', 'E_ddag_rel_vs_pot_quad_100', 'vibs_100', 'vibs_ddag_100', 'free_en_corr_100', 'G_vs_pot_100', 'free_en_corr_ddag_100', 'G_ddag_vs_pot_100', 'E_111', 'E_C_111', 'ne_111', 'E_C_vs_pot_111', 'E_abs_vs_pot_111', 'E_ddag_111', 'E_vs_pot_111', 'E_ddag_vs_pot_111', 'E_ddag_rel_111', 'E_ddag_rel_vs_pot_111', 'E_ddag_vs_pot_quad_111', 'E_ddag_rel_vs_pot_quad_111', 'vibs_111', 'vibs_ddag_111', 'free_en_corr_111', 'G_vs_pot_111', 'vibs_from_IS_111', 'free_en_corr_ddag_111', 'G_ddag_vs_pot_111', 'E_211', 'E_C_211', 'ne_211', 'E_C_vs_pot_211', 'E_abs_vs_pot_211', 'E_ddag_211', 'E_vs_pot_211', 'E_ddag_vs_pot_211', 'E_ddag_rel_211', 'E_ddag_rel_vs_pot_211', 'E_ddag_vs_pot_quad_211', 'E_ddag_rel_vs_pot_quad_211', 'vibs_211', 'vibs_ddag_211', 'free_en_corr_211', 'G_vs_pot_211', 'free_en_corr_ddag_211', 'G_ddag_vs_pot_211', 'E_110', 'E_C_110', 'ne_110', 'E_C_vs_pot_110', 'E_abs_vs_pot_110', 'E_ddag_110', 'E_vs_pot_110', 'E_ddag_vs_pot_110', 'E_ddag_rel_110', 'E_ddag_rel_vs_pot_110', 'E_ddag_vs_pot_quad_110', 'E_ddag_rel_vs_pot_quad_110', 'vibs_110', 'vibs_ddag_110', 'free_en_corr_110', 'G_vs_pot_110', 'free_en_corr_ddag_110', 'G_ddag_vs_pot_110'])

#alldata['OC-CO']['100']['barrier'].keys()
#dict_keys(['E', 'pot_v_q', 'E_v_pot', 'Pot_v_q', 'Eb', 'Eb_v_pot', 'Eb_v_pot_quad'])
    #print(alldata.keys())
    ##CO2gas energy is added manually as it is 0 throughout

    add_gas_phase_data_to_dict(alldata)
    for ads in alldata:
        #print(ads)
        if ads not in ads_and_electron: # and 'barrier' not in alldata[ads]['100']:
            print(f'WARNING! {ads} is not in ./scripts/intermediate_dict.py')
            ads_and_electron[ads]={}
        if ads[-2:] == '_g':
            for quant in alldata[ads]:
                ads_and_electron[ads][quant]=alldata[ads][quant]
            continue
        for facet in alldata[ads]:
            #print(facet)
            if ads == 'clean_slab':
                for quant in alldata[ads][facet]['slab']:
                    ads_and_electron[ads][quant+'_'+facet]=alldata[ads][facet]['slab'][quant]
                continue

            #Translating barriers from being there own entries to being connected to the IS
            IS_FS={'OC-CO': ['COCO','OCCO'],
                    'OCCO-H': ['OCCO','OCCOH'],
                    'H-HCCO': ['HCCO','H2CCO'],
                    'HCCO-H': ['HCCO','HCCOH'],
                    'HOCCO-H':['OCCOH','HOCCOH'],
                    'OCC-OH':['OCCOH','CCO']}

            if ads in IS_FS.keys():
                if ads == 'CO-CO': dd
                out_to_in_quants={f'E_ddag_vs_pot_{facet}':'Eb_v_pot',
                        f'E_ddag_rel_vs_pot_{facet}': 'E_ddag_rel'
                        }
                for outquant in [f'E_ddag_vs_pot_{facet}',f'E_ddag_rel_vs_pot_{facet}']:
                    if outquant not in ads_and_electron[IS_FS[ads][0]]:
                        ads_and_electron[IS_FS[ads][0]][outquant]={}
                    ads_and_electron[IS_FS[ads][0]][outquant].update({IS_FS[ads][1]:
                        {'base':alldata[ads][facet]['barrier'][out_to_in_quants[outquant]]}})
                    if facet == '310':
                        print(facet,IS_FS[ads])
                        testprint(ads_and_electron[IS_FS[ads][0]][outquant])
    #                dsa

#                if f'E_ddag_rel_vs_pot_{facet}' not in ads_and_electron[IS_FS[ads][0]]:
#                    ads_and_electron[IS_FS[ads][0]][f'E_ddag_vs_pot_{facet}']={}
            if 'barrier' in alldata[ads][facet]:
                if ads=='H-CO2':
                    #print(alldata['H-CO2'][facet].keys())
                    ads_and_electron['CO$_{2(g)}$']['E_ddag_vs_pot_%s'%facet]={'HCOO$^-_{(aq)}$':
                            {'base':alldata['H-CO2'][facet]['barrier']['Eb_v_pot']}}
#                    print(ads_and_electron['CO$_{2(g)}$']['E_ddag_vs_pot_%s'%facet])
 #                   das
 #                  # ads_and_electron['CO2']['G_ddag_vs_pot_%s'%facet]={'HCOO':
                   #         {'base':alldata['H-CO2']['Eb_v_pot']}
                   # print(ads)
                   # sss
                continue
            #print(ads,facet)
            for inquant in alldata[ads][facet]['most_stable'][scanning_potential]:
                outquant=inquant
                if inquant == 'Eb_v_pot': outquant='E_vs_pot'
                if inquant == 'E_v_pot': outquant='E_abs_vs_pot'
                if inquant == 'Pot_v_q': outquant='ne'
                if inquant == 'Eb': outquant='E'

                ads_and_electron[ads][outquant+'_'+facet]=\
                        alldata[ads][facet]['most_stable'][scanning_potential][inquant]
#            if ads=='CO2':
#                print(ads_and_electron['CO2'])
            #    asd
            #print(ads_and_electron['CH4_g'])
            #gasdict={}


def parse_dirs(home,system='clean_slab',oldstyle=False,barrier=False,alldata={},specific_facets=None):
    rundir = pj(home,'data')
    for ifac,facet in enumerate(os.listdir(rundir)):
        if facet in ['gas_phase']: continue
        if specific_facets is not None:
            if isinstance(specific_facets,str):
                specific_facets=[specific_facets]
            if facet not in specific_facets: continue

        facdir = pj(rundir,facet)
        if not os.path.isdir(facdir): continue
        sysdir = pj(facdir,system)
        if not os.path.isdir(sysdir): continue
        if barrier:
            finaldir=sysdir
            site='barrier'
            if read_trajectories(finaldir,alldata,facet,system,site,barrier):
                calculate_binding_energy(alldata,facet,system,site,barrier)
        else:
            conformers=[]
            for sitedir0 in os.listdir(sysdir):
                if sitedir0 in ['vacuum','vibs','gas_phase']:
                    continue
                #print(sysdir,sitedir0)
                conformer=sitedir0.split('_')[1][:-7]
                if conformer not in conformers:
                    conformers.append(conformer)
            if not ifac:
                print(f'{system}s conformers are: ',conformers)

            for sitedir0 in os.listdir(sysdir):
                if sitedir0 in ['vacuum','vibs']:
                    continue

                conformer=sitedir0.split('_')[1][:-7]
                finaldir=pj(sysdir,sitedir0)
                #print(sitedir0)
                if not os.path.isdir(finaldir): continue

                if system == 'clean_slab':
                    if sitedir0.split('_')[1] == 'with-water':
                        finaldir1=finaldir
                        for site_ang in os.listdir(finaldir1):
                            finaldir=pj(finaldir1,site_ang)
                            site=site_ang
                            read_trajectories(finaldir,alldata,facet,system,site,barrier)
                            calculate_binding_energy(alldata,facet,system,site,barrier)
                            continue
                    else:
                        site=sitedir0.split('_')[1]
                else:
                    if len(conformers) == 1:
                        site=sitedir0.split('_')[-1]
                    else:
                        site = conformer+'_'+sitedir0.split('_')[-1]

                if barrier:
                    site='barrier'
                #print('a',system)
                if read_trajectories(finaldir,alldata,facet,system,site,barrier):
                    #dd
                    #print(system)
                    calculate_binding_energy(alldata,facet,system,site,barrier)
#                print(sitedir0,finaldir)


def read_trajectories(path,alldata,facet,system,site=None,barrier=False):
   E_v_pot,Erel_v_pot,q_v_pot=[],[],[]

   for files in os.listdir(path):
       if barrier:
           if not (files[:4] == 'pot_'): continue
           infile=files+'/neb_GC_final_climbed.traj'
           if 'neb_GC_final_climbed.traj' not in os.listdir(path+'/'+files):
               print(path+'/'+infile+ 'has a directory, but the converged traj is missing')
               continue
           if system in ['HCCO-H'] and facet == '100' and files in ['pot_2.40','pot_2.65']:
                baratoms=read(path+'/'+infile,index='6:')
           elif system in ['H-HCCO'] and facet == '100':
                baratoms=read(path+'/'+infile,index=':6')
           else:
               baratoms = read(path+'/'+infile,index=':')

           # Only keep the TS as the atoms object for getting the absolute energy
           atoms = baratoms[np.argsort([i.get_potential_energy() for i in baratoms])[-1]]
       else:
            if not (files.split('_')[0] == 'CO2' and files.split('.')[-1] == 'traj'): continue
            #aa
            infile=files
            try:
                atoms = read(path+'/'+infile)
            except:
               print(path+'/'+files, 'is broken, because of disk quota mess')
               continue

       #Check if OCCO split
       Cpos=[]
       if system == 'OCCO':
           for iatm,atom in enumerate(atoms.get_chemical_symbols()):
               if atom == 'C': Cpos.append(atoms.positions[iatm])
           if  np.linalg.norm(Cpos[1]-Cpos[0]) > 2:
                #print(system, facet, site, ' split!')
                continue
       try:
           E_v_pot.append([atoms.calc.results['electrode_potential'],atoms.get_potential_energy()])
       except KeyError:
           raise KeyError(f'Electrode potential not found in {path},\
                   most likely the wrong ase version is used. Did you load_catmap?')

       #For barriers also get the energies relative to the IS of the band
       if barrier:
               Erel_v_pot.append([atoms.calc.results['electrode_potential'],
                   atoms.get_potential_energy()-baratoms[0].get_potential_energy()])

       Acell=np.product(np.diag(atoms.cell[:2,:2]))
       q_v_pot.append([atoms.calc.results['electrode_potential'],-atoms.calc.results['ne']*1.6022*1e3/Acell])


   if not len(E_v_pot): return None
   E_v_pot=np.array(E_v_pot)
   if len(E_v_pot) < 3:
      return None

   coeff,dummy = curve_fit(quad_fun,E_v_pot[:,0],E_v_pot[:,1])
       #coeff,dummy = curve_fit(lin_fun,E_v_pot[:,0],E_v_pot[:,1])
   q_v_pot=np.array(q_v_pot)
   #Cap_coeff,dummy = curve_fit(quad_fun,q_v_pot[:,0],q_v_pot[:,1])
   Cap_coeff,dummy = curve_fit(lin_fun,q_v_pot[:,0],q_v_pot[:,1])

   if system  not in alldata:
      alldata[system] = {}
   if facet not in alldata[system]:
      alldata[system][facet]={}
   if site not in alldata[system][facet]:
      alldata[system][facet][site]={}
   aldat=alldata[system][facet][site]

   aldat['E'] = E_v_pot
   aldat['pot_v_q'] = q_v_pot
   aldat['E_v_pot']=coeff
   aldat['Pot_v_q']=Cap_coeff#-aldat['Pot_v_q'][facet]['slab']['Pot_v_q']

   # Ading the barriers relative to the IS to the dict
   if len(Erel_v_pot):
        Erel_v_pot=np.array(Erel_v_pot)
        aldat['E_ddag_rel'] = Erel_v_pot
        if len(E_v_pot) < 3:
            return True

        aldat['E_ddag_rel_vs_pot'],dummy = curve_fit(lin_fun,Erel_v_pot[:,0],Erel_v_pot[:,1])

   return True

def calculate_binding_energy(alldata,facet,system,site=None,barrier=False):
   if system == 'clean_slab' and site == 'slab': return None
   Eb=[]

   if site not in alldata[system][facet]: return None
   #print(facet,site)
   aldat=alldata[system][facet][site]

   refname=system
   if barrier:
       for rep in (('-H',''),('H-','')):
        refname=refname.replace(*rep)

   for rep in (('H2','HH'),('H3','HHH'),('md-',''),('bd-',''),('bdo-',''),('-',''),('O2','OO')):
        refname=refname.replace(*rep)

   # This was implemented for the study of the best water cluster,
   # For referencing to the slab with water use E not Eb of it
   if refname == 'clean_slab':
        refname = 'HHHHHHHHOOOO'
   #print(system,refname,site)

   if barrier and system not in ['OC-CO','H-CO2']:
       if facet == '100':
        slabsite='0.167-0.125_60'
       elif facet == '310':
        slabsite='0-0.2_90'
       elif facet == '320':
        slabsite='0.083-0.25_30'

   else:
       slabsite='slab'
#   slabsite='slab'
   slab_coeff=alldata['clean_slab'][facet][slabsite]['E_v_pot']
   #print(slab_coeff)
   for pot_E in aldat['E']:#E_v_pot:
        slabE=slab_coeff[0]*pot_E[0]**2+slab_coeff[1]*pot_E[0]+slab_coeff[2]
        #print(system)
        if system == 'clean_minus_H':
            eb = pot_E[1] - slabE +get_reference_energies('H',code='GPAW',references={'C':'CO2'})
        elif system == 'H-CO2':
            eb = pot_E[1] - slabE -get_reference_energies('COO',code='GPAW',references={'C':'CO2'})
        else:
            eb = pot_E[1] - slabE - get_reference_energies(refname,code='GPAW',references={'C':'CO2'})
        Eb.append([pot_E[0],eb])

   aldat['Eb']=np.array(Eb)
   if len(aldat['Eb']) >=3:
       lincoeff,dummy = curve_fit(lin_fun,aldat['Eb'][:,0],aldat['Eb'][:,1])
   else:
       print(f'Could not fit the binding energy of {system} on {facet}, check {path}')
       return

   aldat['Eb_v_pot']=lincoeff
   aldat['Eb_v_pot_quad']=aldat['E_v_pot']-alldata['clean_slab'][facet]['slab']['E_v_pot']
#   print(system,site,get_reference_energies(system,code='GPAW'))
   if system == 'clean_minus_H':
       aldat['Eb_v_pot_quad'][-1]+=get_reference_energies('H',code='GPAW',references={'C':'CO2'})
   elif system == 'H-CO2':
       aldat['Eb_v_pot_quad'][-1]-=get_reference_energies('COO',code='GPAW',references={'C':'CO2'})
   else:
       aldat['Eb_v_pot_quad'][-1]-=get_reference_energies(refname,code='GPAW',references={'C':'CO2'})
   #if system == 'H-CO2':
   #    print(system,aldat)
   #    dsa
   return



def get_most_stable_binding(alldata,outfile='results/most_stable_sites.list'):
    outstring=''
    out=open(outfile,'w')
    for system in alldata:
       if system == 'clean_slab': continue
       if '-' in system: continue
       for ifac,facet in enumerate(alldata[system]):
           facdat=alldata[system][facet]
           if 'most_stable' not in facdat:
               facdat['most_stable'] = {}#[minsite,minEb]
           for scanning_potential in [2.90,3.15,3.40,3.65,3.90]:
            minEb=1000
            minsites={}
            for isite,site in enumerate(facdat):
                  sitedat=facdat[site]
                  if 'Eb_v_pot' not in sitedat: continue
                   #print(sitedat['Eb_v_pot'])
                  Eb=sitedat['Eb_v_pot'][0]*scanning_potential+sitedat['Eb_v_pot'][1]
                  if Eb < minEb:
                       minsite=site
                       minEb=Eb
                       #minsite=site
            facdat['most_stable'][scanning_potential] = alldata[system][facet][minsite].copy()
            facdat['most_stable'][scanning_potential]['site'] = minsite
           outstring+=f'{system} most stable site on {facet} @ {scanning_potential:.2f}: {minsite}, {minEb:.4f}\n'
           #print(system, ' most stable site on %s @ %1.2f: %s, %1.4f'%(facet, scanning_potential, minsite,minEb))
    print(f'Most stable sites written into {outfile}')
    out.write(outstring)
    out.close()



if __name__ == '__main__':
    main()






