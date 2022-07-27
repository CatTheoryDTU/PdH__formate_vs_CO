import os
from ase.io import read,write
import numpy as np

def get_reference_vibrational_contribution(adsorbate,CO_as_ref=True,code='GPAW',references={'C':'CO'}):
    E={}
    if code == 'GPAW':
        E['O'] = get_gas_reference_vibrational_contributions('H2O')-get_gas_reference_vibrational_contributions('H2')
        E['H'] = 0.5*get_gas_reference_vibrational_contributions('H2')
        E['CO'] = get_gas_reference_vibrational_contributions('CO')
        Cref=''
        for ilet, letter in enumerate(references['C']):
            if letter.isnumeric():
                if ilet == 0:
                    print ('C reference starts with a number, can not be interpreted')
                    return False
                for nat in range(int(letter)-1):
                    Cref+=references['C'][ilet-1]
            else:
                Cref+=letter
        E['C']=get_gas_reference_vibrational_contributions(references['C'],code='GPAW')
        for letter in Cref:
            if letter == 'C': continue
            E['C']-=E[letter]
    else:
        print('Only GPAW is implemented for vibrational contribution of gas species at the moment')
        sys.exit()

    E_ref=0
    for elem in adsorbate:
        E_ref += E[elem]

    return E_ref

def get_gas_reference_vibrational_contributions(gas,code='GPAW'):
    if code == 'GPAW':
        gases={'CH4': 0.6609615722468455, #1atm
                'H2':-0.05276775930248084, #1atm
                'H2O': 0.0071360210575864835, #0.035atm
                'CO': -0.38954084777493236, #1atm
                'CO2':-0.25199999999999534 #1atm
                }
    elif code == 'VASP':
        gases={'CH4': 0.6609615722468455, #1atm
                'H2':-0.035, #1atm
                'H2O': 0.009, #0.035atm
                'CO': -0.371 #1atm
                }
    else:
        raise InputError('Code not understood')
    return gases[gas]

def get_reference_energies(adsorbate,CO_as_ref=True,code='VASP',use_potential_energy=True,
        references={'C':'CO'}):
    E = {}
    #OCCO (g)
    if code == 'QE':
        E['C'] = get_gas_reference('CH4')-2*get_gas_reference('H2')
        E['O'] = get_gas_reference('H2O')-get_gas_reference('H2')
        E['H'] = 0.5*get_gas_reference('H2')
        E['CO'] = get_gas_reference('CO')
    elif code == 'VASP':
        E['H'] = 0.5*get_gas_reference_VASP('H2')
        E['O'] = get_gas_reference_VASP('H2O')-get_gas_reference_VASP('H2')
        Cref=''
        for ilet, letter in enumerate(references['C']):
            if letter.isnumeric():
                if ilet == 0:
                    print ('C reference starts with a number, can not be interpreted')
                    return False
                for nat in range(int(letter)-1):
                    Cref+=references['C'][ilet-1]
            else:
                Cref+=letter
        E['C']=get_gas_reference_VASP(references['C'])
        for letter in Cref:
            if letter == 'C': continue
            E['C']-=E[letter]

    elif code == 'GPAW':
        E['O'] = get_gas_reference_GPAW('H2O',use_potential_energy)-get_gas_reference_GPAW('H2',use_potential_energy)
        E['H'] = 0.5*get_gas_reference_GPAW('H2',use_potential_energy)
        Cref=''
        for ilet, letter in enumerate(references['C']):
            if letter.isnumeric():
                if ilet == 0:
                    print ('C reference starts with a number, can not be interpreted')
                    return False
                for nat in range(int(letter)-1):
                    Cref+=references['C'][ilet-1]
            else:
                Cref+=letter
        E['C']=get_gas_reference_GPAW(references['C'])
        for letter in Cref:
            if letter == 'C': continue
            E['C']-=E[letter]
        E['CO'] = get_gas_reference_GPAW('CO',use_potential_energy)

    E_ref=0
    for elem in adsorbate:
        E_ref += E[elem]

    return E_ref

def get_gas_reference_GPAW(gas,use_potential_energy=True):
    if use_potential_energy:
        gases = { 'CH4': -34.082626231,
                  'H2': -8.009969,
                  'H2O': -27.231564,
                 'CO': -34.797848,
                 'CO2': -54.767 #+ 0.33 Empirical correction only for OCO backbone
                  }
    else:
        gases={'CH4': -33.422,
                'H2':-8.063,
                'H2O':-27.156,
                'CO': -35.187,
                'CO2': -55.019  #+0.33#Empirical correction for OCO backbone
                }
    return gases[gas]

def lin_fun(x,a,b):
    return a*x+b

def quad_fun(x,a,b,c):
    return a*x*x+b*x+c

