
ads_and_electron={
        'clean_minus_H':{'nHe':-1,'n_sites':1,'leads_to':['clean_slab']},
        'CO$_{2(g)}$':{'nHe':0,'n_sites':1,'leads_to':['CO2','HCOO$^-_{(aq)}']},
        'CO2':{'nHe':0,'n_sites':1,'leads_to':['COOH','HCOO']},
        'COOH':{'nHe':1,'n_sites':1,'leads_to':['CO']},
        'HCOO$^-_{(aq)}$':{'nHe':2,'n_sites':1,'leads_to':[]},
       'CO':{'nHe':2,'n_sites':1,'leads_to':['COCO','CHO','COH','H']},
        'H':{'nHe':1,'n_sites':1,'leads_to':['H2','CHO','COH']},
        'CO$_{g}$':{'nHe':2,'leads_to':[],'n_sites':0},
        }
