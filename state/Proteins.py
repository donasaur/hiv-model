# -*- coding: utf-8 -*-
import numpy as np
from mainaux.InitClassVars import *

#This is a type of State Class
class Proteins(object):
    index = initialize_protein_dict()

    def __init__(self):
        #Initialize necessary parameters
        self.proteins_cyt = np.zeros((9), int) 
        #[0]: Vif
        #[1]: Vpr
        #[2]: Tat
        #[3]: Env
        #[4]: Rev
        #[5]: Nef
        #[6]: Gag
        #[7]: Gag/Pro/Pol
        #[8]: Gag_Dimers
        self.proteins_nuc = np.zeros((9), int)
        self.proteins_mem = np.zeros((9), int)
        self.proteins_virion = np.zeros((9), int)
        self.env_misc = dict()
        # self.env_misc['gp160'] = 0
        # self.env_misc['gp160 : ER'] = 0
        # self.env_misc['gp160 : ER : with DOLPP sugar'] = 0 # G1
        # self.env_misc['gp160 : ER : 2 sugar truncate'] = 0 # G2 and G3
        # self.env_misc['gp160 : ER : 2 sugar truncate-folded'] = 0
        # self.env_misc['gp160 : ER : 3 sugar truncate'] = 0 # G$
        # self.env_misc['gp160 : ER : overtruncate'] = 0
        # self.env_misc['gp160 : Golgi'] = 0
        # self.env_misc['gp160-trimer : Golgi'] = 0
        # self.env_misc['gp41'] = 0
        # self.env_misc['gp120'] = 0
        # self.env_misc['gp41 : Mem'] = 0
        # self.env_misc['gp41: Virion'] = 0
        # self.env_misc['gp120 : Mem'] = 0
        # self.env_misc['gp120 : Virion'] = 0
        self.env_misc['Env : cytoplasm'] = 0
        self.env_misc['Env : ER'] = 0
        self.env_misc['Env : ER : G1'] = 0
        self.env_misc['Env : ER : G2'] = 0
        self.env_misc['Env : ER : G3'] = 0
        self.env_misc['Env : ER : G3 : folded'] = 0
        self.env_misc['Env : ER : G4 : folded'] = 0
        self.env_misc['Env : Golgi'] = 0
        self.env_misc['Env : Golgi : G5'] = 0
        self.env_misc['Env : Golgi : G5 : error'] = 0
        self.env_misc['Env : trimers'] = np.zeros((4), int)
        self.env_misc['Env : trimers : cleaved'] = np.zeros((4), int)
        self.env_misc['Env : trimers : membrane'] = np.zeros((4), int)

    @staticmethod
    def initialize_class_vars():
        rtn_dict = dict()
        rtn_dict['Vif'] = 0
        rtn_dict['Vpr'] = 1
        rtn_dict['Tat'] = 2
        rtn_dict['Env'] = 3
        rtn_dict['Rev'] = 4
        rtn_dict['Nef'] = 5
        rtn_dict['Gag'] = 6
        rtn_dict['GagProPol'] = 7
        rtn_dict['Gag_dimers'] = 8
        return rtn_dict       
        
    def record_state(self, record, timestep, max_timesteps):
        record.add_tracking(timestep, max_timesteps, 'proteins_cyt', self.proteins_cyt)
        record.add_tracking(timestep, max_timesteps, 'proteins_nuc', self.proteins_nuc)
        record.add_tracking(timestep, max_timesteps, 'proteins_mem', self.proteins_mem)
        record.add_tracking(timestep, max_timesteps, 'proteins_virion', self.proteins_virion)
        record.add_tracking(timestep, max_timesteps, 'env_cyt', self.proteins_cyt[Proteins.index['Env']])
        record.add_tracking(timestep, max_timesteps, 'env_ER', self.env_misc['Env : ER'])
        record.add_tracking(timestep, max_timesteps, 'env_ER_G1', self.env_misc['Env : ER : G1'])
        record.add_tracking(timestep, max_timesteps, 'env_ER_G2', self.env_misc['Env : ER : G2'])
        record.add_tracking(timestep, max_timesteps, 'env_ER_G3', self.env_misc['Env : ER : G3'])
        record.add_tracking(timestep, max_timesteps, 'env_ER_G3_folded', self.env_misc['Env : ER : G3 : folded'])
        record.add_tracking(timestep, max_timesteps, 'env_ER_G4_folded', self.env_misc['Env : ER : G4 : folded'])
        record.add_tracking(timestep, max_timesteps, 'env_Golgi', self.env_misc['Env : Golgi'])
        record.add_tracking(timestep, max_timesteps, 'env_Golgi_G5', self.env_misc['Env : Golgi : G5'])
        record.add_tracking(timestep, max_timesteps, 'env_Golgi_G5_error', self.env_misc['Env : Golgi : G5 : error'])
        record.add_tracking(timestep, max_timesteps, 'env_trimers', self.env_misc['Env : trimers'])
        record.add_tracking(timestep, max_timesteps, 'env_trimers_cleaved', self.env_misc['Env : trimers : cleaved'])
        record.add_tracking(timestep, max_timesteps, 'env_trimers_membrane', self.env_misc['Env : trimers : membrane'])

    def record_at_end(self, record):
        pass
        
