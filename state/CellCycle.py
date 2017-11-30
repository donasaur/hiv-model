# -*- coding: utf-8 -*-
import numpy as np

from process.RevBinding import RevBinding
from mainaux.InitParamValues import *

#This is a type of State Class
class CellCycle(object):
    def __init__(self, param_dict=None):
        if param_dict==None:
            param_dict = generate_param_dict();
        
        #Initialize necessary parameters
        
        #Initialize Variables
        self.cell_cycle_arrest = 0      
        
    def record_state(self, record, timestep, max_timesteps):
        record.add_tracking(timestep, max_timesteps, 'cell_cycle_arrest', self.cell_cycle_arrest)

    def record_at_end(self, record):
        pass  
