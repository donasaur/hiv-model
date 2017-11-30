# -*- coding: utf-8 -*-
import numpy as np

from process.RevBinding import RevBinding
from mainaux.InitParamValues import *

#This is a type of State Class
class ViralParticles(object):
    def __init__(self, param_dict=None):
        if param_dict==None:
            param_dict = generate_param_dict();
        
        #Initialize Variables
        self.viral_progeny = dict()
        #in each viral progeny:
        #item 0: state
        #   0 = simple nucleate on the membrane
        #   1 = growing puncta adding Gag and other proteins
        #   2 = complete viral particle ready for budding
        #   3 = budded viral particle
        #item 1: # of Gag
        #item 2: growth_constant
        self.viral_progeny_keys = []     
        
    def record_state(self, record, timestep, max_timesteps):
        #record.add_tracking(timestep, max_timesteps, 'viral_progeny', self.viral_progeny)
        record.add_tracking(timestep, max_timesteps, 'viral_progeny_keys', self.viral_progeny_keys)

    def record_at_end(self, record):
        pass
