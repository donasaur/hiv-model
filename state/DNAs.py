# -*- coding: utf-8 -*-
import numpy as np

#This is a type of State Class
class DNAs(object):
    def __init__(self):
        #Initialize necessary parameters
        self.promoter_activity = 0 #0=off, 1=on       
        
    def record_state(self, record, timestep, max_timesteps):
        record.add_tracking(timestep, max_timesteps, 'promoter_activity', self.promoter_activity)

    def record_at_end(self, record):
        pass

