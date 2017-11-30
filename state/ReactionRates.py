# -*- coding: utf-8 -*-
import numpy as np

#This is a type of State Class
class ReactionRates(object):
    def __init__(self):
        #Initialize necessary parameters
        self.Tat_derived_transcription_rate = 0
        self.translation_suppressed = 0    
        
    def record_state(self, record, timestep, max_timesteps):
        record.add_tracking(timestep, max_timesteps, 'Tat_derived_transcription_rate', self.Tat_derived_transcription_rate)    

    def record_at_end(self, record):
        pass
        
