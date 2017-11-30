# -*- coding: utf-8 -*-
import numpy as np

from state.Proteins import Proteins
from state.MRNAs import MRNAs
from state.HostFactors import HostFactors
from state.ReactionRates import ReactionRates
from state.DNAs import DNAs
from state.ViralParticles import ViralParticles
from state.CellCycle import CellCycle
from state.ViralProgeny import ViralProgenyContainer

class State(object,):
    def __init__(self, param_dict=None):
        #initialize all the states that you need
        self.states_dict = {}
        protein_state = Proteins()
        self.states_dict['proteins'] = protein_state
        mRNA_state = MRNAs(param_dict)
        self.states_dict['mRNAs'] =  mRNA_state
        host_factor_state = HostFactors()
        self.states_dict['host_factors'] = host_factor_state
        reaction_rate_state = ReactionRates()
        self.states_dict['reaction_rates'] = reaction_rate_state     
        DNA_state = DNAs()
        self.states_dict['DNAs'] = DNA_state  
        Viral_particles_state = ViralParticles()
        self.states_dict['viral_progeny'] = Viral_particles_state
        Cell_cycle_state = CellCycle()
        self.states_dict['cell_cycle'] = Cell_cycle_state
        Viral_progeny_container_state = ViralProgenyContainer(self)
        self.states_dict['viral_progeny_container'] = Viral_progeny_container_state

        
    def get_state(self, state_name):
        #Get the current value of a given state
        return self.states_dict[state_name]
        
    def set_state(self, state_name, state):
        self.states_dict[state_name] = state
