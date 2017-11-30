# -*- coding: utf-8 -*-
"""
This is a process class whose evolve_state method is called
at each timestep.

An instance of ProteinLocalization is initialized once per Simulation
with the State as input. Only the state 'protein' is modified in this
process.

At each timestep, for a protein x, for all x, evolve_state reads the amount
of x in the nucleus and cytoplasm and writes back new values of x in the
nucleus and cytoplasm depending on how much of x is set to be transferred
from one location to the other.

Translocation abundances are randomly determined by probabilities based on rates 
of translocation found in the literature. 

Note: only Rev and Tat shuffling are currently supported.

Summary of the biology:
This process takes into account the shuffling of proteins between the
nucleus and the cytoplasm. 
The proteins Rev and Tat are synthesized in the cytoplasm, but function in the
nucleus--Rev to aid in mRNA export and Tat to increase the viral transcription
frequency. 
At certain rates, they may transition in and out of the cytoplasm. Note that
in the current framework, these are the only proteins that than translocate to
the nucleus. 

#Todo: Check reference 1 page number
#References:
#1. Pond, S. J. K., Ridgeway, W. K., Robertson, R., Wang, J., & Millar, D. P. (2009). HIV-1 Rev protein assembles on viral RNA one molecule at a time. Proceedings of the National Academy of Sciences of the United States of America, 106(5), 14048. doi:10.1073/pnas.0807388106
#2. Kim, H., Yin, J. (2005) In silico mutagenesis of RNA Splicing in HIV-1. Biotechnology and bioengineering 91: 877-893.
#3. Lamond AI, Sleeman JE. Nuclear substructure and dynamics. Curr Biol. 2003 Oct 28 13(21):R825-8.


"""
import numpy as np
from mainaux.State import State
from mainaux.Process import Process
from state.Proteins import Proteins
from mainaux.InitParamValues import *

#This is a Process Class
class ProteinLocalization(Process):
    def __init__(self, state, param_dict=None):
        self.state = state

        if param_dict==None:
            param_dict = generate_param_dict();
                
        #Constant parameters
        self.PROB_REV_SHUTTLING_IN = param_dict['PROB_REV_SHUTTLING_IN'] #1/min C-->N #Kim/Yin 2005
        self.PROB_REV_SHUTTLING_OUT = param_dict['PROB_REV_SHUTTLING_OUT'] #1/min N-->C #Kim/Yin 2005
        self.PROB_TAT_SHUTTLING_OUT = param_dict['PROB_TAT_SHUTTLING_OUT'] #Kim/Yin 2005
        #self.PROB_TAT_SHUTTLING_OUT = 0.1 
        self.PROB_TAT_SHUTTLING_IN = param_dict['PROB_TAT_SHUTTLING_IN'] #1/min #Kim/Yin 2005
                

    def evolve_state(self, timestep):
        #Rev Binding Timescale:
        #The rev binding and dissociating rates are very fast, and I am unable to see reasonable dynamics 
        #at a 1 min timestep becuase of the multiple events that should be occuring withing the timestep
        #Therefore "zooming in time" to a 1 second timestep for Rev binding--making DETERMINISTIC
        #is there a way to add stochasticity?
        #How to discretize??
    
        #get variables
        protein_state = self.state.get_state('proteins')
        proteins_nuc = protein_state.proteins_nuc
        proteins_cyt = protein_state.proteins_cyt
             
        #evolve state
        
        #Rev Shuttling
        #Arbirary assumption of a poisson distribution
        Rev_shuttling_in = np.min([np.random.poisson(proteins_cyt[Proteins.index['Rev']]*self.PROB_REV_SHUTTLING_IN), proteins_cyt[Proteins.index['Rev']]])
        proteins_cyt[Proteins.index['Rev']]=proteins_cyt[Proteins.index['Rev']]-Rev_shuttling_in
        proteins_nuc[Proteins.index['Rev']]=proteins_nuc[Proteins.index['Rev']]+Rev_shuttling_in
        Rev_shuttling_out = np.min([np.random.poisson(proteins_nuc[Proteins.index['Rev']]*self.PROB_REV_SHUTTLING_OUT), proteins_nuc[Proteins.index['Rev']]])
        proteins_cyt[Proteins.index['Rev']]=proteins_cyt[Proteins.index['Rev']]+Rev_shuttling_out
        proteins_nuc[Proteins.index['Rev']]=proteins_nuc[Proteins.index['Rev']]-Rev_shuttling_out
        ####This may be super inefficient
    #    tempRand = np.random.rand(proteins_cyt[Proteins.index['Rev']])
    #    proteins_cyt[Proteins.index['Rev']]=proteins_cyt[Proteins.index['Rev']]-((tempRand<probRev_shuttling_in).sum())
    #    proteins_nuc[Proteins.index['Rev']]=proteins_nuc[Proteins.index['Rev']]+((tempRand<probRev_shuttling_in).sum())
    #    tempRand = np.random.rand(proteins_nuc[Proteins.index['Rev']])
    #    proteins_cyt[Proteins.index['Rev']]=proteins_cyt[Proteins.index['Rev']]+((tempRand<self.PROB_REV_SHUTTLING_OUT).sum())
    #    proteins_nuc[Proteins.index['Rev']]=proteins_nuc[Proteins.index['Rev']]-((tempRand<self.PROB_REV_SHUTTLING_OUT).sum())
        
        #Tat Shuttling
        #Arbirary assumption of a poisson distribution
        Tat_shuttling_in = np.min([np.random.poisson(proteins_cyt[Proteins.index['Tat']]*self.PROB_TAT_SHUTTLING_IN), proteins_cyt[Proteins.index['Tat']]])
        proteins_cyt[Proteins.index['Tat']]=proteins_cyt[Proteins.index['Tat']]-Tat_shuttling_in
        proteins_nuc[Proteins.index['Tat']]=proteins_nuc[Proteins.index['Tat']]+Tat_shuttling_in
        #Currently the shutling out probability is set to 0, so this code is commented out
        Tat_shuttling_out = np.min([np.random.poisson(proteins_nuc[Proteins.index['Tat']]*self.PROB_TAT_SHUTTLING_OUT), proteins_nuc[Proteins.index['Tat']]])
        proteins_cyt[Proteins.index['Tat']]=proteins_cyt[Proteins.index['Tat']]+Tat_shuttling_out
        proteins_nuc[Proteins.index['Tat']]=proteins_nuc[Proteins.index['Tat']]-Tat_shuttling_out               
                    
        #write back parameters to state object
#        protein_state.protein_nuc = proteins_nuc
#        protein_state.protein_cyt = proteins_cyt
        
        #update state to new values
        # self.state.set_state('proteins', protein_state)
