# -*- coding: utf-8 -*-
"""
This is a process class whose evolve_state method is called
at each timestep.

An instance of Degradation is initialized once per Simulation
with the State as input. However, only the states 'proteins' and
'mRNAs' are modified in this process.

At each timestep, evolve_state reads in the abundance arrays of
protein_state (2) and mRNA_state (6), runs the process of Degradation
to calculate the new abundances, and writes back the values to
the respective buckets of each abundance array.

degrade(self, deg_rate, abundances, released_factors=None, num_constructs)
is the function that is used to calculate the new abundance arrays of the
states of interest. It takes in a deg_rate and abundances array
and outputs a new abundances array and released_factors count.
The released_factors count outputted is only significant if there
is a released_factors input, signifying that certain proteins
were released during the process of degradation. A num_constructs can
also be specified if there are multiple types of a mRNA with the same
Rev count.

Note: only Rev release is currently supported.

Summary of the biology:
This process takes into account the degradation step of the Tat/Rev
replication cycle. The setting for this process is in the nucleus and
cytoplasm. Different proteins as well as multi-spliced, single-spliced,
and full length mRNAs get degraded in this process. This step also introduces
certain proteins into the system as some factors get released prior to
mRNA degradation.

References:
#1. Kim, H., Yin, J. (2005) In silico mutagenesis of RNA Splicing in HIV-1. Biotechnology and bioengineering 91: 877-893.    

"""

import numpy as np
from mainaux.State import State
from mainaux.Process import Process
from state.Proteins import Proteins
from mainaux.InitParamValues import *

#This is a Process Class
class Degradation(Process):
    def __init__(self, state, param_dict=None):
        self.state = state

        if param_dict==None:
            param_dict = generate_param_dict();
                
        #Constant parameters
        self.PROB_mRNA_DEG = param_dict['PROB_mRNA_DEG'] #1/min #Kim, Yin 2005
        self.PROB_PROTEIN_DEG_NUC = param_dict['PROB_PROTEIN_DEG_NUC'] #1/min #Kim, Yin 2005
        self.PROB_PROTEIN_DEG_CYT = param_dict['PROB_PROTEIN_DEG_CYT'] #1/min #Kim, Yin 2005
        self.PROB_PROTEIN_DEG_MEM = param_dict['PROB_PROTEIN_DEG_MEM']
        # Currently, assuming no protein degradation occurs for proteins_virion

    def degrade(self, deg_rate, abundances, released_factors=None, num_constructs = None):
        #Inputs:
        #deg_rate = rate of degradation of X
        #abundances = vector holding starting abundances of different forms of X
        #released_factors = if there is some factor Y that is released when Xi is degraded, the number of Y in a cell
        #num_constructs = number of different constructs in the Abundances vector (does NOT include differing numbers of bound factors) E.g. num_constructs = 7 for singleSpliceTranscripts
    
        #Outputs:
        #Returns a tuple in the form of [new array of counts, new released_factors count]
        #Assumptions:
        #This function assumes that your degradation rate is less than 1
        #Assumes no more than 1 type of released factor
        # np.seed(i)
        if np.sum(abundances)>=1: #if at least one X exists
            for i in range(np.size(abundances)): #for each form of X: X[i]

                temp_rand = np.random.rand(abundances[i]) #generate a rand array for each item that can potentially be degraded
                #count how many are less than the degradation probability and decrement
                decrement_amount = (temp_rand < deg_rate).sum()
                abundances[i]=abundances[i]- decrement_amount
                
                # increase released_factors amount accordingly
                if released_factors != None: #if additional factors are bound to X[i], release and account for them
                    released_factors = released_factors+decrement_amount*(np.divide(i,num_constructs))
        return [abundances, released_factors]
        
    def evolve_state(self, timestep):

        #get variables
        protein_state = self.state.get_state('proteins')
        proteins_nuc = protein_state.proteins_nuc # proteins_nuc[Proteins.index['Rev']] is Rev count in nucleus
        proteins_cyt = protein_state.proteins_cyt # proteins_cyt[Proteins.index['Rev']] is Rev count in cytoplasm
        proteins_mem = protein_state.proteins_mem
        
        mRNA_state = self.state.get_state('mRNAs')
        full_len_transcripts_nuc = mRNA_state.full_len_transcripts_nuc
        full_len_transcripts_cyt = mRNA_state.full_len_transcripts_cyt
        single_splice_transcript_nuc = mRNA_state.single_splice_transcript_nuc
        single_splice_transcript_cyt = mRNA_state.single_splice_transcript_cyt
        multi_splice_transcript_nuc = mRNA_state.multi_splice_transcript_nuc
        multi_splice_transcript_cyt = mRNA_state.multi_splice_transcript_cyt
        
        #evolve state
        #Will any existing transcripts/proteins be degraded?    
        #assume that Rev-bound transcripts degrade at the same rate as Rev-free transcripts
        #A. Full mRNA:
        [full_len_transcripts_nuc, proteins_nuc[Proteins.index['Rev']]] = self.degrade(self.PROB_mRNA_DEG, full_len_transcripts_nuc, proteins_nuc[Proteins.index['Rev']], 1)        
        [full_len_transcripts_cyt, proteins_cyt[Proteins.index['Rev']]] = self.degrade(self.PROB_mRNA_DEG, full_len_transcripts_cyt, proteins_cyt[Proteins.index['Rev']], 1)
        #B. Single-spliced mRNA:
        [single_splice_transcript_nuc, proteins_nuc[Proteins.index['Rev']]] = self.degrade(self.PROB_mRNA_DEG, single_splice_transcript_nuc, proteins_nuc[Proteins.index['Rev']], 7)        
        [single_splice_transcript_cyt, proteins_cyt[Proteins.index['Rev']]] = self.degrade(self.PROB_mRNA_DEG, single_splice_transcript_cyt, proteins_cyt[Proteins.index['Rev']], 7)
        #C. Multi-spliced mRNA:
        [multi_splice_transcript_nuc, temp] = self.degrade(self.PROB_mRNA_DEG, multi_splice_transcript_nuc)        
        [multi_splice_transcript_cyt, temp] = self.degrade(self.PROB_mRNA_DEG, multi_splice_transcript_cyt)
        #D. Proteins:
        [proteins_nuc, temp] = self.degrade(self.PROB_PROTEIN_DEG_NUC, proteins_nuc)        
        [proteins_cyt, temp] = self.degrade(self.PROB_PROTEIN_DEG_CYT, proteins_cyt)
        [proteins_mem, temp] = self.degrade(self.PROB_PROTEIN_DEG_MEM, proteins_mem)
       
        #write back parameters to state object
        protein_state.proteins_nuc = proteins_nuc
        protein_state.proteins_cyt = proteins_cyt
        protein_state.proteins_mem = proteins_mem
        mRNA_state.full_len_transcripts_nuc = full_len_transcripts_nuc
        mRNA_state.full_len_transcripts_cyt = full_len_transcripts_cyt
        mRNA_state.single_splice_transcript_nuc = single_splice_transcript_nuc
        mRNA_state.single_splice_transcript_cyt = single_splice_transcript_cyt
        mRNA_state.multi_splice_transcript_nuc = multi_splice_transcript_nuc 
        mRNA_state.multi_splice_transcript_cyt = multi_splice_transcript_cyt
        
        # I think this might be unnecessary. protein_state, mRNA_state are mutable objects.
        #update state to new values
        #self.state.set_state('proteins', protein_state)
        #self.state.set_state('mRNAs', mRNA_state)
