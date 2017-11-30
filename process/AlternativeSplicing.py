# -*- coding: utf-8 -*-
"""
This is a process class whose evolve_state method is called
at each timestep.

An instance of AlternativeSplicing is initialized once per Simulation
with the State as input. The states 'mRNAs' and 'proteins' are modified
in this process.

At each timestep, evolve_state reads in the abundances of full-length,
singly-spliced, and multi-spliced mRNAs as well as the Rev count and
changes the values of the bins depending on which mRNA type gets spliced.
Each bin has some of its mRNA transferred to another bin, except for the
bins that are in their terminal form and represent mRNAs that cannot be
spliced any further. The Rev count is increased accordingly when a
SD4-SA7 splice event takes place.

The new mRNA abundances plus the free Rev count is then written back
to the corresponding State classes at the end of the timestep.

This event takes place in the nucleus.

Summary of the biology:
This process takes into account the process of splicing among the
different types of mRNAs that result in decreased abundances of some
types of mRNA and increased abundances of other types.
"""

import numpy as np
from mainaux.State import State
from mainaux.Process import Process
from process.RevBinding import RevBinding
from mainaux.InitParamValues import *

#References:
#1. Kim, H., Yin, J. (2005) In silico mutagenesis of RNA Splicing in HIV-1. Biotechnology and bioengineering 91: 877-893.
#2. Hammond, B. J. (1993). Quantitative Study of the control of Hiv-1 gene expression. J. Theor. Biol., 163, 199–221.

#This is a Process Class
class AlternativeSplicing(Process):
    def __init__(self, state, param_dict=None):
        self.state = state

        if param_dict==None:
            param_dict = generate_param_dict();

        #Constant parameters
        #From Hammond 1993 and Kim/Yin 2005
        self.PROB_SPLICE_FULL_TO_SINGLE = param_dict['PROB_SPLICE_FULL_TO_SINGLE'] #Hammond 1993 #probability of a D1-->A1,A2,A3,4abc,A5 splice
        self.PROB_SPLICE_SINGLE_TO_MULTI = param_dict['PROB_SPLICE_SINGLE_TO_MULTI'] #probability of a D4-->A7 splice
        self.PROB_VPR_THIRD_SPLICE = param_dict['PROB_VPR_THIRD_SPLICE'] #probability of a D3--> A3,4abc,5 splice #2/3 chance of a fully spliced vpr being further spliced to tat/rev/nef over its ~4h lifespan ***NEED to fit
        self.PROB_VIF_THIRD_SPLICE = param_dict['PROB_VIF_THIRD_SPLICE'] #probability of a D3--> A3,4abc,5 splice #100% chance of a fully spliced vif being further spliced to tat/rev/nef over its ~4h lifespan ***NEED to fit
        F1 = 0.01 #probability of selecting A1 in first splice event
        F2 = 0.02 #probability of selecting A2 in first splice event
        F3 = 0.1 #probability of selecting A3 in first splice event
        F4 = 0.13 #probability of selecting A4a,b,c in first splice event
        F5 = 0.74 #probability of selecting A5 in first splice event
#        F1 = 0.011 #probability of selecting A1 in first splice event
#        F2 = 0.022 #probability of selecting A2 in first splice event
#        F3 = 0.114 #probability of selecting A3 in first splice event
#        F4 = 0.01 #probability of selecting A4a,b,c in first splice event
#        F5 = 0.843 #probability of selecting A5 in first splice event
        ##F1 = 0.01 #probability of selecting A1 in first splice event #Kim/Yin 2005
        ##F2 = 0.02 #probability of selecting A2 in first splice event #Kim/Yin 2005
        ##F3 = 0.05 #probability of selecting A3 in first splice event #Kim/Yin 2005
        ##F4 = 0.13 #probability of selecting A4a,b,c in first splice event #Kim/Yin 2005
        ##F5 = 0.79 #probability of selecting A5 in first splice event #Kim/Yin 2005
        self.CUMULATIVE_F1_TO_F5 = np.cumsum([F1, F2, F3, F4, F5])
        self.CUMULATIVE_F3_TO_F5 = np.cumsum([(float(F3)/(F3+F4+F5)), (float(F4)/(F3+F4+F5)), (float(F5)/(F3+F4+F5))])
        self.SPLICE_DELAY_FACTOR = param_dict['SPLICE_DELAY_FACTOR'] #kim et al. 2005 = 0.8  
        #print self.PROB_VIF_THIRD_SPLICE
        
        self.MAX_REV_PER_TRANSCRIPT = param_dict['MAX_REV_PER_TRANSCRIPT']        
        
        #OLD CODE--REMOVE??
        #freqSpliceSingleToMulti = 4.5/60.0 #probability of a D4-->A7 splice
        #self.PROB_VPR_THIRD_SPLICE = (2.0/3.0)/(4.0*60.0) #probability of a D3--> A3,4abc,5 splice #2/3 chance of a fully spliced vpr being further spliced to tat/rev/nef over its ~4h lifespan ***NEED to fit
        #self.PROB_VIF_THIRD_SPLICE = (1.0)/(4.0*60.0) #probability of a D3--> A3,4abc,5 splice #100% chance of a fully spliced vif being further spliced to tat/rev/nef over its ~4h lifespan ***NEED to fit
                
    def splice(self, unspliced_abundances, spliced_abundances, splice_probability, cum_prob_splice_forms, splice_form_indexes, indexing_factor, rev_bound_indexes, splice_delay_factor, released_factor=None):
        
        #unspliced_abundances = vector holding starting abundances of different forms of transcripts TO BE spliced
        #spliced_abundances = vector holding starting abundances of different forms of transcripts that HAVE BEEN spliced
        #splice_probability = probability of a splice event
        #cum_prob_splice_forms = a vector of the cumulative probilities of the given splice forms a trasncript can take, given that a splice event occurs
        #splice_form_indexes = vector (same len as cum_prob_splice_forms) holding indexes into spliced_abundances for each resulting splice form 
        #indexing_factor = constant to properly index into vectors such as singleSpliceTranscripts that is ~ numConstructs x numRev
        #rev_bound_indexes = vector of indexes of the startingAbundences vector representing Rev-bound constructs (form = np.array), set to -1 if none
        #splice_delay_factor = factor by which to reduce the probability of a splice event if Rev protein is bound to the transcript
        #released_factor = abundance of protein that is released from transcript upon splicing (Rev)
        if np.sum(unspliced_abundances)>=1: #if there are any transcripts to be spliced
            for j in np.where(unspliced_abundances>0)[0]: #for each non-zero bin (j is the index of each bin)
                temp_rand = np.random.rand(unspliced_abundances[j]) #generate a rand num for each transcript
                #apply a delay factor to Rev bound species...that is Rev inhibits splicing. 
                if np.size(rev_bound_indexes)>0 and j in rev_bound_indexes: 
                    d=splice_delay_factor
                else:
                    d=0
                # temp_value represents the amt of a certain mRNA to be spliced for sure
                temp_value = ((temp_rand<(splice_probability*(1-d))).sum())
                #count how many are less than the splice probability, and do accounting
                unspliced_abundances[j]=unspliced_abundances[j]-temp_value
                
                #if you need to allocate the splicing amoung different forms
                if np.size(cum_prob_splice_forms)>1:
                    #for each new splice event, determine which acceptor site was selected
                    temp_rand2 = np.random.rand(temp_value) 
                    for i in range(temp_value):
                        for k in range(np.size(cum_prob_splice_forms)):
                            if temp_rand2[i] < cum_prob_splice_forms[k]: # if probability of particular mRNA instance is less than threshold probability to take on another form
                                spliced_abundances[(j*indexing_factor)+splice_form_indexes[k]]+=1 #then do accounting
                                break
                #if you don't need to allocate the splicing amoung different forms
                else:
                    spliced_abundances = spliced_abundances + temp_value
                    if released_factor != None:
                        released_factor = released_factor+(temp_value*j)
                        #print('released rev', j)
        return [unspliced_abundances, spliced_abundances, released_factor]

    def evolve_state(self, timestep):    
        #get variables
        protein_state = self.state.get_state('proteins')
        proteins_nuc = protein_state.proteins_nuc
        
        mRNA_state = self.state.get_state('mRNAs')
        full_len_transcripts_nuc = mRNA_state.full_len_transcripts_nuc
        single_splice_transcript_nuc = mRNA_state.single_splice_transcript_nuc
        multi_splice_transcript_nuc = mRNA_state.multi_splice_transcript_nuc
        
        #evolve state
        
        ##########
        #Splice 1#
        ##########     
        #Will any existing transcripts be spliced (first splice event)?
        #splice_form_indexes = [0, 1, 2, 3, 6] = [A1 vif, A2 vpr, A3 tat, A4abc env, A5 env]
        [full_len_transcripts_nuc, single_splice_transcript_nuc, released_factors] = self.splice(full_len_transcripts_nuc, single_splice_transcript_nuc, self.PROB_SPLICE_FULL_TO_SINGLE, self.CUMULATIVE_F1_TO_F5, np.array([0, 1, 2, 3, 6]), 7, np.arange(1,self.MAX_REV_PER_TRANSCRIPT+1), self.SPLICE_DELAY_FACTOR)

        #################
        #Step4: Splice 2#
        #################
        #Will any existing transcripts be spliced (second splice event)?  
        #(1) vif single spliced --> vif multi splice 
        #vif is indexed at arange(0,57,7)
        #This splice event releases bound Rev molecules 
        #resulting form: #D1-A1, D4-A7 #vif                 
        [unspliced_abundances, multi_splice_transcript_nuc[0], released_factors] = self.splice(single_splice_transcript_nuc[np.arange(0,63,7)], multi_splice_transcript_nuc[0], self.PROB_SPLICE_SINGLE_TO_MULTI, 1, 0, 7, np.arange(1,self.MAX_REV_PER_TRANSCRIPT+1), self.SPLICE_DELAY_FACTOR, proteins_nuc[4])               
        single_splice_transcript_nuc[np.arange(0,63,7)] = unspliced_abundances #assign the new vif values
        proteins_nuc[4] = released_factors
        #(2) vpr single spliced --> vpr multi splice 
        #vpr is indexed at arange(1,58,7)
        #This splice event releases bound Rev molecules 
        #resulting form: #D1-A2, D4-A7 #vpr               
        [unspliced_abundances, multi_splice_transcript_nuc[6], released_factors] = self.splice(single_splice_transcript_nuc[np.arange(1,63,7)], multi_splice_transcript_nuc[6], self.PROB_SPLICE_SINGLE_TO_MULTI, 1, 0, 7, np.arange(1,self.MAX_REV_PER_TRANSCRIPT+1), self.SPLICE_DELAY_FACTOR, proteins_nuc[4])               
        single_splice_transcript_nuc[np.arange(1,63,7)] = unspliced_abundances #assign the new vif values
        proteins_nuc[4] = released_factors
        #(3) tat single spliced --> tat double spliced at the frequecy of the 2nd splice event 
        [unspliced_abundances, multi_splice_transcript_nuc[12], released_factors] = self.splice(single_splice_transcript_nuc[np.arange(2,63,7)], multi_splice_transcript_nuc[12], self.PROB_SPLICE_SINGLE_TO_MULTI, 1, 0, 7, np.arange(1,self.MAX_REV_PER_TRANSCRIPT+1), self.SPLICE_DELAY_FACTOR, proteins_nuc[4])               
        single_splice_transcript_nuc[np.arange(2,63,7)] = unspliced_abundances #assign the new tat values
        proteins_nuc[4] = released_factors
        #(4) env single spliced --> rev,nef double spliced at the frequecy of the 2nd splice event
        [unspliced_abundances, multi_splice_transcript_nuc[13], released_factors] = self.splice(single_splice_transcript_nuc[np.arange(3,63,7)], multi_splice_transcript_nuc[13], self.PROB_SPLICE_SINGLE_TO_MULTI, 1, 0, 7, np.arange(1,self.MAX_REV_PER_TRANSCRIPT+1), self.SPLICE_DELAY_FACTOR, proteins_nuc[4])               
        single_splice_transcript_nuc[np.arange(3,63,7)] = unspliced_abundances #assign the new rev values
        proteins_nuc[4] = released_factors        
        [unspliced_abundances, multi_splice_transcript_nuc[16], released_factors] = self.splice(single_splice_transcript_nuc[np.arange(6,63,7)], multi_splice_transcript_nuc[16], self.PROB_SPLICE_SINGLE_TO_MULTI, 1, 0, 7, np.arange(1,self.MAX_REV_PER_TRANSCRIPT+1), self.SPLICE_DELAY_FACTOR, proteins_nuc[4])               
        single_splice_transcript_nuc[np.arange(6,63,7)] = unspliced_abundances #assign the new nef values
        proteins_nuc[4] = released_factors        
        
        #################
        #Step5: Splice 3#
        #################        
        #(1b) vif multi splice  --> tat, rev, nef double spliced 
        #splice_form_indexes = [1, 2, 5] = [A3 tat, A4abc rev, A5 nef]
        #resulting forms: 
        #   #D1-A1, D2-A3, D4-A7 tat
        #   #D1-A1, D2-A4abc, D4-A7 rev (lumping abc for now)
        #   #D1-A1, D2-A5, D4-A7 nef
        [unspliced_abundances, multi_splice_transcript_nuc, released_factors] = self.splice(np.array([multi_splice_transcript_nuc[0]]), multi_splice_transcript_nuc, self.PROB_VIF_THIRD_SPLICE, self.CUMULATIVE_F3_TO_F5, np.array([1,2,5]), 0, [], self.SPLICE_DELAY_FACTOR)
        multi_splice_transcript_nuc[0] = unspliced_abundances[0]

        #(2b) vpr multi splice  --> tat, rev, nef double spliced 
        #splice_form_indexes = [7, 8, 11] = [A3 tat, A4abc rev, A5 nef]
        #resulting forms: 
        #   #D1-A2, D3-A3, D4-A7 tat
        #   #D1-A2, D3-A4abc, D4-A7 rev (lumping abc for now)
        #   #D1-A2, D3-A5, D4-A7 nef
        [unspliced_abundances, multi_splice_transcript_nuc, released_factors] = self.splice(np.array([multi_splice_transcript_nuc[6]]), multi_splice_transcript_nuc, self.PROB_VIF_THIRD_SPLICE, self.CUMULATIVE_F3_TO_F5, np.array([7,8,11]), 0, [], self.SPLICE_DELAY_FACTOR)
        multi_splice_transcript_nuc[6] = unspliced_abundances[0]
         
       
                    
        #write back parameters to state object
        protein_state.protein_nuc = proteins_nuc
        mRNA_state.full_len_transcripts_nuc = full_len_transcripts_nuc
        mRNA_state.single_splice_transcript_nuc = single_splice_transcript_nuc
        mRNA_state.multi_splice_transcript_nuc = multi_splice_transcript_nuc
        
        #update state to new values
        #self.state.set_state('proteins', protein_state)
        #self.state.set_state('mRNAs', mRNA_state)

