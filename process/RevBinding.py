# -*- coding: utf-8 -*-
"""
This is a process class whose evolve_state method is called at
each timestep.

An instance of RevBinding is intialized once per Simulation with
the State as input. The states 'mRNAs' and 'proteins' are modified
in this process.

At each timestep, evolve_state reads in the abundances of full-length
and singly-spliced mRNAs, as well as the amount of free Rev in the
nucleus.

The function partitions the amount of free Rev between the full-length
mRNA and singly-spliced mRNAs; the Rev visible to each type of mRNA is
known as bindable Rev.

The abundances of different mRNA types in the previous timestep +
the bindable Rev count is used as the initial condition to calculate
for an ODE solution that gives the new abundances of each mRNA type
as well as the new abundance of free Rev in the nucleus.

This ODE solution is calculated over 60 timesteps, where each timestep
spans 1 second.

The new mRNA abundances plus the free Rev count is then written back
to the corresponding State classes.

Summary of the biology:
This process takes into account the process of Rev binding and
unbinding from the full-length and singly-spliced mRNAs.

"""

import numpy as np
from scipy.integrate import odeint
from mainaux.Process import Process
from state.Proteins import Proteins
from mainaux.InitParamValues import *

# References:
#1. Pond, S. J. K., Ridgeway, W. K., Robertson, R., Wang, J., & Millar, D. P. (2009). HIV-1 Rev protein assembles on viral RNA one molecule at a time. Proceedings of the National Academy of Sciences of the United States of America, 106(5), 1404–8. doi:10.1073/pnas.0807388106
#2. Kim, H., Yin, J. (2005) In silico mutagenesis of RNA Splicing in HIV-1. Biotechnology and bioengineering 91: 877-893.
#3. Lamond AI, Sleeman JE. Nuclear substructure and dynamics. Curr Biol. 2003 Oct 28 13(21):R825-8.

#This is a Process Class
class RevBinding(Process):
    #Define static variables
    #MAX_REV_PER_TRANSCRIPT = 8; #Pond et al., 2009;  ###MAX_REV_PER_TRANSCRIPT = 12 #Kim and Yin 2005
    def __init__(self, state, param_dict=None):
        self.state = state

        if param_dict==None:
            param_dict = generate_param_dict();
                
        #Constant parameters
        #self.MAX_REV_PER_TRANSCRIPT = 8 #Pond et al., 2009;  ###MAX_REV_PER_TRANSCRIPT = 12 #Kim and Yin 2005
        self.VOLUME_NUC = param_dict['VOLUME_NUC'] #L #Kim et al 2005; ###VOLUME_NUC = 5 * (10**-13) #L #fibroblast #Lamond AI, Sleeman JE. 2003 
        self.REV_BINDING_CONSTANTS = param_dict['REV_BINDING_CONSTANTS'] #Pond et al., 2009; last 4 values not-reported-taken as average of the 1st 4
        self.REV_BINDING_CONSTANTS_SCALED = (self.REV_BINDING_CONSTANTS*(10**8))/(self.VOLUME_NUC * 6.022 *(10**23)) #1/(molecule of Rev * sec)
        self.REV_DISSOCIATION_CONSTANTS = param_dict['REV_DISSOCIATION_CONSTANTS'] #Pond et al., 2009; last 4 values not-reported-taken as average of the 1st 4 #1/sec

        self.MAX_REV_PER_TRANSCRIPT = param_dict['MAX_REV_PER_TRANSCRIPT']

    def Rev_ode(self, R, t):
        #R is vector of length MAX_REV_PER_TRANSCRIPT + 1 + free Rev concatenated at the end
        #R[0] = unbound mRNA
        #R[1] = mRNA - 1Rev
        #R[2] = mRNA - 2Rev
        #...
        #R[8] = mRNA - 8Rev
        #R[9] = free Rev
        f0 = R[1]*self.REV_DISSOCIATION_CONSTANTS[0]-R[0]*self.REV_BINDING_CONSTANTS_SCALED[0]*R[9]
        f1 = R[0]*self.REV_BINDING_CONSTANTS_SCALED[0]*R[9]-R[1]*self.REV_DISSOCIATION_CONSTANTS[0]+R[2]*self.REV_DISSOCIATION_CONSTANTS[1]-R[1]*self.REV_BINDING_CONSTANTS_SCALED[1]*R[9]
        f2 = R[1]*self.REV_BINDING_CONSTANTS_SCALED[1]*R[9]-R[2]*self.REV_DISSOCIATION_CONSTANTS[1]+R[3]*self.REV_DISSOCIATION_CONSTANTS[2]-R[2]*self.REV_BINDING_CONSTANTS_SCALED[2]*R[9]
        f3 = R[2]*self.REV_BINDING_CONSTANTS_SCALED[2]*R[9]-R[3]*self.REV_DISSOCIATION_CONSTANTS[2]+R[4]*self.REV_DISSOCIATION_CONSTANTS[3]-R[3]*self.REV_BINDING_CONSTANTS_SCALED[3]*R[9]
        f4 = R[3]*self.REV_BINDING_CONSTANTS_SCALED[3]*R[9]-R[4]*self.REV_DISSOCIATION_CONSTANTS[3]+R[5]*self.REV_DISSOCIATION_CONSTANTS[4]-R[4]*self.REV_BINDING_CONSTANTS_SCALED[4]*R[9]
        f5 = R[4]*self.REV_BINDING_CONSTANTS_SCALED[4]*R[9]-R[5]*self.REV_DISSOCIATION_CONSTANTS[4]+R[6]*self.REV_DISSOCIATION_CONSTANTS[5]-R[5]*self.REV_BINDING_CONSTANTS_SCALED[5]*R[9]
        f6 = R[5]*self.REV_BINDING_CONSTANTS_SCALED[5]*R[9]-R[6]*self.REV_DISSOCIATION_CONSTANTS[5]+R[7]*self.REV_DISSOCIATION_CONSTANTS[6]-R[6]*self.REV_BINDING_CONSTANTS_SCALED[6]*R[9]
        f7 = R[6]*self.REV_BINDING_CONSTANTS_SCALED[6]*R[9]-R[7]*self.REV_DISSOCIATION_CONSTANTS[6]+R[8]*self.REV_DISSOCIATION_CONSTANTS[7]-R[7]*self.REV_BINDING_CONSTANTS_SCALED[7]*R[9]
        f8 = R[7]*self.REV_BINDING_CONSTANTS_SCALED[7]*R[9]-R[8]*self.REV_DISSOCIATION_CONSTANTS[7]
        f9 = (R[1]*self.REV_DISSOCIATION_CONSTANTS[0]-R[0]*self.REV_BINDING_CONSTANTS_SCALED[0]*R[9]+
              R[2]*self.REV_DISSOCIATION_CONSTANTS[1]-R[1]*self.REV_BINDING_CONSTANTS_SCALED[1]*R[9]+
              R[3]*self.REV_DISSOCIATION_CONSTANTS[2]-R[2]*self.REV_BINDING_CONSTANTS_SCALED[2]*R[9]+
              R[4]*self.REV_DISSOCIATION_CONSTANTS[3]-R[3]*self.REV_BINDING_CONSTANTS_SCALED[3]*R[9]+
              R[5]*self.REV_DISSOCIATION_CONSTANTS[4]-R[4]*self.REV_BINDING_CONSTANTS_SCALED[4]*R[9]+
              R[6]*self.REV_DISSOCIATION_CONSTANTS[5]-R[5]*self.REV_BINDING_CONSTANTS_SCALED[5]*R[9]+
              R[7]*self.REV_DISSOCIATION_CONSTANTS[6]-R[6]*self.REV_BINDING_CONSTANTS_SCALED[6]*R[9])
        return [f0, f1, f2, f3, f4, f5, f6, f7, f8, f9]
        
    def ODE_discretizer(self, soln, prev_mRNA_abundances, prev_protein_abundances):
        #This function discretizes, mass balances, and ensures non-negative values of ODE solutions for integration with the rest of the system
        #soln: ODE solution
        #prev_mRNA_abundances: abundances of mRNAs before applying the ODE (last timestep)
        #prev_protein_abundances: abundance of free Rev before applying the ODE (last timestep)
        soln[-1,:][soln[-1,:] == .5] = 1
        soln_round = np.round(soln[-1,:]) # Rounds all Rev balances to the nearest whole number for the current timestep
        #print soln[-1,:] # Should only be printing the array solution
        soln_round[soln_round<0]=0 # Don't allow negative abundances
        mRNA_before = np.sum(prev_mRNA_abundances)
        mRNA_after = np.sum(soln_round[0:-1])
        
        # Mass balance handling (mRNA)
        while mRNA_after != mRNA_before:
            discrepancy = mRNA_after-mRNA_before # positive if mRNA_after > mRNA_before (mRNA was created, so need to remove mRNA from system)
            temp_index = np.random.randint(self.MAX_REV_PER_TRANSCRIPT+1) # Randomly pick bins to adjust the discrepancy; results are from [0, arg).
            soln_round[temp_index]=soln_round[temp_index]-discrepancy # soln_round[temp_index]; want to bring Rev_mRNA_after to the value Rev_mRNA_before
            soln_round[soln_round<0]=0 # This is the reason why the while loop may run multiple times
            mRNA_after = np.sum(soln_round[0:-1])
        
        #Mass Balance (Rev)
        #Compare prev # of Rev in nucleus to current # of Rev in nucleus
        Rev_before = np.sum((prev_mRNA_abundances*np.arange(self.MAX_REV_PER_TRANSCRIPT+1))) + prev_protein_abundances #bound Rev + free Rev @ prev timestep
        Rev_after = np.sum((soln_round[0:-1]*np.arange(self.MAX_REV_PER_TRANSCRIPT+1))) + soln_round[-1] #bound Rev + free Rev at this timestep
        discrepancy = Rev_after - Rev_before # Rev_after > Rev_before means that Rev was created (need to remove Rev from the system)
        #if there is a dicrepancy, deal with it
        if discrepancy == 0: #no issue
            pass
        elif discrepancy < 0: #not enough Rev accounted for in the current timestep
            #this unaccounted for Rev is the result of 1. rounding error and 2. released free Rev in the mRNA adjustment
            # Want to add Rev to the system

            #So, store it back as free Rev
            soln_round[-1] = soln_round[-1] - discrepancy
        else: #discrepancy > 0: #too much Rev accounted for in the current timestep
            #this over accounted for Rev is the result of 1. rounding error and 2. extra Rev bound in the mRNA adjustment
            # Want to remove Rev from the system

            #So, decrement it from free Rev if possible:
            if discrepancy <= soln_round[-1]: #if decrementing from free will not result in a negative value
                soln_round[-1] = soln_round[-1] - discrepancy
            else: #if decrementing from free will result in a negative value
                #First decrement what you can from free
                discrepancy = discrepancy - soln_round[-1]
                soln_round[-1] = 0
                #Then reduce the rest of the discrepancy by moving mRNAs from a higher to lower Rev occupancy
                temp_counter = 0
                while discrepancy != 0:
                    temp_counter += 1
                    temp_index = np.random.randint(1,self.MAX_REV_PER_TRANSCRIPT+1) #randomly pick mRNA bin to adjust the discrepancy; only choose from bins with Rev binding != 0. That is, only pick mRNA bins with at least one Rev bounded.
                    if soln_round[temp_index]>0:
                        soln_round[temp_index] = soln_round[temp_index] - 1
                        soln_round[temp_index-1] = soln_round[temp_index-1] + 1
                        discrepancy = discrepancy - 1
#                    if temp_counter == 99999999: #just to ensure no infinite loop
#                        print("error, proteins_nuc[Proteins.index['Rev']]<0, violation of mass balance")
#                        break
        new_protein_abundances=soln_round[-1]
        new_mRNA_abundances = soln_round[0:-1]
        return [new_mRNA_abundances, new_protein_abundances]

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
        mRNA_state = self.state.get_state('mRNAs')
        full_len_transcripts_nuc = mRNA_state.full_len_transcripts_nuc
        single_splice_transcript_nuc = mRNA_state.single_splice_transcript_nuc
        
        #evolve state
        
        #bindable_Rev = np.min([np.random.poisson(proteins_nuc[Proteins.index['Rev']]), proteins_nuc[Proteins.index['Rev']]]) #arbitrary noise #removed this!
        
        #Allocate Rev amongst different mRNA types relative to mRNA abundances
        #Each type of mRNA will get to potentially bind to a fraction of the total Rev calculated as
        #((Total Rev)/(Total mRNA))*(Abindance of given mRNA)
        #Bindable Rev 
        bindable_mRNA_sum = np.sum([np.sum(full_len_transcripts_nuc[0:8]), np.sum(single_splice_transcript_nuc[0:56])])
        bindable_Rev_sum = 0 #will hold the total amount of Rev that we allow the mRNAs to see
        
        #TODO come up with something to add noise to this allocation
        
        
        #Start with full length transcripts
        # initial ode conditions
        if np.sum(full_len_transcripts_nuc[0:8]) == 0:
            bindable_Rev = 0
        else:
            bindable_Rev = np.floor(np.divide(np.sum(full_len_transcripts_nuc[0:8]), bindable_mRNA_sum)*proteins_nuc[Proteins.index['Rev']])
            # (sum of all the bindable full len mRNAs)/(sum of all the bindable mRNAs)*# of Rev in the nucleus
        bindable_Rev_sum = bindable_Rev_sum + bindable_Rev
        R = np.concatenate((full_len_transcripts_nuc, np.array([bindable_Rev])),0)
        t_seg_Rev  = np.linspace(0, 59, 60)
        soln = odeint(self.Rev_ode, R, t_seg_Rev)
        #discretize, mass balance
        [full_len_transcripts_nuc, net_Rev] = self.ODE_discretizer(soln, full_len_transcripts_nuc, bindable_Rev)
           
        #Next Singly Spliced
        for i in range(7):
            if np.sum(single_splice_transcript_nuc[np.arange(i,50+i,7)]) == 0:
                bindable_Rev = 0
            else:
                bindable_Rev = np.floor(np.divide(np.sum(single_splice_transcript_nuc[np.arange(i,50+i,7)]), bindable_mRNA_sum)*proteins_nuc[Proteins.index['Rev']])
            bindable_Rev_sum = bindable_Rev_sum + bindable_Rev
            # initial ode conditions
            single_spliced_current = single_splice_transcript_nuc[np.arange(i,57+i,7)]
            R = np.concatenate((single_spliced_current, np.array([bindable_Rev])),0)
            t_seg_Rev  = np.linspace(0, 59, 60)
            soln = odeint(self.Rev_ode, R, t_seg_Rev)
            #discretize, mass balance
            [single_spliced_current, temp_Rev] = self.ODE_discretizer(soln, single_spliced_current, bindable_Rev)
            single_splice_transcript_nuc[np.arange(i,57+i,7)] = np.int32(single_spliced_current)
            net_Rev = net_Rev + temp_Rev
            
        proteins_nuc[Proteins.index['Rev']] = (proteins_nuc[Proteins.index['Rev']]-bindable_Rev_sum) + net_Rev  #first term accounts for rounding in allocation. Not all Rev may have been made available for binding               
                    
        #write back parameters to state object
        protein_state.protein_nuc = proteins_nuc
        mRNA_state.full_len_transcripts_nuc = full_len_transcripts_nuc
        mRNA_state.single_splice_transcript_nuc = single_splice_transcript_nuc   
        
        #update state to new values
        #self.state.set_state('proteins', protein_state)
        #self.state.set_state('mRNAs', mRNA_state)
        
    
