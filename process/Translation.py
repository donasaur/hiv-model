# -*- coding: utf-8 -*-
"""
This is a process class whose evolve_state method is called
at each timestep.

An instance of Translation is initialized once per Simulation
with the State as input. The states 'protein' and 'mRNA' are modified
in this process. In particular, only quantities in the cytoplasm are
considered.

At each timestep, for a mRNA strand x, for all x, evolve_state reads
the amount of x and the amount of protein y that x codes for and
writes back new values the new quantity of y in the cytoplasm
depending on how much more of y is generated from the mRNA strand x
during that timestep.

If a represents the quantity of x and b represents the quantity of y,
then during a particular timestep:
    b = b + summation(ci)
where ci is a Poisson random variable with lambda = 4.5, and the
summation ranges from i=1 to i=a.

Note:
Single-spliced mRNA codes for Env, Vif, Vpr, Tat
Multi-spliced mRNA codes for Vif, Tat, Nef, Vpr, Rev
Full-length mRNA codes for Gag or Gag/Pro/Pol

One of the preconditions of evolve_state is that there are no
<MRNAExport.NUM_OF_REV_REQ_FOR_EXPORT Rev-bound single-spliced/full-length
mRNAs in the cytoplasm.

Summary of the biology:
This process takes into account the translation of proteins of
each particular type of mRNA strand in the cytoplasm.

If a full length transcript is ribeosome bound,
there is a 1/20 chance of a frameshift at the STOP site for Gag,
such that the STOP is passed over and the full length Gag/Pro/Pol is 
translated. --Coffin et al.

"""

import numpy as np
from mainaux.State import State
from mainaux.Process import Process
from process.RevBinding import RevBinding
from state.Proteins import Proteins
from mainaux.InitParamValues import *

#References:
#1. Kim, H., Yin, J. (2005) In silico mutagenesis of RNA Splicing in HIV-1. Biotechnology and bioengineering 91: 877-893.
#2. Coffin, J.M., Hughes, S.H., Varmus, H.E. (1997) Retroviruses. Editor: Cold Spring Harbor Laboratory Press, Cold Spring Harbor (NY); ISBN-10: 0-87969-571-4

#This is a Process Class
class Translation(Process):
    def __init__(self, state, param_dict=None):
        self.state = state

        if param_dict==None:
            param_dict = generate_param_dict();
                
        #Constant parameters
        self.FREQ_TRANSLATION = param_dict['FREQ_TRANSLATION'] #proteins/min #kim/yin 2005
        self.FREQ_TRANSLATION_SUPPRESSED = param_dict['FREQ_TRANSLATION_SUPPRESSED'] #proteins/min #fittable parameter, not found in literature
        self.FREQ_TRANSLATION_IRES = param_dict['FREQ_TRANSLATION_IRES'] #proteins/min #fittable parameter, not found in literature
        self.MAX_REV_PER_TRANSCRIPT = param_dict['MAX_REV_PER_TRANSCRIPT']
        self.FREQ_GAG_PRO_POL_TRANSLATION = param_dict['FREQ_GAG_PRO_POL_TRANSLATION'] #fraction of time a full length transcript with be translated to Gag/Pro/Pol rather than Gag #Coffin et al. 1997

    def evolve_state(self, timestep):
        #Rev Binding Timescale:
        #The rev binding and dissociating rates are very fast, and I am unable to see reasonable dynamics 
        #at a 1 min timestep becuase of the multiple events that should be occuring withing the timestep
        #Therefore "zooming in time" to a 1 second timestep for Rev binding--making DETERMINISTIC
        #is there a way to add stochasticity?
        #How to discretize??
    
        #get variables
        protein_state = self.state.get_state('proteins')
        proteins_cyt = protein_state.proteins_cyt
        mRNA_state = self.state.get_state('mRNAs')
        single_splice_transcript_cyt = mRNA_state.single_splice_transcript_cyt
        multi_splice_transcript_cyt = mRNA_state.multi_splice_transcript_cyt
        full_len_transcripts_cyt = mRNA_state.full_len_transcripts_cyt
        reaction_rates_state = self.state.get_state('reaction_rates')
        translation_suppressed = reaction_rates_state.translation_suppressed        
        
        index = Proteins.index
        
        #evolve state
        
        #Translation proceeds at its maximal rate until the cell gets to cell cycle G2/M arrest
        if translation_suppressed == 1:
            translation_frequency = self.FREQ_TRANSLATION_SUPPRESSED
        else:
            translation_frequency = self.FREQ_TRANSLATION
        
        #assuming abundance of translation machinery
        for j in range(7): #Single splice
            for i in np.arange(0+j,(7*self.MAX_REV_PER_TRANSCRIPT)+1+j,7):
                tempRand = np.sum(np.random.poisson(translation_frequency, single_splice_transcript_cyt[i])) #determine num of translation events per transcript
                if j>=3:
                    proteins_cyt[index['Env']]=proteins_cyt[index['Env']]+tempRand #Env
                else:
                    proteins_cyt[j]=proteins_cyt[j]+tempRand #Vif, Vpr, Tat

        for j in range(np.size(multi_splice_transcript_cyt)): #multi splice
            tempRand = np.sum(np.random.poisson(translation_frequency, multi_splice_transcript_cyt[j])) #determine num of translation events per transcript
            if j==0:
                proteins_cyt[index['Vif']]=proteins_cyt[index['Vif']]+tempRand #Vif
            elif j==1 or j==7 or j==12:
                proteins_cyt[index['Tat']]=proteins_cyt[index['Tat']]+tempRand #Tat
            elif j==5 or j==11 or j==16:
                proteins_cyt[index['Nef']]=proteins_cyt[index['Nef']]+tempRand #Nef
            elif j==6:
                proteins_cyt[index['Vpr']]=proteins_cyt[index['Vpr']]+tempRand #Vpr
            else:
                proteins_cyt[index['Rev']]=proteins_cyt[index['Rev']]+tempRand #Rev

        for j in range(np.size(full_len_transcripts_cyt)): #full
            tempRand = np.sum(np.random.poisson(translation_frequency, full_len_transcripts_cyt[j]))
            #if regular 5'cap based translation is suppressed, then HIV-1 IRES can ensure synthesis of Gag/Pol. 
            #Unsure of the rate at which this occurs. 
            if translation_suppressed == 1:            
                tempRand = np.sum(np.random.poisson(self.FREQ_TRANSLATION_IRES, full_len_transcripts_cyt[j]))
            #If a full length transcript is ribeosome bound,
            #there is a self.FREQ_GAG_PRO_POL_TRANSLATION chance of a frameshift at the STOP site for Gag,
            #such that the STOP is passed over and the full length Gag/Pro/Pol is translated #Coffin et al. 
            tempRand2 = (np.random.rand(tempRand)<self.FREQ_GAG_PRO_POL_TRANSLATION).sum()
            proteins_cyt[index['Gag']]=proteins_cyt[index['Gag']]+tempRand-tempRand2 #Gag
            proteins_cyt[index['GagProPol']]=proteins_cyt[index['GagProPol']]+tempRand2 #Gag
                    
        #write back parameters to state object
        protein_state.protein_cyt = proteins_cyt
        mRNA_state.single_splice_transcript_cyt = single_splice_transcript_cyt
        mRNA_state.multi_splice_transcript_cyt = multi_splice_transcript_cyt 
        mRNA_state.full_len_transcripts_cyt = full_len_transcripts_cyt
        
        #update state to new values
        #self.state.set_state('proteins', protein_state)
        #self.state.set_state('mRNAs', mRNA_state)

