# -*- coding: utf-8 -*-
"""
This is a process class whose evolve_state method is called
at each timestep.

An instance of MRNAExport is initialized once per Simulation
with the State as input. Only the state 'mRNAs' is modified in
this process.

At each timestep, for an mRNA strand x, for all x, evolve_state reads
the amount of x in the nucleus and cytoplasm and writes back new values
of x in the nucleus and cytoplasm depending on how much of the particular
strand is set to be exported from the nucleus.

All types of:
-multi-spliced mRNAs (17)
-8-Rev-bound single-spliced mRNAs (7)
-8-Rev-bound full-length mRNAs (1)
have a chance to be exported.

Summary of the biology:
This process takes into account the export of mRNA strands from the nucleus to 
the cytoplasm at each timestep.

The HIV mRNA transcript has an RRE (Rev Response Element) located between splice
donor 4 (D4) and splice acceptor 7 (A7). This RRE must bind a threshold number of
Rev molecules to enable mRNA export to the cytoplasm. Multi-spliced mRNA transcripts
have excised this RRE element (D4-A7 splice), and therefore may be exported 
immediatly to the cytoplasm. Multi-spliced mRNAs are therefore called 
"Rev-independent" mRNAS. The full length and single-spliced mRNAs still have 
their RRE elements, and therefore are "Rev-dependent". It is unclear from the 
literature, how many Rev molecules must be bound for export. We know that the
value is >1 and <= the max number of bound Rev. The requirement for Rev
binding essentially poses a delay on certain transcripts getting to the 
cytoplasm. The Rev-independent transcripts start translation before the 
Rev-dependent transcripts.  

#References:
#1. Pond, S. J. K., Ridgeway, W. K., Robertson, R., Wang, J., & Millar, D. P. (2009). HIV-1 Rev protein assembles on viral RNA one molecule at a time. Proceedings of the National Academy of Sciences of the United States of America, 106(5), 1404–8. doi:10.1073/pnas.0807388106
#2. Kim, H., Yin, J. (2005) In silico mutagenesis of RNA Splicing in HIV-1. Biotechnology and bioengineering 91: 877-893.

"""

import numpy as np
from mainaux.State import State
from mainaux.Process import Process
from process.RevBinding import RevBinding
from mainaux.InitParamValues import *

#This is a Process Class
class MRNAExport(Process):
    def __init__(self, state, param_dict=None):
        self.state = state

        if param_dict==None:
            param_dict = generate_param_dict();    
        
        #Constant parameters
        self.NUM_OF_REV_REQ_FOR_EXPORT = param_dict['NUM_OF_REV_REQ_FOR_EXPORT'] #fittable value. Pond et al. state it is >1.
        self.PROB_REV_INDEP_EXPORT = param_dict['PROB_REV_INDEP_EXPORT'] #1/min #Kim, H., Yin, J. (2005) 
        self.PROB_REV_DEP_EXPORT = param_dict['PROB_REV_DEP_EXPORT'] #1/min #Kim, H., Yin, J. (2005) 
        
        self.MAX_REV_PER_TRANSCRIPT = param_dict['MAX_REV_PER_TRANSCRIPT']
        
    def nuclear_export(self, what_may_be_exported, abundances_nuc, abundances_cyt, export_rate):
        #what_may_be_exported = array/list of indexes in abundances_nuc/Cyt of constructs to export
        #abundances_nuc = starting abundances of things to be exported
        #abundances_cyt = starting abundances of things in the destination location
        #export_rate = rate of export
        for i in what_may_be_exported:
            temp_rand = np.random.rand(abundances_nuc[i])
            #count how many are less than the export probability, and do accounting
            decrement_amount = (temp_rand<export_rate).sum()
            abundances_cyt[i]=abundances_cyt[i]+decrement_amount
            abundances_nuc[i]=abundances_nuc[i]-decrement_amount
        return [abundances_nuc, abundances_cyt]
        
    def evolve_state(self, timestep):  
        
        #get variables
        mRNA_state = self.state.get_state('mRNAs')
        full_len_transcripts_nuc = mRNA_state.full_len_transcripts_nuc
        full_len_transcripts_cyt = mRNA_state.full_len_transcripts_cyt
        single_splice_transcript_nuc = mRNA_state.single_splice_transcript_nuc
        single_splice_transcript_cyt = mRNA_state.single_splice_transcript_cyt
        multi_splice_transcript_nuc = mRNA_state.multi_splice_transcript_nuc
        multi_splice_transcript_cyt = mRNA_state.multi_splice_transcript_cyt
        
        #Part A. Rev Independent Export
        #The fully spliced transcripts stored in multi_splice_transcript_nuc do not have a RRE (Rev binding) element
        #and therefore can be exported without rev
        what_may_be_exported = range(np.size(multi_splice_transcript_nuc))
        [multi_splice_transcript_nuc, multi_splice_transcript_cyt] = self.nuclear_export(what_may_be_exported, multi_splice_transcript_nuc, multi_splice_transcript_cyt, self.PROB_REV_INDEP_EXPORT)
        
        #Part B. Rev Dependent Export
        #The unspliced and singly spliced transcripts stored in full/single_splice_transcript_nuc have a RRE (Rev binding) element
        #and therefore cannot be exported without rev 
        #start with full transcripts
        what_may_be_exported = np.arange(self.NUM_OF_REV_REQ_FOR_EXPORT, self.MAX_REV_PER_TRANSCRIPT+1) # Contains only one value right now
        [full_len_transcripts_nuc, full_len_transcripts_cyt] = self.nuclear_export(what_may_be_exported, full_len_transcripts_nuc, full_len_transcripts_cyt, self.PROB_REV_DEP_EXPORT)
        #single splice transcripts
        what_may_be_exported = np.arange((self.NUM_OF_REV_REQ_FOR_EXPORT*7),(7*(self.MAX_REV_PER_TRANSCRIPT+1))) # Contains indices of 8-Rev elements
        what_may_be_exported = what_may_be_exported[np.where(single_splice_transcript_nuc[what_may_be_exported]>0)[0]]
        [single_splice_transcript_nuc, single_splice_transcript_cyt] = self.nuclear_export(what_may_be_exported, single_splice_transcript_nuc, single_splice_transcript_cyt, self.PROB_REV_DEP_EXPORT)
        
        #write back parameters to state object
        mRNA_state.full_len_transcripts_nuc = full_len_transcripts_nuc
        mRNA_state.single_splice_transcript_nuc = single_splice_transcript_nuc
        mRNA_state.full_len_transcripts_cyt = full_len_transcripts_cyt
        mRNA_state.single_splice_transcript_cyt = single_splice_transcript_cyt
        mRNA_state.multi_splice_transcript_nuc = multi_splice_transcript_nuc 
        mRNA_state.multi_splice_transcript_cyt = multi_splice_transcript_cyt
        
        #update state to new values
        #self.state.set_state('mRNAs', mRNA_state)
