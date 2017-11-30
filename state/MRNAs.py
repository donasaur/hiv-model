# -*- coding: utf-8 -*-
import numpy as np

from process.RevBinding import RevBinding
from mainaux.InitParamValues import *

#This is a type of State Class
class MRNAs(object):
    def __init__(self, param_dict=None):
        if param_dict==None:
            param_dict = generate_param_dict();
        
        #Initialize necessary parameters
        MAX_REV_PER_TRANSCRIPT = param_dict['MAX_REV_PER_TRANSCRIPT']
        #Initialize Variables
        self.full_len_transcripts_nuc = np.zeros((MAX_REV_PER_TRANSCRIPT+1), int) #an array of length MAX_REV_PER_TRANSCRIPT+1. Index 0 holds full transcript with no Rev. Index MAX_REV_PER_TRANSCRIPT holds full transcript with MAX_REV_PER_TRANSCRIPT Rev.
        self.single_splice_transcript_nuc = np.zeros((7*(MAX_REV_PER_TRANSCRIPT+1)), int) 
            #an array of length 7*MAX_REV_PER_TRANSCRIPT to hold the counts of each type of single-spliced mRNA: 
            #[0-6] (inclusive) D1 to A1, A2, A3, A4a, A4b, A4c, A5; 0 Rev
            #[7-13]  D1 to A1, A2, A3, A4a, A4b, A4c, A5; 1 Rev
            #[14-20] D1 to A1, A2, A3, A4a, A4b, A4c, A5; 2 Rev
            #[21-27] D1 to A1, A2, A3, A4a, A4b, A4c, A5; 3 Rev
            #[28-34] D1 to A1, A2, A3, A4a, A4b, A4c, A5; 4 Rev
            #[35-41] D1 to A1, A2, A3, A4a, A4b, A4c, A5; 5 Rev
            #[42-48] D1 to A1, A2, A3, A4a, A4b, A4c, A5; 6 Rev
            #[49-55] D1 to A1, A2, A3, A4a, A4b, A4c, A5; 7 Rev
            #[56-62] D1 to A1, A2, A3, A4a, A4b, A4c, A5; 8 Rev
        self.multi_splice_transcript_nuc = np.zeros((17), int) 
            #an array of length X to hold the counts of each type of multi-spliced mRNA: 
            #[0]: D1-A1, D4-A7 #vif
            #[1]: D1-A1, D2-A3, D4-A7 #tat
            #[2]: D1-A1, D2-A4a, D4-A7 #rev
            #[3]: D1-A1, D2-A4b, D4-A7 #rev
            #[4]: D1-A1, D2-A4c, D4-A7 #rev
            #[5]: D1-A1, D2-A5, D4-A7 #nef
            #[6]: D1-A2, D4-A7 #vpr
            #[7]: D1-A2, D3-A3, D4-A7 #tat
            #[8]: D1-A2, D3-A4a, D4-A7 #rev
            #[9]: D1-A2, D3-A4b, D4-A7 #rev
            #[10]: D1-A2, D3-A4c, D4-A7 #rev
            #[11]: D1-A2, D3-A5, D4-A7 #nef
            #[12]: D1-A3, D4-A7 #tat
            #[13]: D1-A4a, D4-D7 #rev
            #[14]: D1-A4b, D4-D7 #rev
            #[15]: D1-A4c, D4-D7 #rev
            #[16]: D1-A5, D4-D7 #nef
        self.full_len_transcripts_cyt = np.zeros((MAX_REV_PER_TRANSCRIPT+1), int)
        self.single_splice_transcript_cyt = np.zeros((7*(MAX_REV_PER_TRANSCRIPT+1)), int) 
        self.multi_splice_transcript_cyt = np.zeros((17), int) 
        self.transcripts_synthesized = 0 #count of how many viral transcripts were made (TOTAL)
        self.full_len_transcripts_Gag_bound = np.zeros((16), int)
            #In practice will be reshaped to a 2x2x2x2 matrix 
            #in which each dimension can be a 0 or a 1 for each of SL1, SL2, SL3, and SL4 being
            #Bound or not bound
            #all of these are in the cytoplasm
            #[0]: nothing bound -- held in self.full_len_transcripts_cyt 
            #[1]: SL4
            #[2]: SL3
            #[3]: SL3, SL4
            #[4]: SL2
            #[5]: SL2, SL4
            #[6]: SL2, SL3
            #[7]: SL2, SL3, SL4
            #[8]: SL1
            #[9]: SL1, SL4
            #[10]: SL1, SL3
            #[11]: SL1, SL3, SL4
            #[12]: SL1, SL2
            #[13]: SL1, SL2, SL4
            #[14]: SL1, SL2, SL3
            #[15]: SL1, SL2, SL3, SL4
        self.full_length_transcript_dimers_cyt = np.zeros((3), int)
            #[0]: dimer in cytoplasm with 6 Gag bound (SL1,2,3 on both RNAs, and neither SL4)
            #[1]: dimer in the cytoplasm with 7 Gag bound (SL1,2,3 on both RNAs, and one SL4)
            #[2]: dimer in the cytoplasm with 8 Gag bound (SL1,2,3 on both RNAs, and both SL4)
            
        #Set Keys
        #DO I NEED THIS??
#        self.full_len_transcripts_nuc_key = 'full_len_transcripts_nuc'
#        self.single_splice_transcript_nuc_key = 'single_splice_transcript_nuc'
#        self.multi_splice_transcript_nuc_key = 'multi_splice_transcript_nuc' 
#        self.full_len_transcripts_cyt_key = 'full_len_transcripts_cyt'
#        self.single_splice_transcript_cyt_key = 'single_splice_transcript_cyt'
#        self.multi_splice_transcript_cyt_key = 'multi_splice_transcript_cyt'
#        self.transcripts_synthesized_key = 'transcripts_synthesized'
        
    def record_state(self, record, timestep, max_timesteps):
        record.add_tracking(timestep, max_timesteps, 'full_len_transcripts_nuc', self.full_len_transcripts_nuc)
        record.add_tracking(timestep, max_timesteps, 'single_splice_transcript_nuc', self.single_splice_transcript_nuc)
        record.add_tracking(timestep, max_timesteps, 'multi_splice_transcript_nuc', self.multi_splice_transcript_nuc)
        record.add_tracking(timestep, max_timesteps, 'full_len_transcripts_cyt', self.full_len_transcripts_cyt)
        record.add_tracking(timestep, max_timesteps, 'single_splice_transcript_cyt', self.single_splice_transcript_cyt)
        record.add_tracking(timestep, max_timesteps, 'multi_splice_transcript_cyt', self.multi_splice_transcript_cyt)
        record.add_tracking(timestep, max_timesteps, 'transcripts_synthesized', self.transcripts_synthesized)
        record.add_tracking(timestep, max_timesteps, 'full_len_transcripts_Gag_bound', self.full_len_transcripts_Gag_bound)
        record.add_tracking(timestep, max_timesteps, 'full_length_transcript_dimers_cyt', self.full_length_transcript_dimers_cyt)

    def record_at_end(self, record):
        pass
        
