# -*- coding: utf-8 -*-
"""
This is a process class whose evolve_state method is called at each timestep.

An instance of Trasncription is initialized once per Simulation with the State
as input. The states 'mRNAs' and 'DNAs' are modified in this process.

At each timestep, evolve_state reads in the abundance of the abundance 0-Rev bound
Gag/Pol and writes to it a new abundance value, depending on whether any mRNA transcript
is created during this timestep. The state of the promoter may also change during this
timestep, in which case its new state value is written to the record as well.

Things we consider during Transcription are the following:
-Integration site effects (what part of the host the virus integrates itself with)
-Actual rate = max(Tat derived transcription rate, basal rate)
-Actual rate upper limit of transcription
-When actual rate > threshold rate, keep promoter on indefinitely
-State of the promoter (always on after actual rate > threshold rate)
-Different method of transcribing mRNA depending on magnitude of rate (random vs. Poisson)

Summary of the biology:
This process takes into account the process of mRNA transcription (specifically, 0-Rev
bound Gag/Pol). This performance of this transcription process depends heavily on the site
of HIV integration in the host genome. The transcription process is dependent on the state
of a promoter region.

"""

import numpy as np
from mainaux.State import State
from mainaux.Process import Process
from process.processaux.IntegrationSiteEffects import IntegrationSiteEffects
from mainaux.InitParamValues import *

#References:
#1. Kim, H., Yin, J. (2005) In silico mutagenesis of RNA Splicing in HIV-1. Biotechnology and bioengineering 91: 877-893.

#This is a Process Class
class Transcription(Process):
    def __init__(self, state, param_dict=None):
        self.state = state

        if param_dict==None:
            param_dict = generate_param_dict();
                
        #User Option
        self.turn_integration_site_effects_on = 0 #0 for off, 1 for on   
        #TODO: move the toggle to this option in simulation or runSimulation, where easier for user to specify

        #Constant parameters
        self.MAX_TAT_ENHANCEMENT = param_dict['MAX_TAT_ENHANCEMENT'] #100x the Basal Rate #Kim et al.        
        self.THRESH_TAT_FEEDBACK = param_dict['THRESH_TAT_FEEDBACK'] #Threshold of Tat feedback for constitutive (always ON) promoter activity; #Arbitrary value based on Kim/Yin Basal rate 
                
        #Default Constant parameters
        self.BASAL_TRANSCRIPTION_RATE = param_dict['BASAL_TRANSCRIPTION_RATE'] #Skupsky et al.
        self.PROMOTER_ON_RATE = param_dict['PROMOTER_ON_RATE'] #0.0044
        self.PROMOTER_OFF_RATE =  param_dict['PROMOTER_OFF_RATE'] #0.066
        
        if self.turn_integration_site_effects_on == 1:
            self.IntegrationSiteEffects = IntegrationSiteEffects()
            [self.PROMOTER_ON_RATE, self.PROMOTER_OFF_RATE, self.BASAL_TRANSCRIPTION_RATE] = self.IntegrationSiteEffects.initialize_constants()

        self.PROMOTER_ALWAYS_ON = False
        self.ACTUAL_TRANSCRIPTION_RATE = self.BASAL_TRANSCRIPTION_RATE
            
    # def evolve_state(self, timestep):
        
    #     #Test, delete this
    #     print(self.PROMOTER_ON_RATE, self.PROMOTER_OFF_RATE, self.BASAL_TRANSCRIPTION_RATE)
    
    #     #get variables
    #     DNA_state = self.state.get_state('DNAs')
    #     promoter_activity = DNA_state.promoter_activity #mRNA per min
    #     reaction_rate_state = self.state.get_state('reaction_rates')
    #     Tat_derived_transcription_rate = reaction_rate_state.Tat_derived_transcription_rate # num of mRNA created per min (on avg)
    #     mRNA_state = self.state.get_state('mRNAs')
    #     full_len_transcripts_nuc = mRNA_state.full_len_transcripts_nuc
    #     transcripts_synthesized = mRNA_state.transcripts_synthesized
        
    #     #evolve state
    #     #Case 1: Promoter is off, Tat feedback is less than threshold level for constitutive activity
    #     if promoter_activity == 0 and Tat_derived_transcription_rate < self.THRESH_TAT_FEEDBACK:
    #         #No transcription can occur
    #         #Will the promoter turn on?
    #         if np.random.rand() < self.PROMOTER_ON_RATE:
    #             promoter_activity = 1
    #     #Case 2: Promoter is on, Tat feedback is less than threshold level for constitutive activity
    #     #Case 3: Promoter is on, Tat feedback level is greater than threshold for constitutive activity, but the Tat derived rate is less than the Basal Rate
    #     elif (promoter_activity == 1 and Tat_derived_transcription_rate < self.THRESH_TAT_FEEDBACK) or (promoter_activity == 1 and Tat_derived_transcription_rate > self.THRESH_TAT_FEEDBACK and Tat_derived_transcription_rate < self.BASAL_TRANSCRIPTION_RATE):
    #         #Transcription occurs at the Basal rate
    #         if self.BASAL_TRANSCRIPTION_RATE <1:
    #             if np.random.rand() < self.BASAL_TRANSCRIPTION_RATE:         #if rate less than zero use as a probability
    #                 full_len_transcripts_nuc[0] +=1
    #                 transcripts_synthesized +=1  
    #         else:
    #             tempValue = np.random.poisson(self.BASAL_TRANSCRIPTION_RATE)
    #             full_len_transcripts_nuc[0] = full_len_transcripts_nuc[0] + tempValue #Arbitrary Poisson noise here
    #             transcripts_synthesized = transcripts_synthesized + tempValue
    #         #Will the promoter turn off?
    #         if np.random.rand() < self.PROMOTER_OFF_RATE:
    #             promoter_activity = 0
    #     #Case 4: Promoter is off, Tat feedback level is greater than threshold for constitutive activity, but the Tat derived rate is less than the Basal Rate if promoter were on 
    #         #Use the Tat feedback rate, but still see if promoter will turn on becuase if it did, transcription can occur faster
    #     elif promoter_activity == 0 and Tat_derived_transcription_rate > self.THRESH_TAT_FEEDBACK and Tat_derived_transcription_rate < self.BASAL_TRANSCRIPTION_RATE:
    #         #Transcription occurs at the Tat derived rate
    #         if Tat_derived_transcription_rate < 1:
    #             if np.random.rand() < Tat_derived_transcription_rate:         #if rate less than zero use as a probability
    #                 full_len_transcripts_nuc[0] +=1
    #                 transcripts_synthesized +=1  
    #             else:
    #                 tempValue = np.random.poisson(Tat_derived_transcription_rate)
    #                 full_len_transcripts_nuc[0] = full_len_transcripts_nuc[0] + tempValue #Arbitrary Poisson noise here
    #                 transcripts_synthesized = transcripts_synthesized + tempValue
    #         #Will the promoter turn on?
    #         if np.random.rand() < self.PROMOTER_ON_RATE:
    #             promoter_activity = 1
    #     #Case 5: Tat feedback level is greater than threshold for constitutive activity and greater than the basal rate
    #     elif Tat_derived_transcription_rate > self.THRESH_TAT_FEEDBACK and Tat_derived_transcription_rate < self.THRESH_TAT_FEEDBACK * self.MAX_TAT_ENHANCEMENT:
    #         tempValue = np.random.poisson(Tat_derived_transcription_rate)
    #         full_len_transcripts_nuc[0] = full_len_transcripts_nuc[0] + tempValue #Arbitrary Poisson noise here
    #         transcripts_synthesized = transcripts_synthesized + tempValue
    #         promoter_activity = 1
    #     else: #Also don't allow transcription to process faster than the upper limit (essentially limit on transcription machinery) 
    #         tempValue = np.random.poisson(self.BASAL_TRANSCRIPTION_RATE*self.MAX_TAT_ENHANCEMENT)
    #         full_len_transcripts_nuc[0] = full_len_transcripts_nuc[0] + tempValue
    #         transcripts_synthesized = transcripts_synthesized + tempValue

    #     #write back parameters to state object
    #     mRNA_state.full_len_transcripts_nuc = full_len_transcripts_nuc
    #     mRNA_state.transcripts_synthesized = transcripts_synthesized
    #     DNA_state.promoter_activity = promoter_activity

    def evolve_state(self, timestep):
        #get variables
        DNA_state = self.state.get_state('DNAs')
        promoter_activity = DNA_state.promoter_activity #mRNA per min
        reaction_rate_state = self.state.get_state('reaction_rates')
        Tat_derived_transcription_rate = reaction_rate_state.Tat_derived_transcription_rate # num of mRNA created per min (on avg)
        mRNA_state = self.state.get_state('mRNAs')
        full_len_transcripts_nuc = mRNA_state.full_len_transcripts_nuc
        transcripts_synthesized = mRNA_state.transcripts_synthesized

        # Deal with what happens if derived rate is greater than threshold rate
        if Tat_derived_transcription_rate > self.THRESH_TAT_FEEDBACK and self.PROMOTER_ALWAYS_ON == False:
            # print('ON')
            # print(Tat_derived_transcription_rate)
            # print(self.THRESH_TAT_FEEDBACK)
            self.PROMOTER_ALWAYS_ON = True
            promoter_activity = 1

        # Case where Tat_derived_transcription_rate <= self.THRESH_TAT_FEEDBACK (did not reach threshold rate yet)
        # Assign the value of promoter_activity here
        if self.PROMOTER_ALWAYS_ON == False:
            if promoter_activity == 0:
                if np.random.rand() < self.PROMOTER_ON_RATE:
                    promoter_activity = 1
            else: # promoter_activity == 1 at the beginning of this timestep
                if np.random.rand() < self.PROMOTER_OFF_RATE:
                    promoter_activity = 0            

        # Deal with what happens if derived rate is greater than basal rate
        if Tat_derived_transcription_rate > self.BASAL_TRANSCRIPTION_RATE:
            self.ACTUAL_TRANSCRIPTION_RATE = Tat_derived_transcription_rate
        else:
            self.ACTUAL_TRANSCRIPTION_RATE = self.BASAL_TRANSCRIPTION_RATE

        # Deal with what happens if transcription rate is greater than upper limit
        # The difference in the factor being multiplied with self.MAX_TAT_ENHANCEMENT is for consistency
        # with own evolve_state method
        if self.ACTUAL_TRANSCRIPTION_RATE > self.THRESH_TAT_FEEDBACK * self.MAX_TAT_ENHANCEMENT:
            self.ACTUAL_TRANSCRIPTION_RATE = self.THRESH_TAT_FEEDBACK *self.MAX_TAT_ENHANCEMENT

        # Deal with the actual creation of mRNA
        if promoter_activity == 1 or self.PROMOTER_ALWAYS_ON:
            if self.ACTUAL_TRANSCRIPTION_RATE < 1:
                if np.random.rand() < self.ACTUAL_TRANSCRIPTION_RATE:
                    full_len_transcripts_nuc[0] += 1
                    transcripts_synthesized += 1
            else:
                tempValue = np.random.poisson(self.ACTUAL_TRANSCRIPTION_RATE)
                full_len_transcripts_nuc[0] = full_len_transcripts_nuc[0] + tempValue
                transcripts_synthesized = transcripts_synthesized + tempValue

        # write back parameters to state object
        DNA_state.promoter_activity = promoter_activity
        mRNA_state.transcripts_synthesized = transcripts_synthesized
        
        #update state to new values
        #self.state.set_state('mRNAs', mRNA_state)
