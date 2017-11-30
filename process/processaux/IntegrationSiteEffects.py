# -*- coding: utf-8 -*-
"""
This is an initialization class whose initialize_state method is called
before transcription begins.

An instance of IntegrationSiteEffects is initialized once per Simulation by Transcription.

This class sets the RATE_PROMOTER_ON and RATE_PROMOTER_OFF and BASAL_TRANSCRIPTION_RATE
at the start of the simulation. These values then remian constant for the rest of the 
simulation. 

Summary of Biology:
It has been discovered that the site of HIV intregration in the host genome
affects the activity of the HIV promoter. In a simple model posed by Skupsky et
al., the HIV promoter can be in an ON or OFF rate. The transition between these
states is controlled by two parameters: RATE_PROMOTER_ON and RATE_PROMOTER_OFF.
If the promoter is in an ON state, and Tat feedback has not yet started, then
transcription proceeds at the BASAL_TRANSCRIPTION_RATE. If the promoter is in 
an OFF state, and Tat feedback has not yet started, then transcription proceeds
at rate = 0. For now, we assume that after Tat feedback starts, the promoter is
constitutively in an ON state, and the rate of transcription is determined
by the Tat levels in the nucleus. 
The RATE_PROMOTER_ON is very low, resulting in rare bursts of viral transcription.
Skupsky et al. measured the gene expression levels, transcription burst size, 
and transcription burst frequency in 30 clones with different integration sites.
They computationally fit the values of RATE_PROMOTER_ON and RATE_PROMOTER_OFF
and BASAL_TRANSCRIPTION_RATE for each clone. We will usedt hese derived values 
for the constants in our simulation.
We make the major assumption that these 30 clones represent the variation of
these parameters across all integration sites in correct proportion. For each 
individual simulation, we will randomly pick a number 1-30, and that simulation
will run using the Skupsky et al. constants for the corresponding clone.  

Clonal data was obtained from Skupsky et al. Figure 4 A and B.
Step 1: Correlate clones between A and B based on log10(mu) mean
Step 2: Take 10^ of values in graphs to get b, ka/kt-
Step 3: Skupsky et al set the transcript degradation rate to be 0.0033/min.
This is not too different from the Kim/Yin et al 0.0029/min. kt- = 0.0033
as this was used for the Skupsky et al. fitting. 
Step 4: Skupsky et al. set the RATE_PROMOTER_OFF = kR = 20*kt- = 0.066
Step 5: RATE_PROMOTER_ON = ka = (ka/kt-)*0.0033
Step 6: BASAL_TRANSCRIPTION_RATE = kt+ = b*kR

References:
#1. Skupsky, R., Burnett, J. C., Foley, J. E., Schaffer, D. V, & Arkin, A. P.
(2010). HIV promoter integration site primarily modulates transcriptional burst
size rather than frequency. PLoS computational biology, 6(9).
doi:10.1371/journal.pcbi.1000952

"""

import numpy as np

#This is an initialization class that Transcription will call
class IntegrationSiteEffects(object):
    def __init__(self):
        #Constant parameters
        #Skupsky et al. 
        #1/min                     
        self.CLONES_RATE_PROMOTER_ON = np.array([0.002082159, 0.010435516, 0.000791615, 0.003702661, 0.001898952, 0.002446324, 0.000828923, 
                                                 0.001855726, 0.003151476, 0.006005013, 0.004350247, 0.007916149, 0.005351973, 0.007735955,
                                                 0.002744820, 0.001813485, 0.002808756, 0.002621283, 0.004059887, 0.006584366, 0.007735955,
                                                 0.003455524, 0.006434487, 0.005734743, 0.003536014, 0.002130659, 0.006005013, 0.005868322, 
                                                 0.006584366, 0.004451578])
        self.CLONES_BASAL_TRANSCRIPTION_RATE = np.array([3.158958609, 0.104602951, 0.371145275, 0.117366441, 0.250925002, 0.199316814,
                                                         0.630295107, 0.301678205, 0.288100449, 0.165784504, 0.281542482, 0.239631516,
                                                         0.406952701, 0.338488513, 0.870049447, 1.316873128, 1.046029507, 1.201002567,
                                                         0.830890772, 0.512323097, 0.489264759, 1.201002567, 0.675373375, 0.757781390,
                                                         1.228977510, 1.903460792, 0.793494527, 0.811977389, 0.793494527, 1.286897436])
        self.CLONES_RATE_PROMOTER_OFF = 0.066 
    
    def initialize_constants(self):
        clone = np.random.randint(0,30) #equal probability of selecting any of the 30 clones
        self.RATE_PROMOTER_ON = self.CLONES_RATE_PROMOTER_ON[clone]
        self.RATE_PROMOTER_OFF = self.CLONES_RATE_PROMOTER_OFF
        self.BASAL_TRANSCRIPTION_RATE = self.CLONES_BASAL_TRANSCRIPTION_RATE[clone]
        return [self.RATE_PROMOTER_ON, self.RATE_PROMOTER_OFF , self.BASAL_TRANSCRIPTION_RATE]
