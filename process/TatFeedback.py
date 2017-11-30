# -*- coding: utf-8 -*-
"""
This is a process class whose evolve_state method is called
at each timestep.

An instance of TatFeedback is initialized once per Simulation
with the State as input. The states 'host_factors', 'proteins',
and 'reaction_rates' are modified in this process.

At each timestep, evolve_state reads in the amount of Tat in the
nucleus, the amount of pTEFb in the nucleus, the amount of
deacetylated Tat-pTEFb complex, and the amount of acetylated
Tat-pTEFb complex.

The amount of pTEFb in the nucleus during this timestep is given
by the following equation:
    pTEFb = pTEFb(from prev. timestep) + pTefb0*exp(R*(t+1))
                - pTefb0*exp(R*t)
where pTefb0 is the initial pTefb amount and R is the pTEFb
doubling rate.

Tat and pTEFb moves between its complex form and its free form
in this timestep(composed of 60 second steps), and thus the
quantities of interest above vary during this timestep.

At the end of the timestep, the new amounts of Tat in the
nucleus + pTEFb in the nucleus + deacetylated Tat-pTEFb complex
+ acetylated Tat-pTEFb complex are written back to the State class.

The derived Tat description rate is also changed during this time step
and its new value is written back as well.

Summary of the biology:
This process takes into account the effect of Tat and pTEFb in
the nucleus on the transcription rate of mRNA.
"""

import numpy as np
from scipy.integrate import odeint
from mainaux.Process import Process
from mainaux.InitParamValues import *

#References:
#1. Weinberger, L.S., Burnett, J.C., Toettcher, J.E., Arkin, A.P., Schaffer, D.V. (2005). Stochastic Gene Expression in a Lentiviral Positive-Feedback Loop: HIV-1 Tat Fluctuations Drive Phenotypic Diversity. Cell 122: 169-182. 
#2. Kim, H., Yin, J. (2005) In silico mutagenesis of RNA Splicing in HIV-1. Biotechnology and bioengineering 91: 877-893.

#This is a Process Class
class TatFeedback(Process):
    def __init__(self, state, param_dict=None):
        self.state = state

        if param_dict==None:
            param_dict = generate_param_dict();  
        #Constant parameters
        self.pTEFb_DOUBLING_RATE =  param_dict['pTEFb_DOUBLING_RATE']
        self.pTEFb_NUC_INIT = param_dict['pTEFb_NUC_INIT']
        
        self.RATE_TAT_pTEFb_BIND = param_dict['RATE_TAT_pTEFb_BIND']  #1/(molecules*sec) #Weinberger 2005
        self.RATE_TAT_pTEFb_UNBIND = param_dict['RATE_TAT_pTEFb_UNBIND'] #1/sec #Weinberger 2005
        self.RATE_TAT_pTEFb_ACETYL = param_dict['RATE_TAT_pTEFb_ACETYL'] #1/(molecules*sec) #Weinberger 2005
        self.RATE_TAT_pTEFb_DEACETYL = param_dict['RATE_TAT_pTEFb_DEACETYL'] #1/sec #Weinberger 2005
        self.RATE_TAT_ACT_TRANSCRIPTION = param_dict['RATE_TAT_ACT_TRANSCRIPTION'] #/sec #Weinberger 2005
        
    # solve the system dy/dt = f(y, t)
    def TatODE(self, y, t):
        Tat_nuc = y[0]
        pTEFb_nuc = y[1]
        Tat_pTEFb_deacetyl = y[2]
        Tat_pTEFb_acetyl = y[3]
        #"mRNA" = y[4]        
        
        # the model equations
        f0 = self.RATE_TAT_ACT_TRANSCRIPTION*Tat_pTEFb_acetyl - self.RATE_TAT_pTEFb_BIND*Tat_nuc*pTEFb_nuc
        f1 = self.RATE_TAT_ACT_TRANSCRIPTION*Tat_pTEFb_acetyl - self.RATE_TAT_pTEFb_BIND*Tat_nuc*pTEFb_nuc
        f2 = self.RATE_TAT_pTEFb_BIND*Tat_nuc*pTEFb_nuc - self.RATE_TAT_pTEFb_ACETYL*Tat_pTEFb_deacetyl
        f3 = self.RATE_TAT_pTEFb_ACETYL*Tat_pTEFb_deacetyl - self.RATE_TAT_ACT_TRANSCRIPTION*Tat_pTEFb_acetyl
        f4 = self.RATE_TAT_ACT_TRANSCRIPTION*Tat_pTEFb_acetyl
        return [f0, f1, f2, f3, f4]
        
    def ODE_discretizer(self, soln, free_Tat, pTEFb_nuc, Tat_pTEFb_deacetyl, Tat_pTEFb_acetyl):     
        #This function discretizes, mass balances, and ensures positive values of ODE solutions for integration with the rest of the system
        #soln: ODE solution
        #prev_mRNA_abundances: abundances of mRNAs before applying the ODE (last timestep)
        #prev_protein_abundances: abundance of free Rev before applying the ODE (last timestep)
        soln[-1,:][soln[-1,:] == .5] = 1
        soln_round = np.around(soln[-1,:]) #discretize
        soln_round[soln_round<0]=0 #don't allow negatives
        Tat_before = free_Tat + Tat_pTEFb_deacetyl + Tat_pTEFb_acetyl 
        Tat_after = np.sum(soln_round[np.array([0,2,3])])
        temp_counter = 0
        while Tat_after != Tat_before: #mass balance (Tat)
            temp_counter +=1
            discrepancy = Tat_after-Tat_before # positive if Tat_after > Tat_before (Tat was created, so need to remove Tat from system)
            array_of_indices_of_interest = [0,2,3]
            temp_index = array_of_indices_of_interest[np.random.randint(0,3)] #randomly pick bins to adjust the discrepancy
            soln_round[temp_index]=soln_round[temp_index]-discrepancy
            soln_round[soln_round<0]=0
            Tat_after = np.sum(soln_round[np.array([0,2,3])])
            if temp_counter > 9999999999999:
                print('ERROR: Error in Tat mass balance.')
                break
        pTEFb_after = soln_round[1]+soln_round[2]+soln_round[3]
        pTEFb_before = Tat_pTEFb_deacetyl + Tat_pTEFb_acetyl + pTEFb_nuc # Keep in mind pTEFb_nuc has already been incremented; pTEFb_nuc does not represent amount from previous timestep
        if pTEFb_after != pTEFb_before: #mass balance (pTEFb)...care less about this than Tat because shoudl be in abundance from Host
            discrepancy = pTEFb_after - pTEFb_before
            soln_round[1] = soln_round[1]-discrepancy
            if soln_round[1]<0:
                soln_round[1] = 0
                print('ERROR: Error in pTEFb mass balance. Amt of pTEFb in nucleus went below zero')
        free_Tat=soln_round[0]
        pTEFb_nuc=soln_round[1]
        Tat_pTEFb_deacetyl=soln_round[2]
        Tat_pTEFb_acetyl=soln_round[3]
        return [free_Tat, pTEFb_nuc, Tat_pTEFb_deacetyl, Tat_pTEFb_acetyl]
            
    def evolve_state(self, timestep):
        # print('The timestep is:')
        # print(timestep)        
        
        #get variables
        host_factor_state = self.state.get_state('host_factors')
        pTEFb_nuc = host_factor_state.pTEFb_nuc        
        Tat_pTEFb_deacetyl = host_factor_state.Tat_pTEFb_deacetyl
        Tat_pTEFb_acetyl = host_factor_state.Tat_pTEFb_acetyl
        
        protein_state = self.state.get_state('proteins')
        proteins_nuc = protein_state.proteins_nuc
        
        reaction_rate_state = self.state.get_state('reaction_rates')
        #Tat_derived_transcription_rate = reaction_rate_state.Tat_derived_transcription_rate
        
        #evolve state
        
        #replenish pTFEb -- exponential doubling as Tcell grows
        pTEFb_nuc = pTEFb_nuc + np.around(self.pTEFb_NUC_INIT*np.exp(self.pTEFb_DOUBLING_RATE*(timestep+1))) - np.around(self.pTEFb_NUC_INIT*np.exp(self.pTEFb_DOUBLING_RATE*(timestep)))    
        
        #determine the effect of Tat feedback...dependent on the abundance of Tat in the nucleus
        y0 = [proteins_nuc[2], pTEFb_nuc, Tat_pTEFb_deacetyl, Tat_pTEFb_acetyl, 0]       # initial condition vector 
        # solve the ODEs
        t_seg_Tat = np.linspace(0, 59, 60)   # time grid for Tat feedback integration
        soln = odeint(self.TatODE, y0, t_seg_Tat) #use if you have scipy otherwise runge kutta
        #soln = matplotlib.mlab.rk4(TatODE, y0, tsegTat)
        
        #Accounting and discretizing and mass balance
        [free_Tat, pTEFb_nuc, Tat_pTEFb_deacetyl, Tat_pTEFb_acetyl] =  self.ODE_discretizer(soln, proteins_nuc[2], pTEFb_nuc, Tat_pTEFb_deacetyl, Tat_pTEFb_acetyl)
        proteins_nuc[2] = free_Tat

        Tat_derived_transcription_rate = np.max([0, soln[-1,4]]) #allow no negatives
        
        ##NOTE here, the ODE moves items into bin 4 = "mRNA" to indicate the number of mRNA made in the given minute
        #However, in the Transcription Class, mRNA cannot be made beyond a threshold MAX_TAT_ENHANCEMENT*BASAL_TRANSCRIPTION_RATE
        #When the ODE moves mass out of the Tat_pTEFb_acetyl bin into the mRNA bin, it recycles back Tat and pTEFb to their free forms
        #They must now re-work their ways up to the acetyl state
        #This is a source of error, however, once Tat_feedback gets past this threshold, it usually stays ON, and therefore, this error is likely not associated
        #with a big impact on system dynamics.      
                            
        #write back parameters to state object
        protein_state.protein_Nuc = proteins_nuc # Update the free Tat value
        host_factor_state.pTEFb_nuc = pTEFb_nuc # Update pTEFb value in the nucleus
        host_factor_state.Tat_pTEFb_deacetyl = Tat_pTEFb_deacetyl # Update the pTEFb deacetyl value
        host_factor_state.Tat_pTEFb_acetyl = Tat_pTEFb_acetyl # Update the pTEFb acetyl value
        reaction_rate_state.Tat_derived_transcription_rate = Tat_derived_transcription_rate       
        
        #update state to new values
        #self.state.set_state('proteins', protein_state)
        #self.state.set_state('host_factors', host_factor_state)
