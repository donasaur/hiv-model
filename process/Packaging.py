# -*- coding: utf-8 -*-
"""
Biology Summary:
1. Both Vpr and Gag levels in a cell have been implicated in the switch that
transitions the infected cell from a state of "translating" to a state of
"packaging." Here, we include effects of both signals.
A. First, a Vpr threshold is reached that drives the cell into a G2/M arrest.
There are two forms of translation in the the cell, one that is 5' cap 
dependent and one that is IRES dependent. Once in G2/M arrest, normal translation
is shut down, and only IRES dependent translation occurs. At this point, Gag
monomers may start binding to HIV RNA --Darlix et al. 
B. Once the Gag concentration in the cell reached 500nM, it may start dimerizing,
and localizing at the Plasma Membrane. --Fogarty et al. This finding is due to 
detection limits, and does not eliminate the posibility of dimers forming at a
lower concentration. What is explained by this finding is that the probability 
of dimerization increases at higher Gag concentrations. I have chosen to 
include a model of dimerization based on the velocity of Gag particles
diffusing the cytoplasm and the concentration of Gag particle in the 
cytoplasm. 

2. Once in G2/M arrest, cytoplasmic monomeric Gag (not membrane bound) can start to 
bind to full length mRNAs in the cytoplasm. There are four stem-loops (SL) in the
viral RNA, that Gag can bind to at different binding constants. Once Gag has 
bound to SL1, SL2, SL3, the RNA has a conformational change, enabling it to
bind to another RNA that also has that conformational change. --Darlix et al. 

3. The whole "nucleated" complex of two RNAs and at least 6 Gag molecules
translocates to the Plasma Membrane. The translocation time is not known,
and is therefore inferred by the diffusion time-scale. 
Assuming a sphereical stucture (this is not true--approximation!!), and 
that the "nucleated" complex has ~7x the mass of the Gag monomer, the diffusion rate 
of the dimer would be slower (factor inversely proportional to the cube root
of the mass ratio). 

[the next two steps occur throughout the entire cell cycle, not only during
G2/M arrest]

4. Monomeric cytoplasmic Gag can homodimerize. While trimers can be found,
primarity only monomers and dimers are found, and therefore I only model 
these two forms.--Balasubramaniam and Freed 2011
Gag, upon translation, reaches the Plasma Membrane in ~10mins --Balasubramaniam and Freed 2011, Ono 2009
GAG_VELOCITY is approximated as CELL_RADIUS/10min
Gag_concentration is [# Gag cytoplasmic monomers]/VOLUME_CYTOPLASM
rate_of_Gag_Gag_collision = (1/2)(pi)(GAG_DIAMETER^2)(2^(1/2))(GAG_VELOCITY)(Gag_concentration^2) --Ronis 2013
Therefore, the probability of dimerization increases with Gag cytoplasmic concentration

5. Plasma membrane localization. Monomers can translocate to the plasma membrane 
in ~10mins --Balasubramaniam and Freed 2011, Ono 2009
Assuming a sphereical stucture (this is not true--approximation!!), and 
that the dimer has twice the mass of the monomer, the diffusion rate 
of the dimer would be slower (factor inversely proportional to the cube root
of the mass ratio). 

6. It is known that a full (~1500 Gag) --Briggs et al.
particle forms in 4.5+9mins = 13.5mins. --Jouvenet et al.   
There is debate on the 1500 Gag value. One source thinks the particle size is
5000 (3000-11000) Gag --Briggs et al. And others that it is 2000 --Balasubramaniam and Freed 2011

We are modeling exponentially growing puntca, where the rate of growth in a given
timestep is dependent on the cytoplasmic Gag concentration, and the size of the puncta 
(cooperative binding). It has been found that higher Gag concentrations yeild
faster growing punta, and therefore in the model, the exponential rate constant is in part
dependent on the Gag concentration at the start of pucta growth. --Jouvenet et al.  

Gag monomers and dimers can move from the cytoplasm to the growing punta
or laterally from the membrane to the growing puncta. The majority of the Gag
added to a puncta comes from the cytoplasmic population. 

The prob of lateral movement to a puncta should be:
750/([# of Gag-dimers on Mem][13.5])   TEMP VALUE = 0.00009022
Will run model many times to get the average value of [# of Gag-dimers on Mem]
at the time of first puncta to fit this parameter. 
Once set, the higher the conc of Gag at the membrane, the faster puncta will form. 
Computationally, the size of the puncta will be determined (Poisson Dist)
when it is formed, and it will "bud off" when this size has been reached. 

7. Other proteins
A. Vif:
i. It is known that 20-100 Vif molecules exist in each viron (Lake et al. 2003)
And that Vif associates with Gag. Therefore, we assume that Vif associates with
assembling Gag at roughly the ratio of Vif:Gag in the final viron. 
ii. Vif also associates iwth Gag and the RNA to facilitate this interaction (Lake et al. 2003).
Quantitative data on this interaction is not available, so we simply check that
a non-zero amount of Vif is in the cell to allow Gag/RNA binding. 
"""

import numpy as np
from scipy.integrate import odeint
from mainaux.Process import Process
from state.Proteins import Proteins
from state.ViralProgeny import *
from mainaux.InitParamValues import *
from mainaux.TestHelpers import count_total_Gag
import math as math

#References:
#1. Krombach F, Münzing S, Allmeling AM, Gerlach JT, Behr J, Dörger M. Cell size of alveolar macrophages: an interspecies comparison. Environ Health Perspect. 1997 
#2. Darlix, J.L., Lastra, M.L., Mely, Y., Roques, B. Nucleocapsid protein chaperoning of Nucleic Acods at the heart of HIV structure assembly and cDNA synthesis. HIV sequence database. In HIV Sequence Compendium 2002. Publication no. LA-UR 02–2877. Theoretical Biology and Biophysics Group, Los Alamos National Laboratory, Los Alamos
#3. Ivanchenko, S., Godinez, W.J., Lampe, M., Kräusslich, H.G., Eils, R., Rohr, K., Bräuchle, C., Müller, B., Lamb, D.C. Dynamics of HIV-1 assembly and release. PLoS Pathog. 2009 Nov;5(11):e1000652. 
#4. Ono, A. (2009) HIV-1 Assembly at the Plasma Membrane: Gag Trafficking and Localization. Future Virol. 4: 241–257.
#5. Datta, S.A.K., Curtis, J.E., Ratcliff, W., Clark, P.K., Crist, R.M., Lebowitz, J., Krueger, S., Rein, A. (2007) Conformation of the HIV-1 Gag Protein in Solution. J Mol Biol 365: 812–824.
#6. Ronis, D. (2013) Collisions, Chemical Reactions, & Transport. Chimistry 223. McGill University. Lecture. <ronispc.chem.mcgill.ca/ronis/chem223/hnd4.pdf>
#7. Balasubramanium, M., Freed, E.O. (2011) New insights into HIV assembly and trafficing. Physiology 26: 236-251.
#8. Briggs, A., Simon, M.N., Gross, I., Krausslich, H.G., Fuller, S.D., Vogt, V.M. (2004) The stoichiometry of gag protein in HIV-1. Nat. Struct. Mol. Biol. 11: 672–675.
#9. Jouvenet, N., Simon, S.M., Bieniasz, P.D. (2011) Visualizing HIV-1 assembly. J Mol Biol. 410:501-511.
#10. Lake, J.A., Carr, J., Feng, F., Mundy, L., Burrell, C., Li, P. (2003) The role of Vif during HIV-1 infection: interaction with novel host cellular factors. J Clin Virol. 26:143-52.




#This is a Process Class
class Packaging(Process):
    def __init__(self, state, param_dict=None):
        self.state = state
        if param_dict==None:
            param_dict = generate_param_dict();         
        #Constant parameters
        self.VPR_G2ARREST_THRESH = param_dict['VPR_G2ARREST_THRESH'] #Fittable parameter, have not yet found value in the literature
        self.VOLUME_CYTOPLASM = param_dict['VOLUME_CYTOPLASM'] #L #Cell volume from Krombach et al. - Nucleus Volume (See RevBinding)
        self.AVOGADRO_NUM = param_dict['AVOGADRO_NUM']
        self.BINDING_CONSTANT_SL1 = param_dict['BINDING_CONSTANT_SL1'] #M^(-1) #Darlix et al
        self.BINDING_CONSTANT_SL2 = param_dict['BINDING_CONSTANT_SL2'] #M^(-1) #Darlix et al
        self.BINDING_CONSTANT_SL3 = param_dict['BINDING_CONSTANT_SL3'] #M^(-1) #Darlix et al
        self.BINDING_CONSTANT_SL4 = param_dict['BINDING_CONSTANT_SL4'] #M^(-1) #Darlix et al
        self.GAGNC_DISS_RATE = param_dict['GAGNC_DISS_RATE'] #Fittable made-up parameter
        self.GAG_DIFFUSION_PROB = param_dict['GAG_DIFFUSION_PROB'] #/min #Ono 2009; Balasubramaniam and Freed 2011
        self.GAG_DIMER_DIFFUSION_PROB = param_dict['GAG_DIMER_DIFFUSION_PROB'] #Assuming a sphereical stucture (this is not true--approximation!!), and that the dimer has twice the mass of the monomer, the diffusion rate of the dimer would be slower (factor inversely proportional to the cube root of the mass ratio). 
        self.GAG_VELOCITY = param_dict['GAG_VELOCITY']
        #self.GAG_VELOCITY = (((float(3)/(4*math.pi))*(self.VOLUME_CYTOPLASM*0.001))**(float(1)/3)) / (float(1)/self.GAG_DIFFUSION_PROB) #m/min
        self.GAG_DIAMETER = param_dict['GAG_DIAMETER'] #m #Datta et al. Gag is NOT globular/spherical. It is 34A x 41A, so just using mean here as an approximation. 
        self.PROB_GAG_BOUND_RNA_DIMERS = param_dict['PROB_GAG_BOUND_RNA_DIMERS'] #fittable parameter, currently  set to be not limiting at all
        self.viral_progeny_container = self.state.get_state('viral_progeny_container')
        #These are temporary parameters that will be used to obtain the parameters required to simulate viron growth
        #These parameters will be swapped out for fit "physiological" parameters after model is appropriately fit             
        
        ##PARAMETERS USED IN OLDER VERSIONS OF CODE##   
        #self.PROB_GAG_DIMERIZING_LOW_CONC =  0.008 #/min #fittable parameter
        #self.PROB_GAG_DIMERIZING_HIGH_CONC =  0.1 #/min #fittable parameter
        #self.GLOBULE_EXP_GROWTH_CONSTANT = 0.294 #1/min #k = 4.9 * 10^-3/sec #Ivanchenko et al. (average of TIRF and SDCM measurements) 
        #self.GAG_CONC_FOR_DIMERS = 500 #nM #Fogarty et al. 2013
        #self.GAG_LEVEL_FOR_DIMERS = self.GAG_CONC_FOR_DIMERS * (1/(10**9)) * (self.AVOGADRO_NUM) * (self.VOLUME_CYTOPLASM)   #(nmol/L)(mol/nmol)(molecules/mol)(L/cell)
        #self.GAG_DIFFUSION_PROB = 0.0427 #/min 
            #Calculated as: t=(x^2)/2D
            #D = 0.04um^2/sec #Liu, Minoz, Alicea (2012) Proc of Computer Science
            #V = 4990um^3 #Krombach et al
            #x = cell radius = 10.6um
            #t = 1404.7 sec = 23.4min
        ##   
        #self.PROB_GAG_GAG_INTERACT = 3.3 * (10**-8)
            #Calculated as: SA of cell / circular area of Gag molecule
            #SA = 4pir^2 = 1412.162um^2 #Krombach et al
            #Diameter of Gag = 7.8nm #1447a.a. http://www.calctool.org/CALC/prof/bio/protein_size
            #Circular area = 4.8 * 10^-5 um^2
            #Prob = 4.8 * 10^-5 um^2/1412.162um^2 = 3.3*10^-8
    
    def evolve_state(self, timestep, active_proc=[]):       
    
        #get variables
        protein_state = self.state.get_state('proteins')
        proteins_cyt = protein_state.proteins_cyt
        proteins_mem = protein_state.proteins_mem
        proteins_virion = protein_state.proteins_virion
        #gag_globule_size = protein_state.gag_globule_size
        cell_cycle_state = self.state.get_state('cell_cycle')
        cell_cycle_arrest = cell_cycle_state.cell_cycle_arrest
        #viral_globule_started = cell_cycle_state.viral_globule_started
        reaction_rates_state = self.state.get_state('reaction_rates')
        translation_suppressed = reaction_rates_state.translation_suppressed        
        mRNA_state = self.state.get_state('mRNAs')
        full_len_transcripts_cyt = mRNA_state.full_len_transcripts_cyt
        full_len_transcripts_Gag_bound = mRNA_state.full_len_transcripts_Gag_bound
        full_len_transcripts_Gag_bound = full_len_transcripts_Gag_bound.reshape((2,2,2,2))
        full_length_transcript_dimers_cyt = mRNA_state.full_length_transcript_dimers_cyt
        viral_particles_state = self.state.get_state('viral_progeny')
        viral_progeny = viral_particles_state.viral_progeny

        if active_proc == []:
            active_proc = ["cell_cycle_arrest", "dimerize_Gag_monomers", "Gag_and_Gag_dimer_diffusion", "dimerize_gagbound_transcripts", "virion_growth"]


        #Helper functions
        def site_locator():
            one_index_dict = {}
            zero_index_dict = {}
            zero_index_dict['SL1_site'] = []
            zero_index_dict['SL2_site'] = []
            zero_index_dict['SL3_site'] = []
            zero_index_dict['SL4_site'] = []
            one_index_dict['SL1_site'] = []
            one_index_dict['SL2_site'] = []
            one_index_dict['SL3_site'] = []
            one_index_dict['SL4_site'] = []            
            for i in range(2):
                for j in range(2):
                    for k in range(2):
                        for l in range(2):
                            if i == 1:
                                one_index_dict['SL1_site'] += [(i,j,k,l)] # Stores all tuples where i == 1 (should be 8)
                                zero_index_dict['SL1_site'] += [(0,j,k,l)] # Stores all tuples where i == 0 (should be 8)
                            if j == 1:
                                one_index_dict['SL2_site'] += [(i,j,k,l)]
                                zero_index_dict['SL2_site'] += [(i,0,k,l)]
                            if k == 1:
                                one_index_dict['SL3_site'] += [(i,j,k,l)]
                                zero_index_dict['SL3_site'] += [(i,j,0,l)]
                            if l == 1:
                                one_index_dict['SL4_site'] += [(i,j,k,l)]
                                zero_index_dict['SL4_site'] += [(i,j,k,0)]
            return one_index_dict, zero_index_dict
        one_index_dict, zero_index_dict = site_locator()

        def gen_binding_rates():
            rtn_dict = {}
            rtn_dict['SL1_site'] = (long(self.BINDING_CONSTANT_SL1) * long(proteins_cyt[Proteins.index['Gag']]))/(long(self.AVOGADRO_NUM) * float(self.VOLUME_CYTOPLASM))
            rtn_dict['SL2_site'] = (long(self.BINDING_CONSTANT_SL2) * long(proteins_cyt[Proteins.index['Gag']]))/(long(self.AVOGADRO_NUM) * float(self.VOLUME_CYTOPLASM))
            rtn_dict['SL3_site'] = (long(self.BINDING_CONSTANT_SL3) * long(proteins_cyt[Proteins.index['Gag']]))/(long(self.AVOGADRO_NUM) * float(self.VOLUME_CYTOPLASM))
            rtn_dict['SL4_site'] = (long(self.BINDING_CONSTANT_SL4) * long(proteins_cyt[Proteins.index['Gag']]))/(long(self.AVOGADRO_NUM) * float(self.VOLUME_CYTOPLASM))
            return rtn_dict
        binding_rates_dict = gen_binding_rates()

        def process_SL_site():
            list_of_sites = ['SL1_site', 'SL2_site', 'SL3_site', 'SL4_site']
            rate_SL1_binding = binding_rates_dict['SL1_site']
            rate_SL2_binding = binding_rates_dict['SL2_site']
            rate_SL3_binding = binding_rates_dict['SL3_site']
            rate_SL4_binding = binding_rates_dict['SL4_site']
            
            for i in range(2):
                for j in range(2):
                    for k in range(2):
                        for l in range(2):
                            if full_len_transcripts_Gag_bound[i,j,k,l] > 0:
                                sampling_values = np.random.rand(full_len_transcripts_Gag_bound[i,j,k,l])
                                if i == 0:
                                    tempRand = np.sum(sampling_values < rate_SL1_binding)
                                    delta_Gag = np.min([tempRand, proteins_cyt[Proteins.index['Gag']]])
                                    full_len_transcripts_Gag_bound[1,j,k,l] += delta_Gag
                                    full_len_transcripts_Gag_bound[i,j,k,l] -= delta_Gag
                                    proteins_cyt[Proteins.index['Gag']] -= delta_Gag
                                elif i == 1:
                                    tempRand = np.sum(sampling_values < self.GAGNC_DISS_RATE)
                                    full_len_transcripts_Gag_bound[i,j,k,l] -= tempRand
                                    full_len_transcripts_Gag_bound[0,j,k,l] += tempRand
                                    proteins_cyt[Proteins.index['Gag']] += tempRand

                            if full_len_transcripts_Gag_bound[i,j,k,l] > 0:
                                if j == 0:
                                    tempRand = np.sum(sampling_values < rate_SL2_binding)
                                    delta_Gag = np.min([tempRand, proteins_cyt[Proteins.index['Gag']]])
                                    full_len_transcripts_Gag_bound[i,1,k,l] += delta_Gag
                                    full_len_transcripts_Gag_bound[i,j,k,l] -= delta_Gag
                                    proteins_cyt[Proteins.index['Gag']] -= delta_Gag
                                elif j == 1:
                                    tempRand = np.sum(sampling_values < self.GAGNC_DISS_RATE)
                                    full_len_transcripts_Gag_bound[i,j,k,l] -= tempRand
                                    full_len_transcripts_Gag_bound[i,0,k,l] += tempRand
                                    proteins_cyt[Proteins.index['Gag']] += tempRand

                            if full_len_transcripts_Gag_bound[i,j,k,l] > 0:
                                if k == 0:
                                    tempRand = np.sum(sampling_values < rate_SL3_binding)
                                    delta_Gag = np.min([tempRand, proteins_cyt[Proteins.index['Gag']]])
                                    full_len_transcripts_Gag_bound[i,j,1,l] += delta_Gag
                                    full_len_transcripts_Gag_bound[i,j,k,l] -= delta_Gag
                                    proteins_cyt[Proteins.index['Gag']] -= delta_Gag
                                elif k == 1:
                                    tempRand = np.sum(sampling_values < self.GAGNC_DISS_RATE)
                                    full_len_transcripts_Gag_bound[i,j,k,l] -= tempRand
                                    full_len_transcripts_Gag_bound[i,j,0,l] += tempRand
                                    proteins_cyt[Proteins.index['Gag']] += tempRand

                            if full_len_transcripts_Gag_bound[i,j,k,l] > 0:                                    
                                if l == 0:
                                    tempRand = np.sum(sampling_values < rate_SL4_binding)
                                    delta_Gag = np.min([tempRand, proteins_cyt[Proteins.index['Gag']]])
                                    full_len_transcripts_Gag_bound[i,j,k,1] += delta_Gag
                                    full_len_transcripts_Gag_bound[i,j,k,l] -= delta_Gag
                                    proteins_cyt[Proteins.index['Gag']] -= delta_Gag
                                elif l == 1:
                                    tempRand = np.sum(sampling_values < self.GAGNC_DISS_RATE)
                                    full_len_transcripts_Gag_bound[i,j,k,l] -= tempRand
                                    full_len_transcripts_Gag_bound[i,j,k,0] += tempRand
                                    proteins_cyt[Proteins.index['Gag']] += tempRand                                                                                                              

        def process_cyt_proteins():
            bins_of_interest = [(1,0,0,0), (0,1,0,0), (0,0,1,0), (0,0,0,1)]
            rate_SL1_binding = binding_rates_dict['SL1_site']
            rate_SL2_binding = binding_rates_dict['SL2_site']
            rate_SL3_binding = binding_rates_dict['SL3_site']
            rate_SL4_binding = binding_rates_dict['SL4_site']            
            rate_dict = {(1,0,0,0): rate_SL1_binding, (0,1,0,0): rate_SL2_binding, (0,0,1,0): rate_SL3_binding, (0,0,0,1): rate_SL4_binding}
            for i in range(9):
                for bin in bins_of_interest:
                    tempRand = np.sum(np.random.rand(full_len_transcripts_cyt[i]) < rate_dict[bin])
                    delta_num_of_transcripts = np.min([tempRand, proteins_cyt[Proteins.index['Gag']]])
                    full_len_transcripts_Gag_bound[bin] += delta_num_of_transcripts
                    full_len_transcripts_cyt[i] -= delta_num_of_transcripts
                    proteins_cyt[Proteins.index['Gag']] -= delta_num_of_transcripts
                    proteins_cyt[Proteins.index['Rev']] += (tempRand*i)

        def process_invariants():
            full_len_transcripts_cyt[0] += full_len_transcripts_Gag_bound[0,0,0,0]
            full_len_transcripts_Gag_bound[0,0,0,0] = 0

        def dimerize_Gag_monomers():
            Gag_concentration = (proteins_cyt[Proteins.index['Gag']])/(self.VOLUME_CYTOPLASM*.001) #molecules/m^3
            rate_of_Gag_Gag_collision = (float(1)/2)*(math.pi)*(self.GAG_DIAMETER**2)*(2**(float(1)/2))*(self.GAG_VELOCITY)*(Gag_concentration**2)
            #print(rate_of_Gag_Gag_collision)
            tempRand = np.random.poisson(proteins_cyt[Proteins.index['Gag']]*rate_of_Gag_Gag_collision) #number of collisions in timestep
            if (tempRand*2)<=proteins_cyt[Proteins.index['Gag']]:
                dimers_made = tempRand
            else:
                dimers_made = np.floor((proteins_cyt[Proteins.index['Gag']])/float(2))
            proteins_cyt[Proteins.index['Gag']] -= (dimers_made*2) #2 monomers make 1 dimer
            proteins_cyt[Proteins.index['Gag_dimers']] += dimers_made

        def Gag_and_Gag_dimer_diffusion():
            tempRand = np.sum(np.random.rand(proteins_cyt[Proteins.index['Gag']])<self.GAG_DIFFUSION_PROB)
            proteins_cyt[Proteins.index['Gag']] -= tempRand
            proteins_mem[Proteins.index['Gag']] += tempRand 
            tempRand = np.sum(np.random.rand(proteins_cyt[Proteins.index['Gag_dimers']])<self.GAG_DIMER_DIFFUSION_PROB)
            proteins_cyt[Proteins.index['Gag_dimers']] -= tempRand
            proteins_mem[Proteins.index['Gag_dimers']] += tempRand

        def dimerize_gagbound_transcripts():
            tempRand = np.sum(np.random.rand(np.floor(full_len_transcripts_Gag_bound[1,1,1,:].sum()/float(2))) < self.PROB_GAG_BOUND_RNA_DIMERS)

            while tempRand > 0:
                Gag_index = 0

                for i in range(2):
                    prob_SL4_Gag = float(full_len_transcripts_Gag_bound[1,1,1,1])/full_len_transcripts_Gag_bound[1,1,1,:].sum()
                    if prob_SL4_Gag == 1:
                        full_len_transcripts_Gag_bound[1,1,1,1] -= 1
                        Gag_index += 1
                    elif prob_SL4_Gag == 0:
                        full_len_transcripts_Gag_bound[1,1,1,0] -= 1
                    else: # guaranteed at least 1 empty SL4 Gag, 1 nonempty SL4 Gag
                        if np.random.rand()<prob_SL4_Gag:
                            full_len_transcripts_Gag_bound[1,1,1,1] -= 1
                            Gag_index += 1
                        else:
                            full_len_transcripts_Gag_bound[1,1,1,0] -= 1

                self.viral_progeny_container.create_progeny(6 + Gag_index, timestep)
                tempRand -= 1                        


        #in each viral progeny:
        #item 0: state
        #   0 = simple nucleate on the membrane
        #   1 = growing puncta adding Gag and other proteins
        #   2 = complete viral particle ready for budding
        #   3 = budded viral particle
        #item 1: # of Gag
        #item 2: growth_constant
        viral_progeny_keys = viral_particles_state.viral_progeny_keys
        
        #evolve state
        #
        #1. Vpr accumulation to a certain level causes the cell to go into G2/M cell cycle arrest. Darlix et al 
        # Only place where cell_cycle_arrest gets set to 1
        if proteins_cyt[Proteins.index['Vpr']] > self.VPR_G2ARREST_THRESH:
            cell_cycle_arrest = 1
        #2. From the start of the simulation, translation of viral mRNAs progresses as normal using hold translation machinery. 
        #   In G2/M cell cycle arrest, 5'cap dependent translation is suppressed but the HIV-1 IRES sequence can ensure continued systhesis of Gap and Gag/Pol
        if cell_cycle_arrest == 1:
            translation_suppressed = 1
        #3. The nucleocapsid portion of a Gag protein can bind to a full_len_transcript_cyt (gRNA). 
        #There are 4 sites on the gRNA that can bind to a Gag molecule:
        #SL1, SL2, SL3, SL4 each with their own binding constant
        #It is unclear if there is a specified order of binding.            
        #Calculate the bindng rates
        if cell_cycle_arrest == 1 and "cell_cycle_arrest" in active_proc and proteins_cyt[Proteins.index['Vif']] > 0:
            # TODO: removed dissociation rate check
            # rate_SL1_binding = (long(self.BINDING_CONSTANT_SL1) * long(proteins_cyt[Proteins.index['Gag']]))/(long(self.AVOGADRO_NUM) * self.VOLUME_CYTOPLASM)
            # rate_SL2_binding = (long(self.BINDING_CONSTANT_SL2) * long(proteins_cyt[Proteins.index['Gag']]))/(long(self.AVOGADRO_NUM) * self.VOLUME_CYTOPLASM)
            # rate_SL3_binding = (long(self.BINDING_CONSTANT_SL3) * long(proteins_cyt[Proteins.index['Gag']]))/(long(self.AVOGADRO_NUM) * self.VOLUME_CYTOPLASM)
            # rate_SL4_binding = (long(self.BINDING_CONSTANT_SL4) * long(proteins_cyt[Proteins.index['Gag']]))/(long(self.AVOGADRO_NUM) * self.VOLUME_CYTOPLASM)
            # print(rate_SL1_binding, rate_SL2_binding, rate_SL3_binding, rate_SL4_binding)
            #First calculate SL binding in the 4D matrix
            process_SL_site()


            #Second deal with full completely unbound ones
            #Assume any bound Rev is released. This is done for accounting reasons, no evidence for this in the literature. 
            process_cyt_proteins()

            #move everything in the [0,0,0,0] bin back to full_transcript_cyt so that it gets accounted for in 1 place only. 
            process_invariants()
            #Binding of Gag-NC to the SL sites also renders the transcripts unable to translate, 
            #therefore the full_len_transcripts_Gag_bound do not have to be passed into translation
        
        #4. Allow Gag monomers in the cytoplasm to dimerize
        #Not generating a rand # per Gag in the cell to save computational time
        if "dimerize_Gag_monomers" in active_proc:
            dimerize_Gag_monomers()
       
            
        #5. Determine how many Gag proteins diffuse to the cell membrane. 
        #Both monomers and dimers can move to this compartment
        if "Gag_and_Gag_dimer_diffusion" in active_proc:
            Gag_and_Gag_dimer_diffusion()
           

        #TODO: ADD SOME RATE OF BREAKING OFF THE PM AND GETTING BACK INTO THE CYT
        #TODO: ADD PROTEIN DEGRADATION FOR THE MEM COMPARTMENT
            
        #6. Model gRNA dimerization
        #full length RNA with SL1, SL2, and SL3 bound
        num_dimerable_gRNA = np.sum(full_len_transcripts_Gag_bound[1,1,1,:])

        if num_dimerable_gRNA > 2: #if there are at least two full length RNA with SL1, SL2, and SL3 bound,
            if "dimerize_gagbound_transcripts" in active_proc:
                dimerize_gagbound_transcripts()
               

        if "virion_growth" in active_proc:
            self.viral_progeny_container.virion_growth(timestep)
          
        #9. Package other proteins with Gag
                    #TODO
                    #1. model Gag/Pro/Pol Cleavage --as separate module
                    #Processing sites in the human immunodeficiency virus type 1 (HIV-1) Gag-Pro-Pol precursor are cleaved by the viral protease at different rates
                    #2. include proteins in virons
                        #NOTES:
                        #From Sundquist, W.I., Krausslich, H.G. (2012) HIV-1 Assembly, Budding, and Maturation. Cold Springs Harbor Perspectives in Medicine 2:a006924
                        #Proteins in viron: Gag (2500), GPP, Env (7-14 trimers), tRNA primer, Vpu, Vpr (1Vpr:7Gag), Nef
                        #Env and Vpu are synthesized on the Rough ER
                        #2500 Gag/viron according to Carlson 2008
                        #GPP moves to the membrane, promoted by Gag assembly, needs to dimer for PR to act
                        #weak PR activity until budding
                        #PR cleaves 5 Gag sites and 5 GPP sites
                    #3. Env transport -- as separate module
        #NOTES ON HOW MANY VIRONS SHOULD BE MADE:
                    # Numer of virons per T cell = 10^3 to 5*10^4   (Chen et al.)
                    #These can be quickly cleared by the body
                    #There are bout 10^8 virons in the human body during infection
                    #The infected T cell lives about 1 day
                    #De Boer, R.J., Ribiero, R.M., Perelson, A.S. (2010) Current estimates for HIV-1 production imply rapid viral clearance in Lymphoid Tissues. PLoS Comput. Biol. 6:e1000906
        #10. Trigger end of simulation
                    #TODO
                    #The end of the simulation will actually be triggered by the Reverse Transcription module
                    #During RT, incomplete cytosolic viral DNA transcripts are formed.
                    #These transcripts trigger an immune response and therefore apoptotic cell death
                    #Doitsh, G., Galloway, N.K.L., Gen, X., Yang, Z., Monroe, K.M., Zepeda, O., Hunt, P.W., Hatano, H., Sowinski, S., Munoz-Arias, I., Greene, W.C. (2014) Cell death by pyroptosis drives CD4 T-cell depletion in HIV-1 infection. Nature 505: 509-514. 
                    
        #11. Model Budding
                    #TODO
                    #Baumgärtel, V., Ivanchenko, S., Dupont, A., Sergeev, M., Wiseman, P.W., Kräusslich, H.G., Bräuchle, C., Müller, B., Lamb, D.C. (2011) Live-cell visualization of dynamics of HIV budding site interactions with an ESCRT component. Nature Cell Biology 13: 469–474 
                    #A host  protein called VSP4A binds in bursts to the assembly sites
                    #2-5 VSPA4 dodecamers bind at each burst
                    #The burst takes 5s to form, 9s stable, 21s decay
                    #THe bursts occur after the 8-20min of Gag assembly in Phase II that is 47 minutes long. Phase III is the budding.
                    
        
        #write back parameters to state object
        cell_cycle_state.cell_cycle_arrest = cell_cycle_arrest
        reaction_rates_state.translation_suppressed = translation_suppressed

        # for each new protein added to ViralProgeny, also need to update proteins_virion here
        # need to find an efficient way to update num of Env
        proteins_virion[Proteins.index['Gag']] = self.viral_progeny_container.count_num_of_Gag()
        proteins_virion[Proteins.index['Vif']] = self.viral_progeny_container.count_num_of_Vif()
        proteins_virion[Proteins.index['GagProPol']] = self.viral_progeny_container.count_num_of_GagProPol()
        proteins_virion[Proteins.index['Vpr']] = self.viral_progeny_container.count_num_of_Vpr()
        proteins_virion[Proteins.index['Nef']] = self.viral_progeny_container.count_num_of_Nef()
