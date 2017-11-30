from mainaux.State import State
from mainaux.Process import Process
import numpy as np
from state.Proteins import Proteins
from mainaux.InitParamValues import *
from mainaux.ProcessHelpers import * # this is where roll_dice, transfer_buckets is defined
from state.ViralProgeny import *

class EnvProcessing(Process):
    def __init__(self, state, param_dict=None):
        self.state = state

        if param_dict==None:
            param_dict = generate_param_dict()

        #Vars made up in this module
        self.Oligosaccharyltransferase = 1000
        self.GlucosidaseI = 1000
        self.GlucosidaseII = 1000

        #Parameters (Note: may need to fit and/or look up in paper)
        #Unknown
        self.rate_of_ER_localization = 1
        self.rate_of_Golgi_localization = 1
        self.rate_glucosidaseII = 87000
        self.rate_golgi_glycosylation = 1

        #Known
        self.rate_glucosidaseI = 87000 #1/min #Heterologous expression and characterization of processing alpha-glucosidase I from Aspergillus brasiliensis ATCC 9642 Miyazaki, T.; Matsumoto, Y.; Matsuda, K.; Kurakata, Y.; Matsuo, I.; Ito, Y.; Nishikawa, A.; Tonozuka, T.; Glycoconj. J. 28, 563-571 (2011)
        self.rate_Env_folding = 0.00333 #1/min #Vol. 17, 2003, The FASEB Journal, Folding of HIV-1 Envelope glycoprotein involves extensive isomerization of disulfide bonds and conformation-dependent leader peptide cleavage, AAFKE LAND, DUCO ZONNEVELD, AND INEKE BRAAKMAN; 1058-1067
        self.rate_oligosaccharyltransferase = 1.2 #Molecular Cell, Vol. 12, 101-111, 2003, Oligosaccharyltransferase Isoforms that Contain Different Catalytic STT3 Subunits Have Distinct Enzymatic Properties, Daniel J. Kelleher, Denise Karaoglu, Elisabet C. Mandon, and Reid Gilmore
        self.prob_golgi_glycosylation_error = 0.5 #Nature 446, 1038-1045 (2007) Exploiting the defensive sugars of HIV-1 for drug and vaccine design. Christopher N. Scanlan, John Offer, Nicole Zitzmann & Raymond A. Dwek
        self.rate_trimerization = .0083 #Land, A., and Braakman, I. (2001) FOlding of the human immunodeficiency virus type 1 envelope glycoprotein in the endoplasmic reticulum. Biochimie 83, 783-790. 17.
        self.rate_membrane_localization = 180000 #Virology, Volume 324, 2004, Pages 90-102, HIV-1 acute infection env glycomutants designed from 3D model: effects on processing, antigenicity, and neutralization sensitivity; Frederic Reynard, Ahmed Fatmi, Bernard Verrier, Frederic Bedin
        self.rate_cleavage = 0.00333 #1/min #Approx from Figure 2A of: #Inhibition of furin mediated cleavage activation of HIV_1 glycoprotein Env. Hallenbuerger, S., Bosch, V., Angliker, H., Shaw, E., Klenk, H.D., Garten, W. (1992) Nature Vol 360, 358-361.
        #There should be around 1 successful trimer on each viron and 12 (7-14) unsuccessful trimers on each viron: Klasse, 2007, Virology, 369 260-264
        #rate_viron_incorporation = 42 #1/min #~42 Env trimer binding events per min. In other words of the many virons, in a given minute 42 will get an env. 
        self.rate_viron_incorporation = 84 #fit to get average of 12 trimers per viron

    def evolve_state(self, timestep, final_step=10, single_step=None):
        # Assume it is
        protein_state = self.state.get_state('proteins')
        proteins_cyt = protein_state.proteins_cyt
        proteins_virion = protein_state.proteins_virion
        env_misc = protein_state.env_misc
        viral_progeny_container = self.state.get_state('viral_progeny_container')

        env_misc['Env : cytoplasm'] = proteins_cyt[Proteins.index['Env']]

        # Notation: -> means Env turns to Env : ER

        #Step 1. Env proteins from the cytoplasm are brought into the ER
        #assumes rate_of_ER_localization <= 1
        if (single_step == None and final_step >= 1) or (single_step == 1):
            transfer_buckets(env_misc, 'Env : cytoplasm', 'Env : ER', self.rate_of_ER_localization)

        #Step 2. Env proteins in the ER are glycosylated by oligosaccharyltransferase
        if (single_step == None and final_step >= 2) or (single_step == 2):
            rate_of_step = self.Oligosaccharyltransferase*self.rate_oligosaccharyltransferase
            transfer_buckets(env_misc, 'Env : ER', 'Env : ER : G1', rate_of_step)

        #Step 3. Env proteins in the ER are glycosylated by glucosidase I and II
        #assumes rate is very high
        if (single_step == None and final_step >= 3) or (single_step == 3):
            if self.GlucosidaseI * self.rate_glucosidaseI >= env_misc['Env : ER : G1']:
                env_misc['Env : ER : G2'] = env_misc['Env : ER : G2'] + env_misc['Env : ER : G1']
                env_misc['Env : ER : G1'] = 0
            else:
                rate_of_step = self.GlucosidaseI*self.rate_glucosidaseI
                transfer_buckets(env_misc, 'Env : ER : G1', 'Env : ER : G2', rate_of_step)

            if self.GlucosidaseII * self.rate_glucosidaseII >= env_misc['Env : ER : G2']:
                env_misc['Env : ER : G3'] = env_misc['Env : ER : G3'] + env_misc['Env : ER : G2']
                env_misc['Env : ER : G2'] = 0
            else:
                rate_of_step = self.GlucosidaseII*self.rate_glucosidaseII
                transfer_buckets(env_misc, 'Env : ER : G2', 'Env : ER : G3', rate_of_step)

        #Step 4. Folding by Calnexin chaperone
        if (single_step == None and final_step >= 4) or (single_step == 4):
            transfer_buckets(env_misc, 'Env : ER : G3', 'Env : ER : G3 : folded', self.rate_Env_folding)

            #Step 4.5. Second glucosidase II reaction
            if self.GlucosidaseII * self.rate_glucosidaseII >= env_misc['Env : ER : G3 : folded']:
                env_misc['Env : ER : G4 : folded'] = env_misc['Env : ER : G4 : folded'] + env_misc['Env : ER : G3 : folded']
                env_misc['Env : ER : G3 : folded'] = 0
            else:        
                rate_of_step = self.GlucosidaseII*self.rate_glucosidaseII
                transfer_buckets(env_misc, 'Env : ER : G3 : folded', 'Env : ER : G4 : folded', rate_of_step)

            #Need to add error rate and degradation pathway -- cannot add because no parameters found

        #Step 5. Transport to Golgi
        if (single_step == None and final_step >= 5) or (single_step == 5):
            transfer_buckets(env_misc, 'Env : ER : G4 : folded', 'Env : Golgi', self.rate_of_Golgi_localization)

        #Step 6. Golgi glycosylation.
        if (single_step == None and final_step >= 6) or (single_step == 6):
            amt_transferred = transfer_buckets(env_misc, 'Env : Golgi', 'Env : Golgi : G5', self.rate_golgi_glycosylation)
            amt_errored = roll_dice(amt_transferred, self.prob_golgi_glycosylation_error)
            move_buckets(env_misc, 'Env : Golgi : G5', 'Env : Golgi : G5 : error', amt_errored)

        #Step 7. Trimerization
        #A. Determine how many trimers will form this second
        # Note: changed around what the index of this size-4 array means (index = # of successes in trimer)
        if (single_step == None and final_step >= 7) or (single_step == 7):
            total_Golgi_G5_Env = env_misc['Env : Golgi : G5'] + env_misc['Env : Golgi : G5 : error']
            num_trimers_created = roll_dice(total_Golgi_G5_Env/float(3), self.rate_trimerization)
            for i in range(num_trimers_created): # iterate through each trimer
                successes_in_trimer = 0
                for j in range(3):
                    total_Golgi_G5_Env = env_misc['Env : Golgi : G5'] + env_misc['Env : Golgi : G5 : error']
                    successful_trimer_cutoff = float(env_misc['Env : Golgi : G5'])/total_Golgi_G5_Env
                    if np.random.rand()<successful_trimer_cutoff: # if Env protein is successful one as part of trimer
                        successes_in_trimer = successes_in_trimer + 1
                        env_misc['Env : Golgi : G5'] -= 1
                    else:
                        env_misc['Env : Golgi : G5 : error'] -= 1
                    # index of env_misc['Env : trimers'] is the # of successful Envs in trimer
                env_misc['Env : trimers'][successes_in_trimer] += 1

        #Step 8. Cleavage and non-covelent complexation (gp160 breaking up into gp120 and gp140, and gp120 and gp140 reassociating)
        #since we do not have rates of non-covenlent complexation, will assume that it is quick and happens in the same step as cleavage
        if (single_step == None and final_step >= 8) or (single_step == 8):
            if sum(env_misc['Env : trimers']) > 0:
                for i in range(4):
                    temp_rand = roll_dice(env_misc['Env : trimers'][i], self.rate_cleavage)
                    env_misc['Env : trimers : cleaved'][i] += temp_rand
                    env_misc['Env : trimers'][i] -= temp_rand

        #Step 9. Membrane localization
        #Todo: (Address this) Currently, env trimers in the membrane are not included as part of proteins_mem
        if (single_step == None and final_step >= 9) or (single_step == 9):
            if sum(env_misc['Env : trimers : cleaved']) > 0:
                for i in range(4):
                    if self.rate_membrane_localization >= env_misc['Env : trimers : cleaved'][i]:
                        env_misc['Env : trimers : membrane'][i] += env_misc['Env : trimers : cleaved'][i]
                        env_misc['Env : trimers : cleaved'][i] = 0
                    else:
                        tempRand = np.random.poisson(self.rate_membrane_localization)
                        env_misc['Env : trimers : membrane'][i] += np.min([tempRand, env_misc['Env : trimers : cleaved'][i]])
                        env_misc['Env : trimers : cleaved'][i] -= np.min([tempRand, env_misc['Env : trimers : cleaved'][i]])

        #Step 10. Addition to Virons
        if (single_step == None and final_step >= 10) or (single_step == 10):
            num_of_progeny = viral_progeny_container.count_progeny()
            if sum(env_misc['Env : trimers : membrane']) > 0 and num_of_progeny > 0:
                num_trimers_to_bind_virons = np.random.poisson(self.rate_viron_incorporation)
                for i in range(min([num_trimers_to_bind_virons,int(sum(env_misc['Env : trimers : membrane']))])):
                    env_misc['Env : trimers : membrane'] = np.array(env_misc['Env : trimers : membrane']).astype('float') # TODO: think about the logic here further
                    cum_trimer_vector = np.cumsum((env_misc['Env : trimers : membrane'])/float(sum(env_misc['Env : trimers : membrane'])))
                    temp_rand = np.random.rand() 
                    viron_to_bind = np.random.randint(0, num_of_progeny)
                    progeny_of_interest = viral_progeny_container.list_of_progeny[viron_to_bind]
                    if temp_rand < cum_trimer_vector[0]:
                        progeny_of_interest.update_num_of_Env_t(0, 1)
                        env_misc['Env : trimers : membrane'][0] -= 1
                    elif temp_rand < cum_trimer_vector[1]:
                        progeny_of_interest.update_num_of_Env_t(1, 1)
                        env_misc['Env : trimers : membrane'][1] -= 1
                    elif temp_rand < cum_trimer_vector[2]:
                        progeny_of_interest.update_num_of_Env_t(2, 1)
                        env_misc['Env : trimers : membrane'][2] -= 1
                    else:
                        progeny_of_interest.update_num_of_Env_t(3, 1)
                        env_misc['Env : trimers : membrane'][3] -= 1

        # Write back to Env bucket in proteins_cyt
        proteins_cyt[Proteins.index['Env']] = env_misc['Env : cytoplasm']
        env_misc['Env : cytoplasm'] = 0 # 'Env : cytoplasm' is just a bin to play with for the duration of the process
        proteins_virion[Proteins.index['Env']] = 3 * viral_progeny_container.count_num_of_Env_t_array().sum()
