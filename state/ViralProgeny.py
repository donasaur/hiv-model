from state.Proteins import Proteins
import numpy as np
import random
from mainaux.InitParamValues import *
from mainaux.ProcessHelpers import *
#import pdb
#from mainaux.TestHelpers import count_total_Gag

class ViralProgenyContainer(object):
    def __init__(self, state):
        param_dict = generate_param_dict();

        self.helper_state_dict = {
            "NUCLEATE_CYT": ViralProgeny.NUCLEATE_CYT,
            "NUCLEATE_MEM": ViralProgeny.NUCLEATE_MEM,
            "GROWING_VIRION": ViralProgeny.GROWING_VIRION,
            "VIRION_PREBUDDING": ViralProgeny.VIRION_PREBUDDING,
            "BUDDED_VIRION": ViralProgeny.BUDDED_VIRION           
        }

        # maps protein name to a list of progeny's protein counts (for that protein name)
        self.helper_protein_dict = {}    

        self.MAX_NUM_OF_PROGENY = param_dict['MAX_NUM_OF_PROGENY'] # max num of progeny
        self.init_special_instance_vars()                
        self.init_progeny_state_count_dict()
        self.init_progeny_protein_dict()
        self.list_of_progeny = []
        self.list_of_progeny_in_creation_order = []
        self.state = state
        self.proteins_state = state.get_state('proteins')
        self.proteins_cyt = self.proteins_state.proteins_cyt
        self.proteins_mem = self.proteins_state.proteins_mem
        self.GAG_DIMER_DIFFUSION_FOLD_CHANGE = param_dict['GAG_DIMER_DIFFUSION_FOLD_CHANGE']
        self.GAG_LATERAL_DIFFUSION_FOLD_CHANGE = param_dict['GAG_LATERAL_DIFFUSION_FOLD_CHANGE'] #fittable parameter
        self.THRESH_NUCLEATE_TO_STICK_TO_MEM = param_dict['THRESH_NUCLEATE_TO_STICK_TO_MEM'] #fittable parameter #below this number of Gags, a puncta can fall off the membrane and dissociate
        self.PROB_RNA_NUCLEATE_TRANSLOCATION = param_dict['PROB_RNA_NUCLEATE_TRANSLOCATION'] #/min #fittable parameter, currently approx from the gag diffusion value        

        self.AVE_GAG_PER_VIRON = param_dict['AVE_GAG_PER_VIRON'] #molecules #Briggs et al. 
        self.AVE_VIF_PER_VIRON = param_dict['AVE_VIF_PER_VIRON'] #Frankel and Young 1998
        self.AVE_GAGPROPOL_PER_VIRON = param_dict['AVE_GAGPROPOL_PER_VIRON'] #Welker et al. 1996
        self.AVE_VPR_PER_VIRON = param_dict['AVE_VPR_PER_VIRON'] #Muller et al. 2000
        self.AVE_NEF_PER_VIRON = param_dict['AVE_NEF_PER_VIRON'] #Welker et al. 1996

        self.AVE_ENV_T_PER_VIRON = 50 # Number should be around 3-81 (Todo: find source)
        
        #These are temporary parameters that will be used to obtain the parameters required to simulate viron growth
        #These parameters will be swapped out for fit "physiological" parameters after model is appropriately fit     
        self.VIRON_EXPONENTIAL_GROWTH_CONSTANT = param_dict['VIRON_EXPONENTIAL_GROWTH_CONSTANT'] #fittable parameter
        self.NUCLEATE_DISS_RATE = param_dict['NUCLEATE_DISS_RATE'] #1/mins #Jouvenet et al.        
        
    def init_progeny_state_count_dict(self):
        for state_name in self.helper_state_dict:
            state = self.helper_state_dict[state_name]
            self.progeny_state_count[state] = 0

    # Assumes that there exists self.MAX_NUM_OF_PROGENY progeny at most
    def init_progeny_protein_dict(self):
        for protein_name in ['Gag', 'Vif', 'GagProPol', 'Vpr', 'Nef', 'Env_t', 'successful_Env_t', 'unsuccessful_Env_t']:
            self.helper_protein_dict[protein_name] = np.zeros((self.MAX_NUM_OF_PROGENY), int)

    def init_special_instance_vars(self):
        self.num_of_Gag = 0
        self.num_of_Vif = 0
        self.num_of_GagProPol = 0
        self.num_of_Vpr = 0
        self.num_of_Nef = 0
        self.num_of_Env_t = np.array([0,0,0,0])
        self.progeny_count = 0
        self.progeny_state_count = dict()
        self.prebudding_creation_time = []
        self.prebudding_elapsed_time = []

    def create_progeny(self, num_of_Gag, timestep=0):
        self.progeny_count += 1
        if self.progeny_count > self.MAX_NUM_OF_PROGENY:
            raise Exception("The number of progeny created has exceeded MAX_NUM_OF_PROGENY. Please increase the value of this param in parameters.csv")
        new_progeny = ViralProgeny(self, num_of_Gag, timestep) 

        self.list_of_progeny += [new_progeny]
        self.list_of_progeny_in_creation_order += [new_progeny] 

    def randomize_progeny(self):
        random.shuffle(self.list_of_progeny)

    def count_progeny(self):
        return self.progeny_count

    def count_num_of_protein_with_filter(self, protein_name, percent, operator):
        count = 0
        for progeny in self.list_of_progeny:
            if ops[operator](getattr(progeny, "num_of_" + protein_name), percent/float(100) * getattr(self, "AVE_" + protein_name.upper() + "_PER_VIRON") and progeny.state == ViralProgeny.VIRION_PREBUDDING):
                count += 1
        return count
        
    def count_num_of_all_protein_with_filter(self, percent1, percent2, operator):
        count = 0
        Env_t = self.count_num_of_Env_t_array()
        for progeny in self.list_of_progeny:
            if (ops[operator](getattr(progeny, "num_of_Gag"), percent1/float(100) * getattr(self, "AVE_GAG" + "_PER_VIRON")) and
               ops[operator](getattr(progeny, "num_of_Vif"), percent1/float(100) * getattr(self, "AVE_VIF" + "_PER_VIRON")) and
               ops[operator](getattr(progeny, "num_of_GagProPol"), percent1/float(100) * getattr(self, "AVE_GAGPROPOL" + "_PER_VIRON")) and
               ops[operator](getattr(progeny, "num_of_Vpr"), percent1/float(100) * getattr(self, "AVE_VPR" + "_PER_VIRON")) and
               ops[operator](getattr(progeny, "num_of_Nef"), percent2/float(100) * getattr(self, "AVE_NEF" + "_PER_VIRON")) and
               progeny.num_of_Env_t[3]>=1 
               and progeny.state == ViralProgeny.VIRION_PREBUDDING):
                count += 1
        return count

    # Todo: check to make sure getattr is different than __getattr__ (i.e. count_num_of_all_proteins_with_filter works normally)
    def __getattr__(self, attr_name):
        if attr_name.startswith("count_num_of_") and attr_name[13:] in ['Gag', 'Vif', 'GagProPol', 'Vpr', 'Nef']:
            def count_num_of_protein():
                return getattr(self, "num_of_" + attr_name[13:])
            return count_num_of_protein
        return super(object, self).__getattribute__(attr_name)
        
    # # Make more efficient
    # def count_num_of_Gag(self):
    #     count = 0
    #     for progeny in self.list_of_progeny:
    #         count += progeny.num_of_Gag
    #     return count

    # # Make more efficient
    # def count_num_of_Vif(self):
    #     count = 0
    #     for progeny in self.list_of_progeny:
    #         count += progeny.num_of_Vif
    #     return count

    # # Make more efficient
    # def count_num_of_GagProPol(self):
    #     count = 0
    #     for progeny in self.list_of_progeny:
    #         count += progeny.num_of_GagProPol
    #     return count

    # # Make more efficient
    # def count_num_of_Vpr(self):
    #     count = 0
    #     for progeny in self.list_of_progeny:
    #         count += progeny.num_of_Vpr
    #     return count

    # # Make more efficient
    # def count_num_of_Nef(self):
    #     count = 0
    #     for progeny in self.list_of_progeny:
    #         count += progeny.num_of_Nef
    #     return count     

    def count_num_of_Env_t_array(self):
        return self.num_of_Env_t

    def count_total_num_of_Env_t(self):
        return sum(self.num_of_Env_t)

    def count_num_of_Env_t(self, num_of_successes_in_trimer):
        return self.num_of_Env_t[num_of_successes_in_trimer]

    def count_progeny_in_different_state(self):
        state_vector = np.zeros((len(self.helper_state_dict)), int)
        for state_name in self.helper_state_dict:
            state = self.helper_state_dict[state_name]
            i = state - 1
            state_vector[i] = self.progeny_state_count[state]
        return state_vector

    # Returns a vector of the states of progeny in order of creation
    # Assumes that there exists self.MAX_NUM_OF_PROGENY progeny at most
    # TODO: allocate more space as necessary
    def progeny_state_vector(self):
        rtn_vector = np.zeros((self.MAX_NUM_OF_PROGENY), int)
        for i, progeny in enumerate(self.list_of_progeny_in_creation_order):
            rtn_vector[i] = progeny.state
        return rtn_vector


    def progeny_protein_vector(self, protein_name, state=None):
        if state != None:
            rtn_vector = []
            for progeny in self.list_of_progeny_in_creation_order:
                if progeny.state == state:
                    if protein_name == 'Env_t':
                        rtn_vector.append(progeny.num_of_Env_t.sum())
                    elif protein_name == 'successful_Env_t':
                        rtn_vector.append(progeny.num_of_Env_t[3])
                    elif protein_name == 'unsuccessful_Env_t':
                        rtn_vector.append(progeny.num_of_Env_t[:3].sum())
                    else:
                        rtn_vector.append(getattr(progeny, "num_of_" + protein_name))
        else:
            rtn_vector = self.helper_protein_dict[protein_name]
        return rtn_vector

    def progeny_elapsed_time_until_prebudding_vector(self):
        rtn_vector = np.zero((self.MAX_NUM_OF_PROGENY), int)

    def consider_actual_progeny_only(self, vector):
        progeny_count = self.count_progeny()
        return vector[0:progeny_count]

    def record_state(self, record, timestep, max_timesteps):
        record.add_tracking(timestep, max_timesteps, 'progeny_count', self.count_progeny())
        record.add_tracking(timestep, max_timesteps, 'progeny_state_count', self.count_progeny_in_different_state())     
        record.add_tracking(timestep, max_timesteps, 'state_of_diff_progeny', self.progeny_state_vector())

        record.add_tracking(timestep, max_timesteps, 'total_num_of_virion_Gag', self.count_num_of_Gag())
        record.add_tracking(timestep, max_timesteps, 'num_of_Gag_of_diff_progeny', self.progeny_protein_vector('Gag'))


        record.add_tracking(timestep, max_timesteps, 'total_num_of_virion_Vif', self.count_num_of_Vif())
        record.add_tracking(timestep, max_timesteps, 'num_of_Vif_of_diff_progeny', self.progeny_protein_vector('Vif'))

        record.add_tracking(timestep, max_timesteps, 'total_num_of_virion_GagProPol', self.count_num_of_GagProPol())
        record.add_tracking(timestep, max_timesteps, 'num_of_GagProPol_of_diff_progeny', self.progeny_protein_vector('GagProPol'))

        record.add_tracking(timestep, max_timesteps, 'total_num_of_virion_Vpr', self.count_num_of_Vpr())
        record.add_tracking(timestep, max_timesteps, 'num_of_Vpr_of_diff_progeny', self.progeny_protein_vector('Vpr'))

        record.add_tracking(timestep, max_timesteps, 'total_num_of_virion_Nef', self.count_num_of_Nef())
        record.add_tracking(timestep, max_timesteps, 'num_of_Nef_of_diff_progeny', self.progeny_protein_vector('Nef'))

        record.add_tracking(timestep, max_timesteps, 'total_num_of_virion_Env_trimer', self.count_num_of_Env_t_array())
        record.add_tracking(timestep, max_timesteps, 'num_of_Env_t_of_diff_progeny', self.progeny_protein_vector('Env_t'))

        record.add_tracking(timestep, max_timesteps, 'num_virion_Gag_40_ge', self.count_num_of_protein_with_filter("Gag", 40, ">="))
        record.add_tracking(timestep, max_timesteps, 'num_viable_virions', self.count_num_of_all_protein_with_filter(10., 1.5, ">="))

            # def print_helper(protein_name):
            #     protein_of_interest_vector = getattr(self, "progeny_num_of_" + protein_name + "_vector")(ViralProgeny.VIRION_PREBUDDING)
            #     if len(protein_of_interest_vector) == 0:
            #         return                
            #     print("Avg num of " + protein_name + " per Virion: " + str(np.average(protein_of_interest_vector)))
            #     print("Std dev for " + protein_name + " per Virion: " + str(np.std(protein_of_interest_vector)))
            #     print("Min " + protein_name + " per Virion: " + str(np.min(protein_of_interest_vector)))
            #     print("Max " + protein_name + " per Virion: " + str(np.max(protein_of_interest_vector)))

#            print_helper("Gag")
#            print_helper("Vif")
#            print_helper("GagProPol")
#            print_helper("Vpr")
#            print_helper("Nef")
#            print_helper("Successful_Env_t")
#            print_helper("Unsuccessful_Env_t")

    def record_at_end(self, record):
        for protein_name in ['Gag', 'Vif', 'GagProPol', 'Vpr', 'Nef', 'Env_t', 'successful_Env_t', 'unsuccessful_Env_t']:
            record.hist_dict['num_of_' + protein_name + '_of_diff_progeny'] = self.consider_actual_progeny_only(self.progeny_protein_vector(protein_name))

    def virion_growth(self, timestep=0):
        self.randomize_progeny()
        for progeny in self.list_of_progeny: # each iteration deals with 1 progeny
            if progeny.state == ViralProgeny.NUCLEATE_CYT:
                tempRand = np.random.rand() < self.PROB_RNA_NUCLEATE_TRANSLOCATION
                if tempRand:
                    progeny.update_state(ViralProgeny.NUCLEATE_MEM)

            if progeny.state == ViralProgeny.NUCLEATE_MEM or progeny.state == ViralProgeny.GROWING_VIRION:
                if progeny.growth_const == None:
                    # print "self.AVE_GAG_PER_VIRON: " + str(self.AVE_GAG_PER_VIRON)
                    # print "self.VIRON_EXPONENTIAL_GROWTH_CONSTANT: " + str(self.VIRON_EXPONENTIAL_GROWTH_CONSTANT)
                    # print "Gag protein in protein cyt: " + str(self.proteins_cyt[Proteins.index['Gag']])
                    # print "Gag dimers: " + str(self.proteins_cyt[Proteins.index['Gag_dimers']])
                    denominator_of_growth_const = (self.VIRON_EXPONENTIAL_GROWTH_CONSTANT * (self.proteins_cyt[Proteins.index['Gag']]+(2 * self.proteins_cyt[Proteins.index['Gag_dimers']])))
                    if denominator_of_growth_const > 0:
                        progeny.growth_const = (float(self.AVE_GAG_PER_VIRON)/denominator_of_growth_const)
                    else:
                        progeny.growth_const = 0
                if progeny.final_Gag_count == 0:
                    #TODO: Add a better representation of distribution of punta sizes
                    progeny.final_Gag_count = np.random.poisson(self.AVE_GAG_PER_VIRON)

                #this rate of growth exponenetially grows as virion grows   
                rate_of_Gag_mon_cyt_binding = self.proteins_cyt[Proteins.index['Gag']] * progeny.num_of_Gag * progeny.growth_const
                # if timestep >= 4500:
                #     print "progeny.growth_const: " + str(progeny.growth_const)
                #     print "progeny.num_of_Gag: " + str(progeny.num_of_Gag)
                #     print "self.proteins_cyt[Proteins.index['Gag']]: " + str(self.proteins_cyt[Proteins.index['Gag']])
                #     print "rate_of_Gag_mon_cyt_binding: " + str(rate_of_Gag_mon_cyt_binding)
                rate_of_Gag_dim_cyt_binding = self.proteins_cyt[Proteins.index['Gag_dimers']] * progeny.num_of_Gag * progeny.growth_const * self.GAG_DIMER_DIFFUSION_FOLD_CHANGE
                rate_of_Gag_mon_mem_binding = self.proteins_mem[Proteins.index['Gag']] * progeny.num_of_Gag * progeny.growth_const * self.GAG_LATERAL_DIFFUSION_FOLD_CHANGE
                if rate_of_Gag_mon_mem_binding < 0: # todo: remove block after bug found and fixed
                    print "proteins_mem[Gag]: " + str(self.proteins_mem[Proteins.index['Gag']])
                    print "num_of_Gag: " + str(progeny.num_of_Gag)
                    print "growth_constant: " + str(progeny.growth_const)
                    print "GAG_LATERAL_DIFFUSION_FOLD_CHANGE: " + str(self.GAG_LATERAL_DIFFUSION_FOLD_CHANGE)
                # if timestep >= 4500:
                #     print "progeny.growth_const: " + str(progeny.growth_const)
                #     print "rate_of_Gag_mon_mem_binding: " + str(rate_of_Gag_mon_mem_binding)
                rate_of_Gag_dim_mem_binding = self.proteins_mem[Proteins.index['Gag_dimers']] * progeny.num_of_Gag * progeny.growth_const * self.GAG_DIMER_DIFFUSION_FOLD_CHANGE * self.GAG_LATERAL_DIFFUSION_FOLD_CHANGE

                #a. movement of Gag cytoplasmic monomers
                number_of_Gag_still_needed = progeny.final_Gag_count - progeny.num_of_Gag                
                tempRand = min(np.random.poisson(rate_of_Gag_mon_cyt_binding), self.proteins_cyt[Proteins.index['Gag']], number_of_Gag_still_needed)
                # print "tempRand: " + str(tempRand)
                # print "rate_of_Gag_mon_cyt_binding: " + str(rate_of_Gag_mon_cyt_binding)
                # print "Number of Gag in Proteins Cyt : " + str(self.proteins_cyt[Proteins.index['Gag']])
                self.proteins_cyt[Proteins.index['Gag']] = self.proteins_cyt[Proteins.index['Gag']] - tempRand
                progeny.update_num_of_Gag(tempRand)

                # b. movement of Gag cytoplasmic dimers
                number_of_Gag_dimer_still_needed = int((progeny.final_Gag_count - progeny.num_of_Gag)*0.5)
                tempRand = min(np.random.poisson(rate_of_Gag_dim_cyt_binding), self.proteins_cyt[Proteins.index['Gag_dimers']], number_of_Gag_dimer_still_needed)
                self.proteins_cyt[Proteins.index['Gag_dimers']] = self.proteins_cyt[Proteins.index['Gag_dimers']] - tempRand
                progeny.update_num_of_Gag(2*tempRand)

                # c. lateral movement of Gag monomers
                number_of_Gag_still_needed = progeny.final_Gag_count - progeny.num_of_Gag      
                tempRand = min(np.random.poisson(rate_of_Gag_mon_mem_binding), self.proteins_mem[Proteins.index['Gag']], number_of_Gag_still_needed)
                self.proteins_mem[Proteins.index['Gag']] = self.proteins_mem[Proteins.index['Gag']] - tempRand
                progeny.update_num_of_Gag(tempRand)

                # d. lateral movement of Gag dimers
                number_of_Gag_dimer_still_needed = int((progeny.final_Gag_count - progeny.num_of_Gag)*0.5)
                tempRand = min(np.random.poisson(rate_of_Gag_dim_mem_binding), self.proteins_mem[Proteins.index['Gag_dimers']], number_of_Gag_dimer_still_needed)
                self.proteins_mem[Proteins.index['Gag_dimers']] = self.proteins_mem[Proteins.index['Gag_dimers']] - tempRand
                progeny.update_num_of_Gag(2*tempRand)

                total_rate_of_Gag_binding = rate_of_Gag_mon_cyt_binding + (2*rate_of_Gag_dim_cyt_binding) + rate_of_Gag_mon_mem_binding + (2*rate_of_Gag_dim_mem_binding) # per viron per second

                def binding_helper(protein_name):
                    binding_rate_of_protein = getattr(self, "AVE_" + protein_name.upper() + "_PER_VIRON")/float(self.AVE_GAG_PER_VIRON) * total_rate_of_Gag_binding
                    proteins_cyt = getattr(self, "proteins_cyt")
                    tempRand = min(np.random.poisson(binding_rate_of_protein), proteins_cyt[Proteins.index[protein_name]])
                    proteins_cyt[Proteins.index[protein_name]] -= tempRand
                    getattr(progeny, "update_num_of_" + protein_name)(tempRand)

                binding_helper('Vif')
                binding_helper('GagProPol')
                binding_helper('Vpr')
                binding_helper('Nef')                             

                if progeny.num_of_Gag > self.THRESH_NUCLEATE_TO_STICK_TO_MEM and progeny.state != ViralProgeny.GROWING_VIRION:
                    progeny.update_state(ViralProgeny.GROWING_VIRION, timestep) #progeny is now in the "growing state" can cannot fall off

                #end puncta if it is bigger than the average size. 
                if progeny.num_of_Gag >= progeny.final_Gag_count:
                    progeny.update_state(ViralProgeny.VIRION_PREBUDDING, timestep) # already to bud!

                # If no gag binds in a certain amount of time, nuclate falls off
                #does it dissociate? Dissociation not modelled here. 
                if progeny.state == ViralProgeny.NUCLEATE_MEM:
                    if np.random.rand() < self.NUCLEATE_DISS_RATE:
                        progeny.update_state(ViralProgeny.NUCLEATE_CYT)                                         

# state = ViralProgeny.NUCLEATE_CYT (not bound to membrane), ViralProgeny.NUCLEATE_MEM (bound to membrane), ViralProgeny.GROWING_VIRION, ViralProgeny.VIRION_PREBUDDING, ViralProgeny.BUDDED_VIRION

class ViralProgeny(object):
    # static var
    name = 0
    NUCLEATE_CYT = 1
    NUCLEATE_MEM = 2
    GROWING_VIRION = 3
    VIRION_PREBUDDING = 4
    BUDDED_VIRION = 5    

    def __init__(self, container, num_of_Gag, timestep=0):
        self.container = container
        self.state = None        
        self.update_state(ViralProgeny.NUCLEATE_CYT)
        self.num_of_Gag = 0
        self.num_of_Vif = 0
        self.num_of_GagProPol = 0
        self.num_of_Vpr = 0
        self.num_of_Nef = 0
        self.update_num_of_Gag(num_of_Gag)
        self.growth_const = None
        self.name = ViralProgeny.name
        ViralProgeny.name += 1
        self.growing_timestep = timestep
        self.prebud_timestep = timestep
        self.final_Gag_count = 0
        self.num_of_Env_t = np.array([0,0,0,0])

    def update_state(self, state, timestep=0):
        prev_state = self.state
        self.state = state

        if (prev_state != None):
            self.container.progeny_state_count[self.state] += 1
            self.container.progeny_state_count[prev_state] -= 1
        else: # prev_state == None
            self.container.progeny_state_count[self.state] += 1

        if state == ViralProgeny.GROWING_VIRION:
            self.growing_timestep = timestep

        if state == ViralProgeny.VIRION_PREBUDDING:
            self.prebud_timestep = timestep
            self.container.prebudding_creation_time += [self.growing_timestep]
            self.container.prebudding_elapsed_time += [self.prebud_timestep-self.growing_timestep]

    def __getattr__(self, attr_name):
        protein_name = attr_name[14:]
        if attr_name.startswith("update_num_of_") and protein_name in ['Gag', 'Vif', 'GagProPol', 'Vpr', 'Nef']:
            def update_num_of_protein(protein_num_added):
                setattr(self, "num_of_" + protein_name, getattr(self, "num_of_" + protein_name) + protein_num_added)
                setattr(self.container, "num_of_" + protein_name, getattr(self.container, "num_of_" + protein_name) + protein_num_added)
                self.container.helper_protein_dict[protein_name][self.name] += protein_num_added
            return update_num_of_protein
        return super(object, self).__getattribute__(attr_name)

    def update_num_of_Env_t(self, num_of_successes_in_trimer, num_of_Env_t_added):
        self.num_of_Env_t[num_of_successes_in_trimer] += num_of_Env_t_added
        self.container.num_of_Env_t[num_of_successes_in_trimer] += num_of_Env_t_added
        #pdb.set_trace()
        self.container.helper_protein_dict['Env_t'][self.name] += num_of_Env_t_added
        if num_of_successes_in_trimer in range(3):
            self.container.helper_protein_dict['unsuccessful_Env_t'][self.name] += num_of_Env_t_added
        else:
            self.container.helper_protein_dict['successful_Env_t'][self.name] += num_of_Env_t_added
