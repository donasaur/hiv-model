import numpy as np
import cPickle
import os
from mainaux.State import State
from mainaux.Process import Process

from process.TatFeedback import TatFeedback
from process.Transcription import Transcription
from process.AlternativeSplicing import AlternativeSplicing
from process.RevBinding import RevBinding
from process.MRNAExport import MRNAExport
from process.Translation import Translation
from process.ProteinLocalization import ProteinLocalization
from process.Degradation import Degradation
from process.Packaging import Packaging
from process.EnvProcessing import EnvProcessing
from mainaux.InitParamValues import *

class Simulation(object):
    
    def __init__(self, record, number_of_timesteps=360, modified_param=None, new_val_of_param=None, save_state=False, save_state_timestep_list=None):
        
        self.number_of_timesteps = number_of_timesteps
        self.current_timestep = 0
        self.record = record
        self.param_dict = generate_param_dict()

        if modified_param != None:
            self.param_dict[modified_param] = new_val_of_param

        self.save_state = save_state
        if save_state_timestep_list == None:
            self.save_state_timestep_list = [self.number_of_timesteps - 1]
        else:
            self.save_state_timestep_list = save_state_timestep_list
        
        #initialize
#        self.x = 0
#        self.t = 0
#        self.t_hold = np.zeros((1, number_of_timesteps), int)
#        self.x_hold = np.zeros((1, number_of_timesteps), int)
        #self.processes = []
#        self.storgage_vectors = ['self.x_hold', 'self.t_hold']
#        self.saved_variables = ['self.x', 'self.t']
#        self.p1 = Plot(self.number_of_timesteps, 1)
        #self.saved_variables = 'transcriptsSynthesized' ??????????????????
        #self.plots_to_be_made = []
        
        self.state = State(self.param_dict)
        
        self.process_list = []
        self.state_list= []
        self.init_processes()
        self.init_states()

    def init_processes(self):
        #This instantiates an object for each process class with the current object of the states
        self.Tat_feedback_process = TatFeedback(self.state, self.param_dict)
        self.process_list.append(self.Tat_feedback_process)
        self.Transcription_process = Transcription(self.state, self.param_dict)
        self.process_list.append(self.Transcription_process) 
        self.Alternative_splicing_process = AlternativeSplicing(self.state, self.param_dict)
        self.process_list.append(self.Alternative_splicing_process)
        self.Rev_binding_process = RevBinding(self.state, self.param_dict)
        self.process_list.append(self.Rev_binding_process)
        self.mRNA_export_process = MRNAExport(self.state, self.param_dict)
        self.process_list.append(self.mRNA_export_process)
        self.Translation_process = Translation(self.state, self.param_dict)
        self.process_list.append(self.Translation_process)
        self.Protein_localization_process = ProteinLocalization(self.state, self.param_dict)
        self.process_list.append(self.Protein_localization_process)
        self.Degradation_process = Degradation(self.state, self.param_dict)
        self.process_list.append(self.Degradation_process)
        self.Packaging_process = Packaging(self.state, self.param_dict)
        self.process_list.append(self.Packaging_process)
        self.Env_processing_process = EnvProcessing(self.state, self.param_dict)
        self.process_list.append(self.Env_processing_process)
    
    # TODO: add modified params support to states (in a clean way)
    def init_states(self):
        self.state_list.append('proteins')
        self.state_list.append('mRNAs')
        self.state_list.append('reaction_rates')
        self.state_list.append('DNAs')
        self.state_list.append('viral_progeny')
        self.state_list.append('cell_cycle')
        self.state_list.append('viral_progeny_container')
        #self.state_list.append('host_factors')

    def save_state_to_file(self, timestep):
        new_filename = "sample" + str(timestep)
        path = os.getcwd() + "/samplestates/Batch" + self.record.batch_label + "/Sim" + str(self.record.sim_index) + "/"
        fullpath = os.path.join(path, new_filename)

        main_state_notes = fullpath + "notes.txt"

        if not os.path.exists(path):
            os.makedirs(path)

        fileHandler2 = open(fullpath, 'wb')
        cPickle.dump(self.state, fileHandler2)
        fileHandler2.close()

        fileHandler3 = open(main_state_notes, 'wb')
        fileHandler3.write("Starting time of batch: " + str(self.record.batch_label) + "\n")
        fileHandler3.write("Simulation Num: " + str(self.record.sim_index) + "\n")
        fileHandler3.write("Timestep: " + str(timestep) + "\n")
        fileHandler3.close()

#    def record_for_plots(self, step, current_timestep):
#        for i in range(np.size(self.saved_variables)):
#            print "Saved variables at index %d: %d" % (i, eval(self.saved_variables[i]))
#            self.p1.record(current_timestep, eval(self.storgage_vectors[i]), eval(self.saved_variables[i]))
#            print eval(self.storgage_vectors[i])

    def run(self):
        set_of_relevant_timesteps = np.arange(0, self.number_of_timesteps, self.record.sampling_rate)
        while self.current_timestep < self.number_of_timesteps:
            #Todo: ADD line to randomize order of process list

            for process in self.process_list:
                process.evolve_state(self.current_timestep)
            for state_name in self.state_list:
                curr_state = self.state.get_state(state_name)

                if self.current_timestep in set_of_relevant_timesteps:
                    curr_state.record_state(self.record, self.current_timestep, self.number_of_timesteps)

                if self.current_timestep + 1 == self.number_of_timesteps:
                    curr_state.record_at_end(self.record)

            self.current_timestep +=1
            print('in simulation, the time is', self.current_timestep)

            if self.save_state:
                if self.current_timestep in self.save_state_timestep_list:
                    self.save_state_to_file(self.current_timestep)

        # run at the end
        self.record.generate_data_for_dependent_keys()

#            self.x = self.x + 2
#            self.t = t
#            self.record_for_plots(1, t)
            
            
