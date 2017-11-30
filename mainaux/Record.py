# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt 

class Record(object):
    def __init__(self, list_of_dependent_keys=None, set_of_req_record_key_names=None, combined_tracking_dict=None, sampling_rate=1):
        self.last_const_param_section = True
        self.num_of_const_param_sections = 1
        self.plot_args_dict = {}
        self.plot_args_dictVPV = {}
        self.inc_of_interest = 0
        self.use_default_log_setting = True
        self.scatter_plots_list = []
        self.sim_index = 0
        self.batch_label = ""
        self.sampling_rate = sampling_rate
        self.hist_dict = {}
        self.list_of_dependent_keys = list_of_dependent_keys
        self.set_of_req_key_names = set_of_req_record_key_names
        if combined_tracking_dict == None:
            self.variable_tracking_dict = {}
        else:
            # Uses the argument tracking dictionary as the dictionary that it
            # decides to plot
            # Will only plot dictionaries that map a key to a matrix
            self.variable_tracking_dict = combined_tracking_dict
    def add_tracking(self, time, max_timesteps, key, value):
        if key not in self.set_of_req_key_names:
            return

        if time == 0:
            #assumes all data is saved as an array of ints
            #change code to allow richer diversity of data
            self.variable_tracking_dict[key] = np.zeros((np.size(value), len(np.arange(0, max_timesteps, self.sampling_rate))), int)
            #vars(self)['saved_' + key] = np.zeros((np.size(value), max_timesteps), int) 
        
        #vars(self)['saved_' + key][:,time] = value
        self.variable_tracking_dict[key][:,time/self.sampling_rate] = value
        #self.variable_tracking_list.insert(time, (time, key, value))

    def generate_data_for_dependent_keys(self):
        for key in self.list_of_dependent_keys:
            self.variable_tracking_dict[key.name] = key.generate_key_data(self.variable_tracking_dict)

    def get_recorded_data(self, key):
        #return[vars(self)['saved_' + key]]
        return self.variable_tracking_dict[key]
