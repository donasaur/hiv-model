import numpy as np
import operator
from mainaux.Record import Record
from mainaux.Simulation import Simulation
import mainaux.PlotCompilation as PlotCompilation
import mainaux.PlotHelpers as PlotHelpers
import matplotlib.pyplot as plt
import os
import csv
import datetime

# Dictionary of supported operators
_ops =  {
            "<" : operator.lt,
            "<=": operator.le,
            "==": operator.eq,
            "!=": operator.ne,
            ">=": operator.ge,
            ">" : operator.gt
        }

# Returns a dictionary where the key is a key in variable_tracking_dict
# and the value is a list of matrices, where each matrix corresponds to
# a particular simulation.

# A matrix for a particular simulation is to be added to the list of matrices
# for a particular key of the output dictionary if:
# operator(matrix[row_index, timestep], value) == True

# Omit the operator argument to have generate_subset_dict act as the identity function
# Also added timestep=None, so can filter out matrices depending on whether or not
# a particular row has all its values meet a certain condition
# e.g. key='proteins_nuc', row_index=Proteins.index['Rev'], timestep=None, value=0, operator="<"
# means that we are only interested in simulations where the amount of Rev in the nucleus < 0
# for all timesteps
def generate_subset_dict(input_dict, key=None, row_index=None, timestep=None, value=None, operator=None):
    subset_dict = dict()

    if operator == None:
        return input_dict

    # input_dict[key] is a list of matrices corresponding to a particular key, e.g. 'full_len_transcripts_nuc'
    # list_of_indices is a list of all the indices of the simulations that meet a certain criteria
    list_of_indices = []
    if timestep != None:
        for sim_index in range(len(input_dict[key])):
            if _ops[operator](input_dict[key][sim_index][row_index, timestep], value) == True:
                list_of_indices.append(sim_index)
    else: # timestep == None, checks all 360 timesteps of a particular row with index row_index
        for sim_index in range(len(input_dict[key])):
            if len(input_dict[key][sim_index][row_index][_ops[operator](input_dict[key][sim_index][row_index], value)]) > 0:
                list_of_indices.append(sim_index)

    # Add relevant matrices to each key in subset_dict
    for key in input_dict:
        subset_dict[key] = []
        for sim_index in list_of_indices:
            subset_dict[key].append(input_dict[key][sim_index])

    return subset_dict

# num_of_sim: Number of simulations to run per const param value (section)
# num_of_timesteps: Number of timesteps to run for all plots
# param_name: Name of the parameter whose value will be varied
# init_val: initial value of the parameter
# increment_val: value at which we increment the parameter value by
# increment_type: options are 'linear', 'exponential'
#   'linear' : init_val + increment_val*curr_increment
#   'exponential' : init_val * (increment_val)^curr_increment
#   where curr_increment ranges between 0 and num_increments (inclusive)
# num_of_increments: how many times should we increment the value of the param
# type_of_plot: determines how all the simulation data for a particular const param section will be displayed
#   options: "standard deviation", "range", "all data", None
#   Note average is always plotted for >1 sim. per const param section
# list_of_key_names: determines what keys will be plotted
def generate_increment_param_plots(num_of_sim, num_of_timesteps, param_name, init_val, increment_val, increment_type, num_increments, export_raw_data_no_plot, group_by_row_num, type_of_plot, log_setting_opt, sampling_rate, list_of_key_names):

    if export_raw_data_no_plot:
        output_data_label = initialize_label("outputdata")
        create_output_data_info_file(output_data_label, group_by_row_num, num_of_sim, param_name, increment_type)

    # each iteration of this for loop corresponds to a parameter cross section in VPV
    for curr_increment in range(num_increments+1):
        if increment_type == 'exponential':
            val_of_param = init_val * (increment_val)^curr_increment
        else: # increment_type == 'linear'
            val_of_param = init_val + increment_val*curr_increment
        master_tracking_dict, list_of_plotting_keys = initialize_runsim_dict(list_of_key_names, num_of_sim, num_of_timesteps, sampling_rate, param_name, val_of_param)
        
#         #save master_tracking_dict here such that it can be opened in excel!
#         f = open("outputfiles/" + datetime.datetime.now().strftime("%m_%d_%H_%M") + "total_" + "_" + str(curr_increment) + ".csv", "w")
# #        w = csv.writer(f, dialect='excel')
#         for key, val in master_tracking_dict.items():
#             if key == 'progeny_state_count':
#                 for i in range(np.size(val,0)):
# #                    w.writerow([key, val[i][3]])
#                     np.savetxt(f, val[i][3], delimiter=",", newline='\n')
#         f.close()
#         #save master_tracking_dict here such that it can be opened in excel!
#         f = open("outputfiles/" + datetime.datetime.now().strftime("%m_%d_%H_%M") + "viable_" + "_" + str(curr_increment) + ".csv", "w")
# #        w = csv.writer(f,  dialect='excel')
#         for key, val in master_tracking_dict.items():
#             if key == 'num_viable_virions':
#                 for i in range(np.size(val,0)):
# #                    w.writerow([key, val[i]])
#                     np.savetxt(f, val[i], delimiter=",")
#         f.close()        
                
        if export_raw_data_no_plot:
            for key_name in list_of_key_names:
                write_key(master_tracking_dict, key_name, output_data_label, group_by_row_num, param_name, val_of_param)
        else:
            # Plots from different const-param cross-sections should overlay on top of ea. other
            # since plots with the same title will be associated with the same figure
            if curr_increment != num_increments:
                plot_tracking_dictionary(master_tracking_dict, num_of_timesteps, list_of_plotting_keys, type_of_plot, inc_of_interest=curr_increment, num_of_const_param_sections=num_increments+1, is_last_const_param_section=False, log_setting=log_setting_opt, sampling_rate = sampling_rate)
            else:
                plot_tracking_dictionary(master_tracking_dict, num_of_timesteps, list_of_plotting_keys, type_of_plot, inc_of_interest=curr_increment, num_of_const_param_sections=num_increments+1, is_last_const_param_section=True, log_setting=log_setting_opt, sampling_rate = sampling_rate)

def plot_tracking_dictionary(input_dict, NUM_OF_TIMESTEPS, list_of_plotting_keys, side_operation=None, inc_of_interest=0, num_of_const_param_sections=1, is_last_const_param_section=True, log_setting=True, sampling_rate=1, export_raw_data_no_plot=False):
    combined_dict = dict() # This is the dictionary that we end up plotting

    for key in list_of_plotting_keys:
        combined_dict[key.name] = []

    sample_key = list_of_plotting_keys[0].name

    if len(input_dict[sample_key]) <= 0:
        raise Exception("There are no simulations remaining in input_dict!! They have been all filtered out!")
    elif len(input_dict[sample_key]) == 1: # Only one simulation in the input_dict
        for key in combined_dict:
            combined_dict[key] = input_dict[key][0]
    else:
        mean_dict = dict()
        if side_operation == "standard deviation":
            std_dev_dict = dict()
            for key in combined_dict:
                mean_dict[key] = np.mean(input_dict[key], axis=0)
                std_dev_dict[key] = np.std(input_dict[key], axis=0)

            mean_plus = dict()
            mean_minus = dict()
            for key in combined_dict:
                mean_plus[key] = np.add(mean_dict[key], std_dev_dict[key])
                mean_minus[key] = np.subtract(mean_dict[key], std_dev_dict[key])
            for key in combined_dict:
                combined_dict[key] = np.concatenate((mean_dict[key], mean_plus[key], mean_minus[key]))
        elif side_operation == "range":
            max_dict = dict()
            min_dict = dict()
            for key in combined_dict:
                mean_dict[key] = np.mean(input_dict[key], axis=0)
                min_dict[key] = np.min(input_dict[key], axis=0)
                max_dict[key] = np.max(input_dict[key], axis=0)
            for key in combined_dict:
                combined_dict[key] = np.concatenate((mean_dict[key], max_dict[key], min_dict[key]))
        elif side_operation == "all data":
            for key in combined_dict:
                mean_dict[key] = np.mean(input_dict[key], axis=0)
            combined_dict = mean_dict
            for key in combined_dict:
                # Can also do "for simulation in input_dict[key]"
                for sim_index in range(len(input_dict[key])):
                    combined_dict[key] = np.concatenate((combined_dict[key], input_dict[key][sim_index]))
        else: # Only plot the mean; side_operation is either None is something gibberish
            for key in combined_dict:
                mean_dict[key] = np.mean(input_dict[key], axis=0)
            combined_dict = mean_dict



    record1 = Record(combined_tracking_dict=combined_dict, sampling_rate=sampling_rate)
    record1.num_of_const_param_sections = num_of_const_param_sections
    record1.inc_of_interest = inc_of_interest
    if not log_setting: # If log_setting is False
        record1.use_default_log_setting = False
    if not is_last_const_param_section:
        record1.last_const_param_section = False    
    if not export_raw_data_no_plot:
        for key in list_of_plotting_keys:
            list_of_args = [record1, NUM_OF_TIMESTEPS, key.name] + key.plotting_options
            PlotHelpers.add_time_var_plot(*list_of_args)
        if is_last_const_param_section:
            plt.show()
    return record1

def initialize_runsim_dict(list_of_key_names, NUM_OF_SIMULATIONS, NUM_OF_TIMESTEPS, SAMPLING_RATE, modified_param=None, new_val_of_param=None, save_state_mode=False, timestep_list=None, batch_label="", rtn_hist_dict=False):
    list_of_plotting_keys, list_of_dependent_keys, set_of_req_record_key_names = PlotCompilation.generate_plotting_keys(list_of_key_names)

    master_tracking_dict = dict()
    for sim_index in range(NUM_OF_SIMULATIONS):
        record1 = Record(list_of_dependent_keys, set_of_req_record_key_names, sampling_rate=SAMPLING_RATE)
        record1.sim_index = sim_index
        record1.batch_label = batch_label
        sim1 = Simulation(record1, NUM_OF_TIMESTEPS, modified_param, new_val_of_param, save_state=save_state_mode, save_state_timestep_list=timestep_list)
        sim1.run()
        if sim_index == 0:
            initial_tracking_dict = record1.variable_tracking_dict
            for key in initial_tracking_dict:
                master_tracking_dict[key] = [initial_tracking_dict[key]]
            hist_dict = record1.hist_dict
        else:
            for key in master_tracking_dict:
                master_tracking_dict[key].append(record1.variable_tracking_dict[key])
    if rtn_hist_dict:
        return master_tracking_dict, list_of_plotting_keys, hist_dict
    else:
        return master_tracking_dict, list_of_plotting_keys

def initialize_label(save_type):
    label = datetime.datetime.now().strftime("%m.%d.%y_%H.%M.%S")

    directory = os.getcwd() + "/" + save_type + "/Batch" + label + "/"

    if not os.path.exists(directory):
        os.makedirs(directory)

    return label

def write_key(input_dict, key_name, output_data_label, group_by_row_num, param_name=None, val_of_param=None):
    new_filename = key_name
    if param_name != None:
        new_filename += "_" + param_name + "_" + str(val_of_param)

    path = os.getcwd() + "/outputdata/Batch" + output_data_label + "/"
    fullpath = os.path.join(path, new_filename)

    with open(fullpath + '.csv', 'wb') as csvfile:
        writer = csv.writer(csvfile)
        list_of_simulations = input_dict[key_name]
        if group_by_row_num:
            list_of_lists = []
            if len(np.shape(list_of_simulations[0])) == 1: # key's data is either numpy array or list
                list_of_lists.append([])
            else:
                for i in range(len(list_of_simulations[0])):
                    list_of_lists.append([])
            for sim_data in list_of_simulations:
                if type(sim_data) is list:
                    list_of_lists[0].append(sim_data)
                elif len(sim_data.shape) == 1: # is numpy array, not numpy matrix
                    list_of_lists[0].append(sim_data.tolist())
                else:
                    for row_index, row_content in enumerate(sim_data):
                        list_of_lists[row_index].append(row_content.tolist())

            for data_grouped_by_row in list_of_lists:
                writer.writerows(data_grouped_by_row)
                writer.writerow("")

        else:
            for sim_data in list_of_simulations:
                if type(sim_data) is list:
                    writer.writerow(sim_data)
                elif len(sim_data.shape) == 1: # is numpy array, not numpy matrix
                    writer.writerow(sim_data.tolist())
                else: # is numpy matrix
                    writer.writerows(sim_data.tolist())
                writer.writerow("")

def create_output_data_info_file(output_data_label, group_by_row_num, num_of_sim, param_name=None, increment_type=None):
    new_filename = "batch_info.txt"
    path = os.getcwd() + "/outputdata/Batch" + output_data_label + "/"
    fullpath = os.path.join(path, new_filename)

    fileHandler = open(fullpath, 'wb')
    fileHandler.write("Starting time of batch: " + output_data_label + "\n")
    fileHandler.write("Number of simulations: " + str(num_of_sim) + "\n")
    fileHandler.write("Data in csv grouped by sim? " + str(not group_by_row_num) + "\n")
    if param_name != None:
        fileHandler.write("Parameter Name: " + param_name + "\n")
        fileHandler.write("Increment Type: " + increment_type + "\n")
    fileHandler.close()
