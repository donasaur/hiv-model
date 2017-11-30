import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as color

def add_time_var_plot(record, max_timesteps, key, legend_colors, legend_names, plot_title, y_axis, log='off', show_legend=True, show_plot=True):
    sub_args_of_interest = [max_timesteps, legend_colors, legend_names, y_axis, log]
    if key not in record.plot_args_dict:
        record.plot_args_dict[key] = sub_args_of_interest
    else:
        record.plot_args_dict[key] = record.plot_args_dict[key] + sub_args_of_interest
    if key not in record.plot_args_dictVPV and record.num_of_const_param_sections != 1:
        record.plot_args_dictVPV[key] = sub_args_of_interest + [plot_title]
        _generate_VPV_subplots(key, record)
    if not show_plot or record.num_of_const_param_sections != 1:
        return
    plt.figure(plot_title)
    t = np.arange(0, max_timesteps, record.sampling_rate)
    orig_legend_length = len(legend_colors)
    orig_num_rows = orig_legend_length
    num_of_lines = record.get_recorded_data(key).shape[0]/orig_num_rows
    legend_colors = legend_colors * num_of_lines
    plt.gca().set_color_cycle(legend_colors)
    value = record.get_recorded_data(key)
    try:
        if log == 'on' and np.any(value>0) and record.use_default_log_setting:
            plt.yscale('log')
        lines = plt.plot(t, value.T[:,:orig_legend_length])
    except ValueError:
        plt.yscale('linear')
        lines = plt.plot(t, value.T[:,:orig_legend_length])            
    if num_of_lines != 1: # If data regarding more than 1 simulation is plotted
        plt.setp(lines, 'linewidth', 2.0)
        plt.plot(t, value.T[:, orig_legend_length:])
        if num_of_lines == 3: # If considering mean + std. dev, or mean + range
            for i in range(orig_legend_length): # i is an index
                plt.fill_between(t, value.T[:,i], value.T[:,i+orig_legend_length], facecolor=legend_colors[i], alpha=0.5)
                plt.fill_between(t, value.T[:,i], value.T[:,i+2*orig_legend_length], facecolor=legend_colors[i], alpha=0.5)
    plt.title(plot_title)
    plt.xlabel('Time(min)')
    plt.ylabel(y_axis)
    if len(legend_names)>0 and show_legend:
        plt.legend(legend_names).draggable()

# Only call this function after all the other plots have been generated
# That is, make sure subplot_helper is on a line that is after the
# function call that generates the original plots
def generate_subplots(list_of_keys, plot_title, record):
    if record.num_of_const_param_sections != 1:
        return
    plt.figure(plot_title)
    plot_args_dict = record.plot_args_dict
    num_of_subplots = 0
    num_of_subplot_specs = 5
    for key in list_of_keys:
        num_of_subplots += len(plot_args_dict[key])/num_of_subplot_specs
    subplot_index = 1
    for i in range(len(list_of_keys)):
        key = list_of_keys[i]
        sub_args = plot_args_dict[key]
        num_of_plots_for_key = len(sub_args)/num_of_subplot_specs

        for j in range(num_of_plots_for_key):
            max_timesteps = sub_args[0 + num_of_subplot_specs*j]
            legend_colors = sub_args[1 + num_of_subplot_specs*j]
            legend_names = sub_args[2 + num_of_subplot_specs*j]
            y_axis = sub_args[3 + num_of_subplot_specs*j]
            log = sub_args[4 + num_of_subplot_specs*j]

            # Creation of a new subplot
            plt.subplot(num_of_subplots, 1, subplot_index)
            orig_legend_length = len(legend_colors)
            orig_num_rows = orig_legend_length
            num_of_lines = record.get_recorded_data(key).shape[0]/orig_num_rows
            legend_colors = legend_colors * num_of_lines
            plt.gca().set_color_cycle(legend_colors)

            #Subplot specific display
            #gca() stands for 'get current axis'
            ax = plt.gca()
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

            t = np.arange(0, max_timesteps, record.sampling_rate)
            value = record.get_recorded_data(key)
            try:          
                if log == 'on' and np.any(value>0) and record.use_default_log_setting:
                    plt.yscale('log')
                lines = plt.plot(t, value.T[:,:orig_legend_length])
            except ValueError:
                plt.yscale('linear')
                lines = plt.plot(t, value.T[:,:orig_legend_length])                    
            if num_of_lines == 3: # If data regarding more than 1 simulation is plotted
                plt.setp(lines, 'linewidth', 2.0)
                plt.plot(t, value.T[:, orig_legend_length:])
                if num_of_lines == 3: # If considering mean + std. dev, or mean + range
                    for k in range(orig_legend_length): # i is an index
                        plt.fill_between(t, value.T[:,k], value.T[:,k+orig_legend_length], facecolor=legend_colors[k], alpha=0.5)
                        plt.fill_between(t, value.T[:,k], value.T[:,k+2*orig_legend_length], facecolor=legend_colors[k], alpha=0.5)                        
            if subplot_index == 1:
                plt.title(plot_title)
            subplot_index += 1 # Increment subplot_index for creation of a future subplot
            plt.ylabel(y_axis)               
    plt.xlabel('Time(min)')

def plot_histogram(hist_dict, list_of_hist_keys, plot_title, x_axis_label, y_axis_label, legend_colors, bins=None):

    plt.figure(plot_title)
    plt.title(plot_title)
    plt.xlabel(x_axis_label)
    plt.ylabel(y_axis_label)

    for i, key in enumerate(list_of_hist_keys):
        sim_data = hist_dict[key]

        if type(sim_data) is list or len(sim_data.shape) == 1:
            # work with the last timestep
            if np.size(sim_data) > 0:
                if bins == None:
                    bins = max(sim_data)
                plt.hist(sim_data, bins, alpha=0.5, color=legend_colors[i])

# A helper function for generating subplots in VPV mode
def _generate_VPV_subplots(key, record):
    sub_args = record.plot_args_dictVPV[key]

    # [max_timesteps, legend_colors, legend_names, y_axis, log, plot_title]
    max_timesteps = sub_args[0]
    legend_colors = sub_args[1]
    legend_names = sub_args[2]
    y_axis = sub_args[3]
    log = sub_args[4]
    plot_title = sub_args[5]

    plt.figure(plot_title)
    t = np.arange(0, max_timesteps, record.sampling_rate)

    num_of_subplots = len(legend_colors)

    num_of_lines = record.get_recorded_data(key).shape[0]/num_of_subplots

    for i in range(num_of_subplots):
        plt.subplot(num_of_subplots, 1, i+1)
        if record.inc_of_interest == 0:
            list_of_colors = []
            for cross_section_num in range(record.num_of_const_param_sections):
                cross_sect_tuple = _color_selector(cross_section_num, record.num_of_const_param_sections)
                list_of_colors += [cross_sect_tuple] * num_of_lines
            plt.gca().set_color_cycle(list_of_colors)
        if record.last_const_param_section:
            ax = plt.gca()
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

        t = np.arange(0, max_timesteps, record.sampling_rate)
        value = record.get_recorded_data(key)[i, :]
        value = np.array((value,))
        for j in range(num_of_lines):
            if j != 0:
                value = np.concatenate((value, (record.get_recorded_data(key)[i+j*num_of_subplots, :],)))
        try:
            if log == 'on' and value.any() and record.inc_of_interest == 0 and record.use_default_log_setting:
                plt.yscale('log')
            lines = plt.plot(t, value.T[:, 0])
        except ValueError:
            plt.yscale('linear')
            lines = plt.plot(t, value.T[:, 0])
        if num_of_lines > 1:
            plt.setp(lines, 'linewidth', 2.0)
            plt.plot(t, value.T[:,1:])
            if num_of_lines == 3:
                color_sect_tuple = _color_selector(record.inc_of_interest, record.num_of_const_param_sections)
                plt.fill_between(t, value.T[:, 0], value.T[:, 1], facecolor=color_sect_tuple, alpha=0.5)
                plt.fill_between(t, value.T[:, 0], value.T[:, 2], facecolor=color_sect_tuple, alpha=0.5)
        min_num = 0
        max_num = num_of_subplots - 1
        mid_num = (min_num + max_num)/2
        if i == 0 and record.last_const_param_section:
            plt.title(plot_title)
        elif i == mid_num and record.last_const_param_section:
            plt.ylabel(y_axis)
    if record.last_const_param_section:
        plt.xlabel('Time(min)')

# Returns an RGB tuple depending on what cross section number is passed in
def _color_selector(cross_section_num, num_of_sections):
    if num_of_sections % 2 == 1:
        if cross_section_num < num_of_sections/2:
            last_param_section_in_lower_half = num_of_sections/2
            slope = float(1)/(last_param_section_in_lower_half)
            num_of_interest = slope * cross_section_num
            return (num_of_interest, 0, 0)
        else:
            first_param_section_in_upper_half = num_of_sections/2
            last_param_section_in_upper_half = num_of_sections - 1
            slope = .7/(last_param_section_in_upper_half - first_param_section_in_upper_half)
            num_of_interest = slope * (cross_section_num - last_param_section_in_upper_half) + .7
            return (1, num_of_interest, num_of_interest)
    else: # num_of_sections is even
        if cross_section_num < num_of_sections/2:
            last_param_section_in_lower_half = num_of_sections/2 - 1
            slope = float(1)/(last_param_section_in_lower_half + 1) # will span from black to (red), where (red) is a color that is not quite red yet
            num_of_interest = slope * cross_section_num
            return (num_of_interest, 0, 0)
        else:
            first_param_section_in_upper_half = num_of_sections/2
            last_param_section_in_upper_half = num_of_sections - 1
            if num_of_sections == 2:
                slope = .7
            else:
                slope = .7/(last_param_section_in_upper_half - first_param_section_in_upper_half)
            num_of_interest = slope * (cross_section_num - last_param_section_in_upper_half) + .7
            return (1, num_of_interest, num_of_interest)
