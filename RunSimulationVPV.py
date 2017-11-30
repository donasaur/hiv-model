import mainaux.SimHelpers as SimHelpers

# Run RunSimulationVPV to plot multiple simulations with varying parameter value

# Number of simulations to run per constant param value (section)
NUM_OF_SIMULATIONS = 2

# Number of timesteps to run for all plots
NUM_OF_TIMESTEPS = 1000

# Parameter that determines the increment of timestep you care about, starting from t = 0
# Default is 1
SAMPLING_RATE = 1

# Name of the parameter whose value will be varied (See parameters.csv for table of param names)
PARAM_NAME = 'THRESH_TAT_FEEDBACK'

# Initial value of the parameter (See parameters.csv for suggested initial values)
# Generally, you want to set INIT_VAL to values listed in parameters.csv,
# which are obtained from the literature or fitted so that the simulation outputs
# data consistent with the literature
INIT_VAL = 0.75

# Value at which we increment the parameter value by
INCREMENT_VAL = 0.0625

# Options are 'linear', 'exponential'
#   'linear' : init_val + increment_val*curr_increment
#   'exponential' : init_val * (increment_val)^curr_increment
#   where curr_increment ranges between 0 and num_increments (inclusive)
INCREMENT_TYPE = 'linear'

# How many times should we increment value of param past INIT_VAL
# A total of (NUM_OF_INCREMENTS + 1) distinct parameter values are considered
NUM_OF_INCREMENTS = 1

# Export raw data from simulations instead of plotting them
EXPORT_RAW_DATA_NO_PLOT = False

# Group set of rows in a csv file by their row number in the
# specified key(s) instead of by their simulation number
# These group of rows are separated by newlines.
# Setting value is only relevant if EXPORT_RAW_DATA_NO_PLOT is True
GROUP_BY_ROW_NUM = True

# Determines how all the simulation data for a particular const param section will be displayed
#   options: "standard deviation", "range", "all data", None
#   Note average is always plotted for >1 sim. per const param section
# Setting value is only relevant if EXPORT_RAW_DATA_NO_PLOT is False
TYPE_OF_PLOT = "standard deviation"

# Use the default plot on log axes setting of each plot
# If False, will turn on nonlog y-axes for ALL plots
# Setting value is only relevant if EXPORT_RAW_DATA_NO_PLOT is False
USE_DEFAULT_LOG_SETTING = False

# Determines what keys will be plotted
LIST_OF_KEY_NAMES = ['proteins_nuc', 'total single spliced mRNA cyt', 'total proteins_cyt']

# Don't need to touch this line; just modify the variables above
SimHelpers.generate_increment_param_plots(NUM_OF_SIMULATIONS, NUM_OF_TIMESTEPS, PARAM_NAME, INIT_VAL, INCREMENT_VAL, INCREMENT_TYPE, NUM_OF_INCREMENTS, EXPORT_RAW_DATA_NO_PLOT, GROUP_BY_ROW_NUM, TYPE_OF_PLOT, USE_DEFAULT_LOG_SETTING, SAMPLING_RATE, LIST_OF_KEY_NAMES)
