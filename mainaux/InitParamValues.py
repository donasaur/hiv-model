import csv
import numpy as np

def generate_param_dict():
    param_dict = dict()
    with open('parameters.csv', 'rbU') as csvfile:
        csvreader = csv.reader(csvfile, delimiter = ',')
        for i, row in enumerate(csvreader):
            if i != 0 and int(row[3]) == 1:
                param_dict[row[0]] = modified_eval(row[1])
    return param_dict

def modified_eval(string_input):
    if string_input[0] != "[":
        return eval(string_input)
    else:
        return eval("np.array(" + string_input + ")")