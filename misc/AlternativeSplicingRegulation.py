# -*- coding: utf-8 -*-
import numpy as np

class AlternativeSplicingRegulation(object):
    def __init__(self):
        self.num_binding_sites = 11
    def convert_binary_to_dec(self, binary_array):
        #say you have a vector of 0s and 1s
        #you want to convert this to a decimal base index
        binary_string = ''.join(str(x) for x in binary_array)
        dec_index = int(binary_string,2)
        return [dec_index]
    def convert_dec_to_binary(self, dec_index):
        #say you have a decimal base index and you want to 
        #convert this to an np.array of binary values that indicate where
        #different regulons are bound
        #this array must have length = self.num_binding_sites 
        binary_int = int(bin(dec_index)[2:])
        binary_array = np.array([int(char) for char in str(binary_int)])
        binary_array = np.pad(binary_array, (self.num_binding_sites-np.size(binary_array),0), 'constant', constant_values=(0,0))
        
