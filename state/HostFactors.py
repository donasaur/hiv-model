# -*- coding: utf-8 -*-
import numpy as np

#This is a type of State Class
class HostFactors(object):
    def __init__(self):
        #Initialize necessary parameters
    
        #Factors associated with Tat feedback
        self.pTEFb_nuc = 500
        self.pTEFb_nuc_init = 500
        self.Tat_pTEFb_deacetyl = 0
        self.Tat_pTEFb_acetyl = 0
