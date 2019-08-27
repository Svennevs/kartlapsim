# -*- coding: utf-8 -*-
"""
Created on Mon May 20 15:31:13 2019

@author: Sven van Koutrik
"""
import numpy as np
try:
    import cPickle as pickle
except ModuleNotFoundError:
    import pickle
    
    
def save_object(obj, filename):
    with open(filename, 'wb') as output:  # Overwrites any existing file.
        pickle.dump(obj, output)

def load_object(filename):  
    with open(filename, 'rb') as pickle_file:
        obj = pickle.load(pickle_file)
    return obj

def shiftdist(s,sshift,stretch):
    s    = s*stretch 
    sn   = s - sshift
    smax = np.max(s)
    
    for i,y in enumerate(sn):
        if y<0:
            sn[i]=sn[i]+smax
    return sn