# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 10:58:58 2019
@author: Sven van Koutrik
"""

import matplotlib.pyplot as plt
import numpy as np
try:
    import cPickle as pickle
except ModuleNotFoundError:
    import pickle

from kartlapsim.utils.helperfunctions import *
from kartlapsim.track import Track
from kartlapsim.lapsim import Lapsim
from kartlapsim.kart import Kart

track = Track('..trackMariembourgMini.csv')

transoptions = np.array([65,66,67,68,69,70,71,72,73])/10
optsplt = ['10/65','10/66','10/67','10/68','10/69','10/70','10/71','10/72','10/73']

lt=np.zeros_like(transoptions)
ts=np.zeros_like(lt)
maxrpm=np.zeros_like(lt)
minrpm=np.zeros_like(lt)
   
kart = load_object('../data/kartMini.pickle')

plt.figure()
for i,tr in enumerate(transoptions):
    kart.rFinal = tr
    ls1 = Lapsim(kart,track)
    ls1.simulate()
    lt[i]=ls1.laptime
    ts[i]=max(ls1.vqss)
    maxrpm[i]=max( ls1.xsol[:,9]*kart.rFinal*60/(2*np.pi) )
    plt.plot(ls1.track.s,ls1.vqss*3.6,label=optsplt[i])
    
plt.ylim((30, 120))
plt.ylabel("speed [kph]")
plt.xlabel("distance [m]")
plt.legend(loc='lower right')
plt.show()

plt.figure()
plt.plot(optsplt,lt, 'ro')
plt.ylabel("lap time [s]")
plt.xlabel("sprocket ratio")
plt.show()    

plt.figure()
plt.plot(optsplt,ts*3.6, 'ro')
plt.ylabel("top speed [kph]")
plt.xlabel("sprocket ratio")
plt.show()

plt.figure()
plt.plot(optsplt,maxrpm, 'ro')
plt.ylabel("max rpm")
plt.xlabel("sprocket ratio")
plt.show()
