# -*- coding: utf-8 -*-
from kartlapsim.kart import Kart, Tire
from kartlapsim.utils.helperfunctions  import save_object
import numpy as np


kart = Kart('./data/torqueCurveMiniDyno.csv')
kart.Izz = 110*0.5**2
kart.m = 110#175
kart.vmax = 150/3.6#180/3.6
kart.hcg = 0.4
kart.l = 1.05
kart.wd = 0.4    
kart.rFinal = 71/10#24/17
kart.rTrans = np.array([1])#np.array([33/13,29/16,27/18,27/22,23/22,25/27])*(75/19) #transmission ratio (nAxle/neng)
kart.tf = 1.1#1.3 #m, track width
kart.tr = 1.3
kart.frontbrakes = False#True
kart.tlltd = 0.3

#tyre
kart.tire = Tire()

kart.mux = 1.72#1.3
kart.muy = 1.7#1.95
kart.betamax = 9 #deg, max beta for slip angle drag calc
kart.rw = 0.132 #m, wheel radius
kart.rwAyCoeff = 0.02 #rel change of rw with ay
kart.sxmax = 0.08 #long slip for max long force

#resistance
kart.rho_air = 1.2
kart.cxa = 0.05#0.5
kart.croll = 0.02#0.03

kart.ini()

save_object(kart, './data/kartMini.pickle')


"""

kart = Kart('./data/torqueCurveKz.csv')
kart.m = 175
kart.hcg = 0.4
kart.l = 1.05
kart.wd = 0.4    
kart.rFinal = 24/17
kart.rTrans = np.array([33/13,29/16,27/18,27/22,23/22,25/27])*(75/19) #transmission ratio (nAxle/neng)
kart.T = 1.3 #m, track width
kart.frontbrakes = True

#tyre
kart.mux = 1.3
kart.muy = 1.95
kart.betamax = 9 #deg, max beta for slip angle drag calc
kart.rw = 0.132 #m, wheel radius
kart.rwAyCoeff = 0.02 #relative change of rw with ay
kart.sxmax = 0.08 #long slip for max long force

#resistance
kart.rho_air = 1.2
kart.cxa = 0.05
kart.croll = 0.03

kart.ini()

save_object(kart, './data/kartKz.pickle')


""" 





