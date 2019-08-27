# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 11:55:52 2019

@author: Sven van Koutrik
"""
import kartlapsim.utils.constants as constants
import numpy as np
import math


class Kart:    
    
    def __init__(self, filenm):
        #meta
        self.m = 110
        self.vmax = 150
        self.hcg = 0.4
        self.l = 1.05
        self.wd = 0.4    
        self._rFinal = 71/10
        self.rTrans = np.array([1])
        self.T = 1 #m, track width
        self.frontbrakes = False
        
        #tyre
        self.mux = 1.72
        self.muy = 1.7
        self.betamax = 9 #deg, max beta for slip angle drag calc
        self.rw = 0.132 #m, wheel radius
        self.rwAyCoeff = 0.02 #rel change of rw with ay
        self.sxmax = 0.08 #long slip for max long force
        
        #resistance
        self.rho_air = 1.2
        self.cxa = 0.05
        self.croll = 0.02
        self.torqueCurveFile = filenm
    
    def ini(self):
        self.processTorqueCurve()
        self.calcAccLims()
     
    def processTorqueCurve(self):
        
        datarray      = np.loadtxt(self.torqueCurveFile, delimiter=';', skiprows=1) 
        
        nEng = datarray[:,0]*2*math.pi/60 #rad/s
        tEng = datarray[:,1]
        
        if len(self.rTrans)==1: #single gear
            self.nAxlePT   = nEng/self.rFinal/self.rTrans
            self.TAxlePT   = tEng*self.rFinal*self.rTrans
        else:    #multiple gears
            self.nAxlePT   = np.linspace(0,self.vmax/self.rw) #rad/s, axle speed range
            TAxleGear = []
            for i,r in enumerate(self.rTrans): #calc axle torque over speed range for each gear
                nEngPT = self.nAxlePT*self.rFinal*r
                TAxleGear.append( np.interp(nEngPT,nEng,tEng)*self.rFinal*r )    
            self.TAxlePT = np.amax( np.array(TAxleGear) ,axis=0)
                   
    def calcAccLims(self):
        self.aymax = constants.g*self.muy
        self.fxmax = self.mux*((1-self.wd)*self.l*self.m*constants.g)/(self.l-self.mux*self.hcg)
        if self.frontbrakes:
            self.fxmin = -self.mux*self.m*constants.g
        else:
            self.fxmin = -self.mux*((1-self.wd)*self.l*self.m*constants.g)/(self.l+self.mux*self.hcg)      

