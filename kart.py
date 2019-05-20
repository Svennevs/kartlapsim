# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 11:55:52 2019

@author: Sven van Koutrik
"""
import settings
import numpy as np
import math

class Kart:
    def __init__(self, filenm=str()):
        #meta
        self.m = 110#175
        self.vmax = 150#180/3.6
        self.hcg = 0.4
        self.l = 1.05
        self.wd = 0.4    
        self.rFinal = 71/10#24/17
        self.rTrans = np.array([1])#np.array([33/13,29/16,27/18,27/22,23/22,25/27])*(75/19) #transmission ratio (nAxle/neng)
        self.T = 1#1.3 #m, track width
        self.frontbrakes = False#True
        
        #tyre
        self.mux = 1.72#1.3
        self.muy = 1.7#1.95
        self.betamax = 9 #deg, max beta for slip angle drag calc
        self.rw = 0.132 #m, wheel radius
        self.rwAyCoeff = 0.02 #rel change of rw with ay
        self.sxmax = 0.08 #long slip for max long force
        
        #resistance
        self.rho_air = 1.2
        self.cxa = 0.05#0.5
        self.croll = 0.02#0.03
        
        self.calcAccLims()
        if filenm:
            self.readTorqueCurve(filenm)       
        
    def readTorqueCurve(self,filenm):
        datarray      = np.loadtxt(filenm, delimiter=';', skiprows=1) 
        
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
        self.aymax = settings.g*self.muy
        self.fxmax = self.mux*((1-self.wd)*self.l*self.m*settings.g)/(self.l-self.mux*self.hcg)
        if self.frontbrakes:
            self.fxmin = -self.mux*self.m*settings.g
        else:
            self.fxmin = -self.mux*((1-self.wd)*self.l*self.m*settings.g)/(self.l+self.mux*self.hcg)      

