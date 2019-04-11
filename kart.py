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
        self.m = 110
        self.vmax = 150/3.6
        self.hcg = 0.4
        self.l = 1.05
        self.wd = 0.4        
        self.rTrans = 10/72#71 #transmission ratio (nwheel/neng)
        self.T = 1 #m, track width
        
        #tyre
        self.mux = 1.3
        self.muy = 1.72
        self.betamax = 9 #deg, max beta for slip angle drag calc
        self.rw = 0.132 #m, wheel radius
        self.rwAyCoeff = 0.02 #rel change of rw with ay
        self.sxmax = 0.08 #long slip for max long force
        
        #resistance
        self.rho_air = 1.2
        self.cxa = 0.1*0.5
        self.croll = 0.02
        
        self.calcAccLims()
        if filenm:
            self.readTorqueCurve(filenm)       
        
    def readTorqueCurve(self,filenm):
        datarray      = np.loadtxt(filenm, delimiter=';', skiprows=1) 
        self.nEngcurv = datarray[:,0]*2*math.pi/60
        self.Tcurve   = datarray[:,1]
        
    def calcAccLims(self):
        self.aymax = settings.g*self.muy
        self.fxmax = self.mux*((1-self.wd)*self.l*self.m*settings.g)/(self.l-self.mux*self.hcg)
        self.fxmin = -self.mux*((1-self.wd)*self.l*self.m*settings.g)/(self.l+self.mux*self.hcg)      

