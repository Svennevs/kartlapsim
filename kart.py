# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 11:55:52 2019

@author: Sven van Koutrik
"""
import numpy as np
import settings
import math

class Kart:
    def __init__(self, filenm=str()):
        #meta
        self.Izz = 110*0.5**2
        self.m = 110#175
        self.vmax = 150/3.6#180/3.6
        self.hcg = 0.4
        self.l = 1.05
        self.wd = 0.4    
        self.rFinal = 71/10#24/17
        self.rTrans = np.array([1])#np.array([33/13,29/16,27/18,27/22,23/22,25/27])*(75/19) #transmission ratio (nAxle/neng)
        self.tf = 1.1#1.3 #m, track width
        self.tr = 1.3
        self.frontbrakes = False#True
        self.tlltd = 0.3
        
        #tyre
        self.tire = tire()
        
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
        
        self.parameterprep()
        self.calcAccLims()
        if filenm:
            self.readTorqueCurve(filenm)       

    def parameterprep(self):
        # Get lengths -> move to precalc func
        self.lr = self.l * self.wd;    #[m] distance cog - front axle
        self.lf = self.l - self.lr;  	#[m] distance cog - rear axle
        
        
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
    

class tire:
    #tmeasy tire parameters
    def __init__(self, filenm=str()):    
        
        self.Fz_L  = 2000;              # [N] 'Low' vertical tire load for normalization
        self.F_M_L = self.Fz_L*1.5;     # Peak friction coefficient at 'low' vertical load 
        self.F_M_H = self.Fz_L*3;       # Peak friction coefficient at 'double' nominal vertical load
        self.F_S_L = self.F_M_L*0.9;    # Full sliding friction coefficient at 'low' vertical load 
        self.F_S_H = self.F_M_H*0.9;    # Full sliding friction coefficient at 'high' vertical load
        
        self.s_M_L = .07;              # slip for peak force at 'low' vertical load 
        self.s_M_H = .07;              # slip for peak force at 'high' vertical load
        self.s_S_L = .3;               # slip defining full sliding at 'low' vertical load
        self.s_S_H = .3;               # slip defining full sliding at 'high' vertical load
        
        self.dF_0_L = self.F_M_L / self.s_M_L * 2; # parameter defining initial stiffness at 'low' vertical load
        self.dF_0_H = self.F_M_H / self.s_M_H * 2; # parameter defining initial stiffness at 'high' vertical load
        
  
    
