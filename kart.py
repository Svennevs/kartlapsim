# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 11:55:52 2019

@author: Sven van Koutrik
"""
import settings
import numpy as np
import math
import matplotlib.pyplot as plt

class Kart:
    def __init__(self):
        self.m = 110
        self.P = 6.4e3
        self.mux = 1.3
        self.muy = 1.75
        self.vmax = 150/3.6
        self.hcg = 0.4
        self.l = 1.05
        self.wd = 0.4
        self.rho_air = 1.2
        self.cxa = 0.5*0.6
        self.croll = 0.03
        self.fxmaxEngine = 320 #N, max engine force
        self.betamax = 9 #deg, max beta for slip angle drag calc
        
        self.calcAccLims()
        
    def calcAccLims(self):
        self.aymax = settings.g*self.muy
        self.fxmax = self.mux*((1-self.wd)*self.l*self.m*settings.g)/(self.l-self.mux*self.hcg)
        self.fxmin = -self.mux*((1-self.wd)*self.l*self.m*settings.g)/(self.l+self.mux*self.hcg)      
        
    def frictionEllipse(self,v,ay,fac): #handles & returns arrays
        fx = np.sqrt(1- np.minimum((ay/self.aymax)**2,1) )*fac
        return fx
    
    def calcDrag(self,ay,v): #handles & returns arrays
        #slip angle drag
        beta=np.radians(self.betamax*np.arcsin( np.minimum( (ay/self.aymax)**2 ,1) )/(math.pi/2)) #slip angle
        fxdragslip = self.m*np.abs(ay)*np.sin(beta)
        
        #aero & rolling resistance
        fxdragroll=self.m*settings.g*self.croll 
        fxdragaero=0.5*self.rho_air*self.cxa*v**2
        
        return fxdragslip + fxdragroll + fxdragaero
        
    def calcAx(self,direction,v,ay): #handles & returns arrays
        #direction==1: accel. direction==-1: braking
        
        fac = (direction==1)*self.fxmax  +  (direction==-1)*self.fxmin
        fxtyres  = self.frictionEllipse(v,ay,fac)  #call friction circle(dep on direction
        fxdrag   = self.calcDrag(ay,v)
        fxengine = np.minimum(self.P/v , self.fxmaxEngine) #engine power/torque limit        
        
        ax = (np.minimum(fxtyres,fxengine) - fxdrag)/self.m #add various drag components
        return ax
        
    
    def plotGGV(self):
        v  = np.linspace(1, self.vmax, 100)
        ay = np.transpose(np.linspace(0, self.aymax, 100))
        
        axpos = np.zeros( (len(v),len(ay)) )
        axneg = np.zeros( (len(v),len(ay)) )
        for i in range(len(v)):
            axpos[i,:] = self.calcAx(1,v[i], ay)
            axneg[i,:] = self.calcAx(-1,v[i], ay)
            plt.plot(ay,axpos[i,:])
            plt.plot(ay,axneg[i,:])
        plt.show() 
