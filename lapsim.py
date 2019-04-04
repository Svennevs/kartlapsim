# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 11:57:22 2019

@author: Sven van Koutrik
"""
import numpy as np
import matplotlib.pyplot as plt
from itertools import chain

class Lapsim:
    def __init__(self, kartobj=None,trackobj=None):
        
        self.kart=kartobj
        self.track=trackobj
        self.idxApex = None
        self.vApex = None
        
    def calcapexes(self):
        #todo:exception if no kart or track defined
        vApexSol  = np.sqrt(self.kart.aymax/ np.abs(self.track.curv) )
        dvApexSol = np.gradient(vApexSol,axis=0)*vApexSol      
        axApexSol = self.kart.calcDrag(self.kart.aymax,vApexSol)/self.kart.m #replace by calcAx and the actual ay in the more general case(max ay at nonzero fxtyre)

        idxApex = np.flatnonzero(np.diff(np.sign(dvApexSol-axApexSol),axis=0))+1 #apex def: intersection of ax at max ay with gradient of v for max ay
        self.vApexSol = vApexSol
        self.idxApex  = idxApex
        self.vApex    = vApexSol[idxApex]
        
        plt.plot(self.track.s,vApexSol)
        plt.plot(self.track.s,self.track.v)
        plt.plot(self.track.s[idxApex],self.vApex, 'ro')
        plt.ylim((0, 40))
        plt.show()

    def calcfwdbwd(self):
        idxApexSlowest=np.argmin(self.vApex)
        idxSlowestPoint=self.idxApex[idxApexSlowest]

        vfwd = np.zeros_like(self.track.s)
        vbwd = np.zeros_like(self.track.s) 
        vfwd[idxSlowestPoint] = self.vApex[idxApexSlowest]
        vbwd[idxSlowestPoint] = self.vApex[idxApexSlowest]
        
        for i in chain( range(idxSlowestPoint+1,len(self.track.s)) , range(idxSlowestPoint)): 
            j= (2*idxSlowestPoint-i) % len(self.track.s) #reverse order bwd
            jn = (j+1) % len(self.track.s) #to make it circular
            vfwd[i]=self.QSSstep(vfwd[i-1],i,i-1,1)
            vbwd[j]=self.QSSstep(vbwd[jn],j,jn,-1)
                      
        self.vfwd=vfwd
        self.vbwd=vbwd
        self.vqss=np.minimum(vfwd,vbwd)
        
    def QSSstep(self,vprev,idx,idxprev,direction):
        
        ayprev = vprev**2*abs(self.track.curv[idxprev])
        axprev = (self.kart.calcAx(direction,vprev,ayprev) + self.track.axIncl[idx])*direction #include inclination
        vnew = vprev + (axprev/vprev)*(self.track.ds[idxprev])
        vnew = min(vnew,self.vApexSol[idx],self.kart.vmax)
        return vnew   
            
            
       
            
            
            
            