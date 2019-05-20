# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 11:57:22 2019

@author: Sven van Koutrik
"""
import numpy as np
import matplotlib.pyplot as plt
from itertools import chain
from scipy.optimize import brentq #fsolve
import settings

class Lapsim:
    def __init__(self, kartobj=None,trackobj=None):
        
        self.kart=kartobj
        self.track=trackobj
        self.idxApex = None
        self.vApex = None
        self.x = State() #state obj
        self.xvec = np.zeros([len(vars(self.x)),1])
        self.xvecchans = vars(self.x).keys()
        self.laptime = None
        
    def calcapexes(self):
        #todo:exception if no kart or track defined
        vApexSol  = np.sqrt(self.kart.aymax/ np.abs(self.track.curv) )
        dvApexSol = np.gradient(vApexSol,axis=0)*vApexSol
        
        axApexSol = np.zeros_like(vApexSol)
        for i in range(len(vApexSol)):
            self.x.ay=self.kart.aymax
            self.x.v=vApexSol[i]
            axApexSol[i] = self.calcDrag()/self.kart.m #replace by calcAx and the actual ay in the more general case(max ay at nonzero fxtyre)

        idxApex = np.flatnonzero(np.diff(np.sign(dvApexSol-axApexSol),axis=0))+1 #apex def: intersection of ax at max ay with gradient of v for max ay
        self.vApexSol = vApexSol
        self.idxApex  = idxApex
        self.vApex    = vApexSol[idxApex]
        
        #plt.figure()
        #plt.plot(self.track.s,vApexSol)
        #plt.plot(self.track.s,self.track.v)
        #plt.plot(self.track.s[idxApex],self.vApex, 'ro')
        #plt.ylim((0, 40))
        #plt.show()

    def calcfwdbwd(self):
        idxApexSlowest=np.argmin(self.vApex)
        idxSlowestPoint=self.idxApex[idxApexSlowest]

        vfwd = np.zeros_like(self.track.s)
        vbwd = np.zeros_like(vfwd) 
        xvecfwd = np.zeros( (len(vars(self.x)),len(self.track.s)) )
        xvecbwd = np.zeros_like(xvecfwd)
        vfwd[idxSlowestPoint] = self.vApex[idxApexSlowest]
        vbwd[idxSlowestPoint] = self.vApex[idxApexSlowest]
        
        for i in chain( range(idxSlowestPoint+1,len(self.track.s)) , range(idxSlowestPoint)): 
            j= (2*idxSlowestPoint-i) % len(self.track.s) #reverse order bwd
            jn = (j+1) % len(self.track.s) #to make it circular
            
            vfwd[i]=self.QSSstep(vfwd[i-1],i,i-1,1)
            xvecfwd[:,i]= np.transpose(self.xvec)  #log state vector
            vbwd[j]=self.QSSstep(vbwd[jn],j,jn,-1)
            xvecbwd[:,j]= np.transpose(self.xvec)  #log state vector
                     
        self.vfwd=vfwd
        self.vbwd=vbwd
        
        self.vqss = np.minimum(vfwd,vbwd)  #pick correct speed
        self.xsol = np.where(vfwd<vbwd, np.transpose(xvecfwd) , np.transpose(xvecbwd) ) #pick correct statevector
        self.laptime = np.sum(self.track.ds / self.vqss)
        
    def QSSstep(self,vprev,idx,idxprev,direction):
        axprev = (self.calcAx(direction,vprev,idxprev) + self.track.axIncl[idx])*direction #include inclination
        vnew = vprev + (axprev/vprev)*(self.track.ds[idxprev])
        vnew = min(vnew,self.vApexSol[idx],self.kart.vmax)
        return vnew          

    def calcAx(self,direction,v,idx): #does not handle arrays because of brentq
        #direction==1: accel. direction==-1: braking
        
        #calc 'state' from v and curv
        self.x.direction = direction
        self.x.v   = v
        self.x.ay  = v**2*abs(self.track.curv[idx])
        self.x.rw  = self.rwdyn()
        self.x.fac = (direction==1)*self.kart.fxmax  +  (direction==-1)*self.kart.fxmin
        self.x.fxtyres = self.frictionEllipse()  #call friction circle(dep on direction
        self.x.fxdrag  = self.calcDrag()
        self.x.vow = v*(1+0.5*self.kart.T*abs(self.track.curv[idx])) #outside wheel 
               
        #solve PT force (unknown is RPM because of slip)
        nmin = self.x.v/self.x.rw
        nmax = self.x.vow*(1+self.kart.sxmax)/self.x.rw
        self.x.nAxle = brentq( self.errornAxle ,nmin,nmax) 
        #self.x.nAxle =  v/self.x.rw #assume no slip
        
        self.x.fxPT = self.calcFxPT(self.x.nAxle)
        ax = (np.minimum(self.x.fxtyres,self.x.fxPT) - self.x.fxdrag)/self.kart.m #add various drag components
        
        statevars = vars(self.x)
        np.put(self.xvec, range(len(statevars)) , list(statevars.values()) )
        
        return ax
    

    def frictionEllipse(self): #handles & returns arrays
        return np.sqrt(1- np.minimum((self.x.ay/self.kart.aymax)**2,1) )*self.x.fac #fx
    
    def rwdyn(self): return  self.kart.rw*(1 - self.kart.rwAyCoeff*abs(self.x.ay/self.kart.aymax))

    def calcDrag(self): #handles & returns arrays
        #slip angle drag
        beta=np.radians(self.kart.betamax*np.arcsin( np.minimum( np.abs(self.x.ay/self.kart.aymax) ,1) )/(np.pi/2)) #slip angle
        fxdragslip = self.kart.m*np.abs(self.x.ay)*np.sin(beta)
        
        #aero & rolling resistance
        fxdragroll = self.kart.m*settings.g*self.kart.croll 
        fxdragaero = 0.5*self.kart.rho_air*self.kart.cxa*self.x.v**2
        
        return fxdragslip + fxdragroll + fxdragaero
    
    
    def errornAxle(self,nAxle):
        fx      = self.calcFxPT(nAxle) #first solve what would fx be at this rpm?
        sx      = self.calcSxTyres(fx)  #then solve what would slip be at this fx        
        nAxlesol = self.x.vow*(1+sx)/self.x.rw #corresponding nAxle
        return (nAxle - nAxlesol)
    
    def calcFxPT(self,nAxle):    
        T = np.interp(nAxle,self.kart.nAxlePT,self.kart.TAxlePT)
        return T/self.x.rw #return fxPT
    
    def calcSxTyres(self,fx): 
        #basic inverse tyre model (TmEasy combined, sin(quad slip arg))
        fxn = np.minimum(fx / self.kart.fxmax,1)
        fyn = np.minimum(self.x.ay / self.kart.aymax,1)
        fn = np.minimum(np.sqrt(fxn**2 + fyn**2),1)
        sxn = fxn * np.arcsin(fn)/(np.pi/2)/fn        
        return sxn*self.kart.sxmax
        
    def plotGGV(self): #not working -> make dummy track and use idx 1 all the time
        v  = np.linspace(1, self.vmax, 100)
        ay = np.transpose(np.linspace(0, self.aymax, 100))
        
        axpos = np.zeros( (len(v),len(ay)) )
        axneg = np.zeros( (len(v),len(ay)) )
        plt.figure()
        for i in range(len(v)):
            for j in range(len(ay)):
                
                axpos[i,j] = self.calcAx(1,v[i], ay[j])
                axneg[i,j] = self.calcAx(-1,v[i], ay[j])
            plt.plot(ay,axpos[i,:])
            plt.plot(ay,axneg[i,:])
        plt.show() 

      
class State:
    def __init__(self): 
        self.direction = None
        self.v = None
        self.ay = None
        self.rw = None   
        self.fac = None
        self.fxtyres = None
        self.fxPT = None
        self.fxdrag = None
        self.vow = None
        self.nAxle = None
        
        
        #if lsobj:
            #self.ls=lsobj #couple to lapsim obj
        