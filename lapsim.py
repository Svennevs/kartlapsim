# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 11:57:22 2019

@author: Sven van Koutrik
"""
import matplotlib.pyplot as plt
from itertools import chain
#from scipy.optimize import brentq #fsolve
from scipy.optimize import root
import settings
import numpy as np
from scipy.optimize import minimize
from scipy.optimize import Bounds
from scipy.optimize import NonlinearConstraint
from scipy.optimize import BFGS
from scipy.optimize import LinearConstraint


class Lapsim:
    def __init__(self, kartobj=None,trackobj=None):
        
        self.kart=kartobj
        self.track=trackobj
        self.idxApex = None
        self.vApex = None
        self.x = State() #state obj
        self.out = EomOutput() #output struct for EoM
        
        self.xvec = np.zeros([len(vars(self.x)),1])
        self.xvecchans = vars(self.x).keys()
        self.laptime = None
        
        #scaling factor from physical -> solver
        self.scl_x       = np.array([1/30, 1/4, 1/1.5, 1/0.4, 1/30,1,1])
        self.scl_conEq   = np.array([1/15, 1/14, 1/0.1, 1/100, 1/100])      
        self.scl_conIneq = np.array([1/500])
        self.scl_J       = np.array([1/4])
        
        
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
        return T/self.kart.rw #return fxPT
    
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


    def EoM(self,X):  #consider making a separate class out of this
        #returns state derivates and some other outputs, dep on state X and control input U
        
        #states
        vx              = X[0]/self.scl_x[0];        # [m/s]   longitudinal speed
        vy              = X[1]/self.scl_x[1];        # [m/s]   lateral speed
        psidot          = X[2]/self.scl_x[2];        # [rad/s] yaw rate
        
        #control inputs (part of state vector)
        delta           = X[3]/self.scl_x[3];
        v_Axle_r        = X[4]/self.scl_x[4];
        sx_fl           = settings.eps#0#X[6]; #no front braking
        sx_fr           = settings.eps#0#X[7]; 
        
        lty             = X[5]/self.scl_x[5]*self.kart.m*settings.g           
        ltx             = X[6]/self.scl_x[6]*self.kart.m*settings.g  
        
        #r_Fz_l          = X[8];
        #r_Fz_f          = X[9];
        
        #V_wheel_fl      = (X(6,:,:)+1).*V;  % wheel speed fl from Vwhl_rel_f
        v = np.sqrt(vx**2 + vy**2)
        beta = np.arctan(vy/vx)
        
        # Calculate aerodynamic drag
        F_ad = 0.5 * self.kart.rho_air * self.kart.cxa * v**2; # [N] aerodynamic drag
        F_adx_B = F_ad*np.cos(beta) #to x body coordinate (still drag-> negative to wrt coordinate)
        F_ady_B = F_ad*np.sin(beta) #to y body coordinate (still drag-> negative to wrt coordinate)
        #assume there is no aerodynamic moment

        TLLTD = 0.3        
                
        # neglect load transfer etc for now
        Fz_fl = self.kart.wd*self.kart.m*settings.g*0.5 + 0.5*ltx - TLLTD*lty
        Fz_fr = self.kart.wd*self.kart.m*settings.g*0.5 + 0.5*ltx + TLLTD*lty
        Fz_rl = (1-self.kart.wd)*self.kart.m*settings.g*0.5 - 0.5*ltx - (1-TLLTD)*lty
        Fz_rr = (1-self.kart.wd)*self.kart.m*settings.g*0.5 - 0.5*ltx + (1-TLLTD)*lty
        
        # Calculate ground speeds at tire contact patch
        vx_fl_B =  vx - self.kart.tf * psidot;      # vehicle body coordinate system                                     
        vx_fr_B =  vx + self.kart.tf * psidot;                             
        vx_rl_B =  vx - self.kart.tr * psidot;                                          
        vx_rr_B =  vx + self.kart.tr * psidot;                             
        vy_f_B  =  vy +  self.kart.lf * psidot;      
        vy_r_B  =  vy -  self.kart.lr * psidot;   
        
        vx_fl_T  =  vx_fl_B * np.cos(delta) + vy_f_B  * np.sin(delta);    # tire coordinate system
        vx_fr_T  =  vx_fr_B * np.cos(delta) + vy_f_B  * np.sin(delta);   
        vx_rl_T  =  vx_rl_B;                   #no rw steering                     
        vx_rr_T  =  vx_rr_B;                                        
        vy_fl_T  =  vy_f_B  * np.cos(delta) - vx_fl_B * np.sin(delta);    
        vy_fr_T  =  vy_f_B  * np.cos(delta) - vx_fr_B * np.sin(delta);    
        vy_rl_T  =  vy_r_B;                                         
        vy_rr_T  =  vy_r_B;                                         
        
        # Calculate longitudinal slip
        sx_rl = -(vx_rl_T-v_Axle_r) / (v_Axle_r + settings.eps*(v_Axle_r==0));
        sx_rr = -(vx_rr_T-v_Axle_r) / (v_Axle_r + settings.eps*(v_Axle_r==0));
        
        # Calculate lateral slip
        sy_fl   = -vy_fl_T / vx_fl_T / (1 + sx_fl + settings.eps*(sx_fl==1) ) ;
        sy_fr   = -vy_fr_T / vx_fr_T / (1 + sx_fr + settings.eps*(sx_fr==1)) ;
        sy_rl   = -vy_rl_T  / vx_rl_T / (1 + sx_rl + settings.eps*(sx_rl==1)) ;
        sy_rr   = -vy_rr_T  / vx_rr_T / (1 + sx_rr  + settings.eps*(sx_rr==1)) ;
        
        
        # Tire Forces
        Fx_fl_T, Fy_fl_T = self.CalcTire_TMeasy( sx_fl, sy_fl, Fz_fl);
        Fx_fr_T, Fy_fr_T = self.CalcTire_TMeasy( sx_fr, sy_fr, Fz_fr);
        Fx_rl_T, Fy_rl_T = self.CalcTire_TMeasy( sx_rl, sy_rl, Fz_rl);
        Fx_rr_T, Fy_rr_T = self.CalcTire_TMeasy( sx_rr, sy_rr, Fz_rr);
        
        # Tire force in body coordinate system
        Fx_fl_B = np.cos(delta)*Fx_fl_T - np.sin(delta)*Fy_fl_T
        Fx_fr_B = np.cos(delta)*Fx_fr_T - np.sin(delta)*Fy_fr_T
        Fx_rl_B = Fx_rl_T
        Fx_rr_B = Fx_rr_T
        Fy_fl_B = np.cos(delta)*Fy_fl_T + np.sin(delta)*Fx_fl_T
        Fy_fr_B = np.cos(delta)*Fy_fr_T + np.sin(delta)*Fx_fr_T
        Fy_rl_B = Fy_rl_T
        Fy_rr_B = Fy_rl_T
        
        Fx_B = Fx_fl_B + Fx_fr_B + Fx_rl_B + Fx_rr_B;  # vehicle body coordinates
        Fy_B = Fy_fl_B + Fy_fr_B + Fy_rl_B + Fy_rr_B;                       
        Mz   = (Fy_fl_B+Fy_fr_B)*self.kart.lf - (Fy_rl_B+Fy_rr_B)*self.kart.lr + (Fx_fr_B-Fx_fl_B)*self.kart.tf/2 + (Fx_rr_B-Fx_rl_B)*self.kart.tr/2  
        
        #planar dynamics
        vxdot    = (1/self.kart.m)*(Fx_B - F_adx_B) + psidot*vy  # [m/s**2] x acceleration
        vydot    = (1/self.kart.m)*(Fy_B - F_ady_B) - psidot*vx 
        psiddot  = (1/self.kart.Izz)*Mz
        
        # Compose vector with state derivative
        Xdot = np.array([vxdot,vydot,psiddot]);
        
        self.out.X        = X
        self.out.vx       = vx
        self.out.vy       = vy
        self.out.psidot   = psidot
        self.out.delta    = delta
        self.out.v_Axle_r = v_Axle_r
        self.out.v        = v
        self.out.Fx_rl_T  = Fx_rl_T
        self.out.Fx_rr_T  = Fx_rr_T
        self.out.vxdot    = vxdot
        self.out.vydot    = vydot
        self.out.psiddot  = psiddot
        self.out.Fy_B     = Fy_B
        self.out.Fx_B     = Fx_B
        self.out.Fz_fl    = Fz_fl
        self.out.Fz_fr    = Fz_fr
        self.out.Fz_rl    = Fz_rl
        self.out.Fz_rr    = Fz_rr
        
        return Xdot
        
            
    #check if Xdot for this x has been calced
    #if not -> do so
    #if yes: calc objective
        
    def CalcTire_TMeasy(self,sx, sy, Fz):
        # Load relevant parameters
        Fz_L_fac = Fz/self.kart.tire.Fz_L; # normalized to lower force in which tire parameters are defined 
        # Forces: expressed with 2nd degree polynom: F = -a*Fz^2+ b*Fz + c
        F_M  = Fz_L_fac * ( 2*self.kart.tire.F_M_L  - 0.5*self.kart.tire.F_M_H  - (self.kart.tire.F_M_L  - 0.5*self.kart.tire.F_M_H)  * Fz_L_fac ); # maximum xy force possible at given Fz
        F_S  = Fz_L_fac * ( 2*self.kart.tire.F_S_L  - 0.5*self.kart.tire.F_S_H  - (self.kart.tire.F_S_L  - 0.5*self.kart.tire.F_S_H)  * Fz_L_fac ); # sliding y force at given Fz
        # Slips: linear interpolation
        s_M  = (self.kart.tire.s_M_L + (self.kart.tire.s_M_H-self.kart.tire.s_M_L) * (Fz_L_fac-1)); # slip in x needed for maximum Fy
        s_S  = (self.kart.tire.s_S_L + (self.kart.tire.s_S_H-self.kart.tire.s_S_L) * (Fz_L_fac-1)); # slip in x at which Fy is sliding force
        # Stiffness: linear interpolation
        dF_0 = self.kart.tire.dF_0_L + (self.kart.tire.dF_0_H-self.kart.tire.dF_0_L) * (Fz_L_fac-1);
        
        
        # Slip
        s = np.sqrt( sx**2 + sy**2 ); # direction independent slip
        """
        # Calculation of combined normalised tire force
        sigma   =    (s <= s_M) * s/s_M  + \
        ((s_M < s) & (s < s_S)) * (s - s_M) / (s_S - s_M);	# sigma 0 -> 1; slip at max force -> pure sliding slip
        
        F       =    (s <= s_M) * s_M * dF_0 * sigma / (1 + sigma * (dF_0* s_M/F_M  - 2 + sigma) )  + \
        ((s_M < s) & (s < s_S)) * ( F_M  - (F_M  - F_S) * sigma**2 * (3 - 2 * sigma) ) + \
        (s_S <= s) * F_S;                                                                           
        """
        F = F_M*np.tanh(s/0.05)-0.2*s
        
         
        
        # Decomposition of combined in Lateral and longitudinal force and denormalisation
        Fy = F * sy/(s + (s==0)*settings.eps); # distribute forces to x and y as slips are distributed
        Fx = F * sx/(s + (s==0)*settings.eps); # distribute forces to x and y as slips are distributed
        
       
        
        return Fx, Fy
        
    
    def objApex(self,X):
        
        if True:#not(np.all(self.out.X == X)): #check if EoM(X) was done already
            self.EoM(X)
        return -self.scl_J*(self.out.psidot)**2 
    
    def conEqApex(self,X):
        if True:#not(np.all(self.out.X == X)): #check if EoM(X) was done already
            self.EoM(X)
        
        eqLty = self.kart.hcg*self.out.Fy_B  -  0.5*self.kart.tf*(self.out.Fz_fr-self.out.Fz_fl) - 0.5*self.kart.tr*(self.out.Fz_rr-self.out.Fz_rl) 
        eqLtx = self.kart.hcg*self.out.Fx_B  -  (self.kart.lr*(self.out.Fz_rr+self.out.Fz_rl) - self.kart.lf*(self.out.Fz_fr - self.out.Fz_fl))
        
        return np.array([self.out.vydot , self.out.psiddot,\
                         (self.out.curv*self.out.v - self.out.psidot),\
                          eqLty,eqLtx])*self.scl_conEq

    def conEqApexPre(self,X):
        if True:#not(np.all(self.out.X == X)): #check if EoM(X) was done already
            self.EoM(X)
        
        return self.out.vydot , self.out.psiddot,(self.out.curv*self.out.v - self.out.psidot)*0.01, self.out.v - np.sqrt(10/self.out.curv) #presolve at ay=10m/s**2
    
    
    def conIneqApex(self,X): #outcome needs to be pos
        if True:#not(np.all(self.out.X == X)): #check if EoM(X) was done already
            self.EoM(X)
        
        fxmax = self.calcFxPT(X[4]/self.scl_x[4]/self.kart.rw) 
        ineq1   = fxmax - (self.out.Fx_rl_T + self.out.Fx_rr_T)
        return ineq1*self.scl_conIneq
    
    
    def apexSolver(self,curv):
        
        self.out.curv = curv        
        x0 = np.array([np.sqrt(10/curv),-3,0,0, np.sqrt(10/curv),0,0])*self.scl_x

        #sol0 = root(self.conEqApexPre, x0, method='broyden1',tol=1e-2)
        #x0 = sol0.x
        #theoretically can check bounds
        
       # bnds = Bounds( [1,0,0,0,0.01]*self.scl_x , [self.kart.vmax,3,3,0.3,1]*self.scl_x)
        bnds = ( (5/30,self.kart.vmax/30),(-3/4,3/4),(0/1.5,2/1.5),(0/0.4,0.5/0.4),(5/30,self.kart.vmax/30),(-1,1),(-1,1))
        #conlin = LinearConstraint([],[],[])
        #connonlin = NonlinearConstraint(self.conEqApex, 0, 0, jac='2-point', hess=BFGS() )
        con1 = {'type': 'ineq', 'fun': self.conIneqApex} 
        con2 = {'type': 'eq', 'fun': self.conEqApex}
        
        cons = ([con1,con2])
        opts = {'maxiter': 1000}

        solution = minimize(self.objApex ,x0,method='SLSQP',\
                            bounds=bnds,constraints=cons,options=opts)                  
        return solution



class EomOutput:
    def __init__(self): 
        self.X        = None
        self.vx       = None
        self.vy       = None
        self.psidot   = None
        self.delta    = None
        self.v_Axle_r = None
        self.v        = None
        self.Fx_rl_T  = None
        self.Fx_rr_T  = None
        self.vxdot    = None
        self.vydot    = None
        self.psiddot  = None
        self.Fy_B = None
        self.Fx_B = None
        self.Fz_fl= None
        self.Fz_fr= None
        self.Fz_rl= None
        self.Fz_rr= None
        
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
        