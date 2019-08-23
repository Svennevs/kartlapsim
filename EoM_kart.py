# -*- coding: utf-8 -*-
"""
Created on Thu May 23 12:45:45 2019

@author: Sven van Koutrik
"""

function [ X_dot, Output ] = EoM_Vehicle( X, U, C, ModelParams )
# EOM_VEHICLE evaluates the equations of motion for a two-track vehicle
# model, with aerodynamic downforce and a limited-slip differential
#
# INPUT:
# X           [dimX x N x:]  [double]    [var]   states
# U           [dimU x N x:]   [double]    [var]   control inputs
# C           [dimC x N x:]   [double]    [par]   track curvature [rad/m]
# ModelParams           [struct]    [-]     all model parameters
#
# OUTPUT:
# States_dot  [12xNx:]  [double]    [var]   distance derivatives of states
# Output      [?xNx:]   [struct]    output data used for constraint function and plots
#
#   Sven van Koutrik, 2014-10-23 


# Reassign variables

#states
vx              = X[1];        # [m/s]   longitudinal speed
vy              = X[2];        # [m/s]   lateral speed
psidot          = X[3];        # [rad/s] yaw rate

#control inputs (part of state vector)
delta           = U[1];
v_Axle_r        = U[2];
sx_fl           = 0#X[6]; #no front braking
sx_fr           = 0#X[7]; 
#r_Fz_l          = X[8];
#r_Fz_f          = X[9];

#V_wheel_fl      = (X(6,:,:)+1).*V;  % wheel speed fl from Vwhl_rel_f
v = np.sqrt(vx**2 + vy**2)
beta = np.atan(vy/vx)

# Get lengths -> move to precalc func
self.lf = self.l * self.wd;    #[m] distance cog - front axle
self.lr = self.l - self.lf;  	#[m] distance cog - rear axle

# Calculate aerodynamic drag
F_ad = 0.5 * self.rho_air * self.cxa * v**2; # [N] aerodynamic drag
F_adx_B = F_d*cos(beta) #to x body coordinate (still drag-> negative to wrt coordinate)
F_ady_B = F_d*sin(beta) #to y body coordinate (still drag-> negative to wrt coordinate)
#assume there is no aerodynamic moment

# neglect load transfer etc for now
F_z_fl = self.wd*self.m*settings.g
F_z_fr = self.wd*self.m*settings.g
F_z_rl = (1-self.wd)*self.m*settings.g
F_z_rr = (1-self.wd)*self.m*settings.g

# Calculate ground speeds at tire contact patch
vx_fl_B =  vx - self.tf * psidot;      # vehicle body coordinate system                                     
vx_fr_B =  vx + self.tf * psidot;                             
vx_rl_B =  vx - self.tr * psidot;                                          
vx_rr_B =  vx + self.tr * psidot;                             
vy_f_B  =  vy +  self.lf * psidot;      
vy_r_B  =  vy -  self.lr * psidot;   

vx_fl_T  =  vx_fl_B * np.cos(delta) + vy_f_B  * np.sin(delta);    # tire coordinate system
vx_fr_T  =  vx_fr_B * np.cos(delta) + vy_f_B  * np.sin(delta);   
vx_rl_T  =  vx_rl_B;                   #no rw steering                     
vx_rr_T  =  vx_rr_B;                                        
vy_fl_T  =  vy_f_B  * np.cos(delta) - vx_fl_B * np.sin(delta);    
vy_fr_T  =  vy_f_B  * np.cos(delta) - vx_fr_B * np.sin(delta);    
vy_rl_T  =  vy_r_B;                                         
vy_rr_T  =  vy_r_B;                                         

# Calculate longitudinal slip
sx_rl = -(vx_rl_T-v_Axle_r) / v_Axle_r*(v_Axle_r==0);
sx_rr = -(vx_rr_T-v_Axle_r) / v_Axle_r*(v_Axle_r==0);

# Calculate lateral slip
s_y_fl   = -vy_fl_T / vx_fl_T / (1 + sx_fl + settings.eps*(sx_fl==1) ) ;
s_y_fr   = -vy_fr_T / vx_fr_T / (1 + sx_fr + settings.eps*(sx_fr==1)) ;
s_y_rl   = -vy_rl_T  / vx_rl_T / (1 + sx_rl + settings.eps*(sx_rl==1)) ;
s_y_rr   = -vy_rr_T  / vx_rr_T / (1 + sx_rr  + settings.eps*(sx_rr==1)) ;

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
Mz   = (Fy_fl_B+Fy_fr_B)*self.lf - (Fy_rl_B+Fy_rr_B)*self.lr + (Fx_fr_B-Fx_fl_B)*self.tf/2 + (Fx_rr_B-Fx_rl_B)*self.tr/2  

#planar dynamics
vxdot    = (1/self.m)*(+psidot*vy + Fx_B - F_adx_B)  # [m/s**2] x acceleration
vydot    = (1/self.m)*(-psidot*vx + Fy_B - F_ady_B)
psiddot  = (1/self.Izz)*Mz

# Compose vector with state derivative
X_dot = np.array([vxdot,vydot,psiddot]);

#cost
J = v #apex solution
J = vxdot #forward/backwards

#equality constraints
psiddot = 0
vydot   = 0     #no spinning -> not entirely correct

#inequality constraints
fxPTmax = calcFxPT(vAxle/self.rw) 
Fx_rl_T + Fx_rr_T < fxPTmax


