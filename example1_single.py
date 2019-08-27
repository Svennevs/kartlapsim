import matplotlib.pyplot as plt
import numpy as np
try:
    import cPickle as pickle
except ModuleNotFoundError:
    import pickle

from kartlapsim.utils.helperfunctions import load_object
from kartlapsim.track import Track
from kartlapsim.lapsim import Lapsim
from kartlapsim.kart import Kart


track = Track('./data/trackMariembourgMini.csv')
kart  = load_object('./data/kartMini.pickle') #create object using makeKart.py
ls1   = Lapsim(kart,track)


# ///// Forward simulation using euler integration  \\\\\\

#initial cond
X = np.array([1,0,0,0.1,1,0,0])

t = np.linspace(0,1,100)
#ls1.fwdSim(X0,U,t)

Xh = np.empty((0,7))
Xd = np.empty((0,3))
for T in t:       
    Xdot = ls1.EoM(X)
    X = X + 0.01* np.multiply( np.append(Xdot,[0,0,0,0]) , ls1.scl_x  )
    Xh = np.vstack( (Xh,X))
    Xd = np.vstack( (Xd,Xdot))
    
plt.figure()
plt.plot(Xh[:,2]) #yaw rate
plt.show()    



# ///// Experimental section: \\\\\\
"""
End summary (20.05): it is difficult to make the optimizer converge
With the right scaling etc, it will find smt when neglecting the rigid axle
when including this it cannot find solutios for high curvature.
Tried to include lateral load transfer, still does not work

"""
from scipy.optimize import fsolve
from scipy.optimize import root
from scipy.optimize import Bounds
from scipy.optimize import minimize
import time

ls1 = Lapsim(kart,track)

#start = time.time()
solution = ls1.apexSolver(0.1)
#end = time.time()
#print(end - start)
solution

ls1.out.psidot * ls1.out.v
solution.x/ls1.scl_x

x0 = np.array([3,0,0,0,3,0,0])*ls1.scl_x
sol0 = root(ls1.conEqApexPre, x0, method='broyden1',tol=1e-2)
sol0
x0 = sol0.x
bnds = Bounds( [1,0,0,0,0.01]*ls1.scl_x , [ls1.kart.vmax,3,3,0.3,1]*ls1.scl_x)
con2 = {'type': 'eq', 'fun': ls1.conEqApex}
cons = con2
opts = {'maxiter': 1000}
solution = minimize(ls1.objApex ,x0,method='SLSQP',\
                    bounds=bnds,constraints=cons,options=opts)   

x=solution.x
ls1.out.curv = 0.05
ls1.objApex(x)
ls1.conEqApex(x)
ls1.conIneqApex(x)

fxmax = ls1.calcFxPT(x[4]/ls1.scl_x[4] /ls1.kart.rw) 
ls1.out.delta*ls1.out.v/ls1.kart.l
ls1.out.psidot



deltaneutral = 1.05*0.2  # l * c
ls1.out.curv = 0.05
x0 = np.array([8,0,0,0,0])

sol = root(ls1.conEqApex, x0, method='broyden1',tol=1e-2)


# plot tyre, see if smooth

fx,fy = ls1.CalcTireForce(sx, sy, Fz)

difffx = np.diff(np.diff(fx))

plt.plot(difffx)
plt.show()


