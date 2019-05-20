import matplotlib.pyplot as plt
import numpy as np
try:
    import cPickle as pickle
except ModuleNotFoundError:
    import pickle

from lapsim_utils import *
from track import Track
from lapsim import Lapsim
from kart import Kart

#kart = Kart('TorqueCurve_luchtgekoeld_testbank.csv')
#kart = Kart('TorqueCurve_KZ.csv')
#save_object(kart, 'mini.pickle')
kart = load_object('KZ.pickle')
#kart = load_object('mini.pickle')

mariembourg = Track('mariembourg_mini.csv')
mbkz = Track('mariembourg_kz.csv') #lap for comparison (kz)

ls1 = Lapsim(kart,mariembourg)
ls1.calcapexes()
ls1.calcfwdbwd()

s  = shiftdist(mbkz.s,260,13.1/14)

throttle = np.maximum(np.minimum(ls1.xsol[:,5]/ls1.xsol[:,6],1),0)

plt.figure()
#plt.plot(ls1.track.s[ls1.idxApex],ls1.vApex*3.6, 'ro',label='apex')
#plt.plot(ls1.track.s,ls1.vApexSol*3.6,label='apex solution')
#plt.plot(ls1.track.s,ls1.track.v*3.6,label='Data Rene')
plt.plot(s,mbkz.v*3.6*13.1/14,label='Data Sven')
plt.plot(ls1.track.s,ls1.vqss*3.6,label='sim')
plt.plot(ls1.track.s,throttle*40, marker='.',markersize=2,label='throttle sim')
plt.ylim((0, 150))
plt.ylabel("speed [kph]")
plt.xlabel("distance [m]")
plt.legend(loc='lower right')
plt.show()

print(ls1.laptime)

#laptime    = np.sum(mariembourg.ds / mariembourg.v)

#Final drive ratio scan for mini starts here

#print(laptime)
transoptions = np.array([65,66,67,68,69,70,71,72,73])/10
optsplt = ['10/65','10/66','10/67','10/68','10/69','10/70','10/71','10/72','10/73']

lt=np.zeros_like(transoptions)
ts=np.zeros_like(lt)
maxrpm=np.zeros_like(lt)
minrpm=np.zeros_like(lt)
   
plt.figure()
for i,tr in enumerate(transoptions):
    mini = load_object('mini.pickle')#Kart('TorqueCurve_luchtgekoeld_testbank.csv')
    mini.rFinal = tr
    mini.readTorqueCurve('TorqueCurve_luchtgekoeld_testbank.csv')
    ls1 = Lapsim(mini,mariembourg)
    ls1.calcapexes()
    ls1.calcfwdbwd()
    lt[i]=ls1.laptime
    ts[i]=max(ls1.vqss)
    maxrpm[i]=max( ls1.xsol[:,9]*60/(2*np.pi) )
    plt.plot(ls1.track.s,ls1.vqss*3.6,label=optsplt[i])
    
plt.ylim((30, 120))
plt.ylabel("speed [kph]")
plt.xlabel("distance [m]")
plt.legend(loc='lower right')
plt.show()


plt.figure()
plt.plot(optsplt,lt, 'ro')
plt.ylabel("lap time [s]")
plt.xlabel("sprocket ratio")
plt.show()    

plt.figure()
plt.plot(optsplt,ts*3.6, 'ro')
plt.ylabel("top speed [kph]")
plt.xlabel("sprocket ratio")
plt.show()

plt.figure()
plt.plot(optsplt,maxrpm, 'ro')
plt.ylabel("max rpm")
plt.xlabel("sprocket ratio")
plt.show()


"""
plt.figure()
plt.plot(mariembourg.s,mariembourg.rpm / (mariembourg.v/ls1.vqss))
plt.plot(mariembourg.s, ls1.xsol[:,8]*60/(2*np.pi)  )
plt.show()
"""

"""
#postprocess / validation
v=ls1.vqss
ax=np.gradient(v, axis=0)*v
ay = v**2*mariembourg.curv
fx=mini.m*ax + mini.calcDrag(ay,v) - ls1.track.axIncl

#variable wheel radius
v=ls1.vqss
ay  = v**2*abs(self.ls.track.curv[idx])
rw  = mini.rw*(1 - mini.rwAyCoeff*abs(ay/mini.aymax))
vow = v*(1+0.5*mini.T*abs(mariembourg.curv)) #outside wheel 

nEngsol = vow*(1+sx)/self.kart.rTrans/self.x.rw       
        
ls1.x.calcState(1,v)

self.x.vow*(1+sx)/self.kart.rTrans/self.x.rw 

rw = mini.rw *(1 - 0.02*abs(ay/mini.aymax))

sx       = mini.calcSxTyres(fx,ay)  #then solve what would slip be at this fx        

nEng     = v*(1+sx)/mini.rTrans/rw
nEngcurv = v*(1+0.5*abs(mariembourg.curv))*(1+sx)/mini.rTrans/rw

fx1 = mini.calcFxEng(nEng)
fx2 = mini.calcFxEng(nEng)*(1-sx)


plt.figure()
plt.plot(mariembourg.s,mariembourg.rpm)
plt.plot(mariembourg.s, nEngcurv*60/(2*np.pi) *mariembourg.v/ls1.vqss )
plt.show()

#plt.plot(fx1)
#plt.plot(fx2)
#plt.plot(sx)



nEngsx=np.zeros_like(ls1.vqss)
nEng=np.zeros_like(ls1.vqss)
errNEng=np.zeros_like(ls1.vqss)

v=25
ay=1

nmin = v/mini.rTrans/mini.rw
nmax = v*(1+mini.sxmax)/mini.rTrans/mini.rw


nEngsx = brentq( mini.errorNEng ,nmin ,nmax ,args=(v,ay))
nEngsxi = fsolve( mini.errorNEng ,nmin, args=(v,ay))


fx      = mini.calcFxEng(nEngsx) #first solve what would fx be at this rpm?
sx      = mini.calcSxTyres(fx,ay)  #then solve what would slip be at this fx        
nEngsol = v*(1+sx)/mini.rTrans/mini.rw #corresponding nEng
(nEngsx - nEngsol)



vwhl   = nEngsx*mini.rTrans*mini.rw
"""
""""
nEngsx[i] = fsolve( mini.errorNEng ,0, args=(v[i],ay[i]))
"""


"""

for i in range(len(ls1.vqss)):
    nEngsx[i] = fsolve( mini.errorNEng ,0, args=(v[i],ay[i]))
    nEng[i] = v[i]/mini.rTrans/mini.rw

#plt.plot(mariembourg.s,nEngsx)
#plt.plot(mariembourg.s,nEng)
plt.plot(mariembourg.s,nEngsx/nEng-1)
plt.plot(mariembourg.s,sx)
plt.plot()
"""

"""


import numpy as np



plt.plot(mariembourg.s,mariembourg.azIncl)
plt.plot()

"""

"""
nEng = mariembourg.v/mini.rw/mini.rTrans
T = np.interp(nEng,mini.nEngcurv,mini.Tcurve)

plt.plot(mariembourg.s,nEng*60/(2*math.pi))
#plt.plot(mariembourg.s,T)
plt.plot()
"""

#next step: make suitable for power curve, drivetrain loss, etc
#display power limited sections
#model steering drag due to no diff (load transfer model needed?)
#GIT
#user interface
#fit parameters


"""
v=np.linspace(1,100/3.6,100)
F=np.minimum(mini.P/v,270)
plt.plot(v*3.6,F)
plt.show()
"""


"""
s = mariembourg.s
v = mariembourg.v

vnew=np.zeros_like(s)
vnew[100] = 30
for i in range(100+1,len(s)):
    vnew[i]=vnew[i-1]+1
    
        
for i in range(idxSlowestPoint,len(self.track.s)):
            vsim[i]=self._calcfwd(i)
        

#next step is fwd, backward calculation including drag
for i in range(10):
    print(i)
"""
    
"""
plt.plot(mariembourg.s,mariembourg.v)
plt.show()



plt.plot(mariembourg.s,1/mariembourg.curv)
plt.ylim((-100,100))
plt.show()

"""



