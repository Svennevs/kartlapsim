import matplotlib.pyplot as plt

from kart import Kart
from track import Track
from lapsim import Lapsim


mini = Kart()
mini.plotGGV()
mariembourg = Track('mariembourg_mini.csv')
ls1 = Lapsim(mini,mariembourg)
ls1.calcapexes()
ls1.calcfwdbwd()


#plt.plot(ls1.track.s[ls1.idxApex],ls1.vApex*3.6, 'ro',label='apex')
#plt.plot(ls1.track.s,ls1.vApexSol*3.6,label='apex solution')
plt.plot(ls1.track.s,ls1.track.v*3.6,label='Data Rene')
plt.plot(ls1.track.s,ls1.vqss*3.6,label='simulatie')
plt.ylim((0, 120))
plt.ylabel("speed [kph]")
plt.xlabel("distance [m]")
plt.legend(loc='lower right')
plt.show()

import numpy as np
laptime = np.sum(mariembourg.ds / ls1.vqss)
print(laptime)

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



