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
#kart  = load_object('./data/kartKz.pickle') #create object using makeKart.py

ls1 = Lapsim(kart,track)
ls1.simulate()

plt.figure()
plt.plot(ls1.track.s,ls1.vApexSol*3.6,label='apex solution')
plt.plot(ls1.track.s,ls1.track.v*3.6,label='Reference data')
plt.plot(ls1.track.s,ls1.vqss*3.6,label='sim')
plt.ylim((0, 150))
plt.ylabel("speed [kph]")
plt.xlabel("distance [m]")
plt.legend(loc='lower right')
plt.show()

print("simulated lap time is " , round(ls1.laptime,2))

