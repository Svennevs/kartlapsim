# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 11:56:46 2019

@author: Sven van Koutrik
"""
import numpy as np
import matplotlib.pyplot as plt
import settings

class Track:
    def __init__(self, filenm=str()):
        self.s = None
        self.v = None
        self.curv = None
        self.incl = None
        self.alt = None
        
        if filenm:
            self.readLap(filenm)
        
    def readLap(self,filenm):
        datarray = np.loadtxt(filenm, delimiter=',', skiprows=17)
        with open(filenm) as fd:
            for n, line in enumerate(fd):
                if n==12: #line of channel names
                    chans =  line.split(',')
        scol    = [i for i,s in enumerate(chans) if "Distance" in s]
        vcol    = [i for i,s in enumerate(chans) if "GPS_Speed" in s]
        dpsicol = [i for i,s in enumerate(chans) if "GPS_Gyro" in s]
        inclcol = [i for i,s in enumerate(chans) if "GPS_Slope" in s]
        rpmcol  = [i for i,s in enumerate(chans) if "RPM" in s]
        
        self.s      = datarray[:,scol]*1000
        self.v      = datarray[:,vcol]/3.6
        dpsi        = np.deg2rad(datarray[:,dpsicol])
        self.curv   = dpsi/self.v
        self.ds     = np.gradient(self.s,axis=0)
        self.axIncl = np.sin( np.radians(-datarray[:,inclcol]))*settings.g
        self.rpm    = datarray[:,rpmcol]
        #self.azIncl = self.v**2*np.gradient(np.sin(np.radians(datarray[:,inclcol])),axis=0)
        
        #plt.plot(self.s,self.curv)
        #plt.show()
        
        