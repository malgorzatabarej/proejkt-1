# -*- coding: utf-8 -*-
"""
Created on Sun Apr 23 21:11:54 2023

@author: USER
"""
import math
import numpy as np 
from math import *

class Transformacje:
       
    def __init__(self, model: str = "WGS84"):
        if model == "WGS84":
            self.a = 6378137.000
            self.b = 6356752.31424518
        elif model == "GRS80":
            self.a = 6378137.0
            self.b = 6356752.31414036
        else:
            raise NotImplementedError("nieobsługiwana elipsoida")
        self.flat = (self.a - self.b)/self.a
        self.e2 = sqrt(2 * self.flat - self.flat ** 2)
        self.ep2 = (2 * self.flat - self.flat ** 2)
        
    """Konwersja stopni na stopnie, minuty, sekundy"""
    def dms(x):     
        znak = ' '
        if x<0:
            znak = '-'
            x = abs(x)
        x = x * 180/np.pi
        d = int(x)
        m = int((x - d) * 60)
        s = (x - d - m / 60) * 3600
        return(znak, "%3d %2d %7.5f" %(d, m, s))
        
    """
    Tranformacja współrzędnych geocentrycznych XYZ na współrzędne elipsoidalne fi, lambda, h
    """
    def XYZ2flh(self, X, Y, Z, output = "dms"):
        """Zastosowano algorytm Hirvonena, transformujący współrzędne prostokątne na współrzędne elipsoidalne. W procesie iteracyjnym, uzyskujemy dokładne wyniki"""
        p = np.sqrt(X**2+Y**2)
        f = np.arctan(Z/(p*(1-self.ep2)))
        while True:
            N = self.a / np.sqrt(1 - self.ep2 * sin(p)**2 )
            h=(p/np.cos(f))-N
            fp=f
            f=np.arctan(Z/(p*(1-e2*N/(N+h))))
            if abs(fp-f)<(0.000001/206265):
                break
        l=np.arctan2(Y,X)
        if output == "dms":
            f = self.dms(degrees(f))
            l = self.dms(degrees(l))
            return(f,l,h)
        else:
            raise NotImplementedError("nieobsługiwana elipsoida")
    
        """
        Definicja Np
        """
        def Np(self,f):
            N=self.a/np.sqrt(1-self.ep2*np.sin(f)**2)
            return(N)

        """
        Transformacja współrzędnych elipsoidalnych fi, lambda, h na współrzędne XYZ
        """
        def flh2XYZ(self,f,l,h):
            while True:
                N=self.Np(f,self.a,self.ep2)
                X=(N+h)*np.cos(f)*np.cos(l)
                Xp=X
                Y=(N+h)*np.cos(f)*np.sin(l)
                Z=(N*(1-self.ep2)+h)*np.sin(f)
                if abs(Xp-X)<(0.000001/206265):
                    break
            return(X,Y,Z)
        
        """
        Tranformacja współrzędnych geocentryczny do współrzędnych topocentrycznych
        """
        
        def saz2neu(self, s,alfa,z):
            alfa = deg2rad(alfa)
            z = deg2rad(z)
            n = s * sin(z) * cos(alfa)
            e = s * sin(z) * sin(alfa)
            u = s * cos(z)
        return(n, e, u)
            
print('komentarz')
        
