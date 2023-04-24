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
            
        """ 
        Transformacja fi, lambda do układu 2000
        """
        
        def sigma(self, f, a, ep2):
            A0 = 1 - ep2/4 - 3 * ep2**2/64 - 5 * ep2**3/256
            A2 = (3/8) * (ep2 + ep2**2/4 + 15 * ep2**3/128)
            A4 = (15/256) * (ep2**2 + 3 * ep2**3/4)
            A6 = 35 * ep2**3/3072
            sig = a * (A0 * f - A2 * np.sin(2 * f) + A4 * np.sin(4 * f) - A6 * np.sin(6 * f))
            return(sig)
        
        def GK2000(self, f, l, a, ep2):
            m=0.999923
            l0 = 0 
            strefa = 0
            if l >np.deg2rad(13.5) and l < np.deg2rad(16.5):
                strefa = 5
                la0 = np.deg2rad(15)
            elif l >np.deg2rad(16.5) and l < np.deg2rad(19.5):
                strefa = 6
                l0 = np.deg2rad(18)
            elif l >np.deg2rad(19.5) and l < np.deg2rad(22.5):
                strefa =7
                l0 = np.deg2rad(21)
            elif l >np.deg2rad(22.5) and l < np.deg2rad(25.5):
                strefa = 8
                l0 = np.deg2rad(24)
            else:
                return("Punkt poza strefami odwzorowawczymi układu PL-2000")        
    
        b2 = (a**2) * (1-e2)   #krotsza polos
        e2p = ( a**2 - b2 ) / b2   #drugi mimosrod elipsy
        dl = l - l0
        t = np.tan(f)
        ni = np.sqrt(e2p * (np.cos(f))**2)
        N = self.Np(f, a, ep2)
        sigma = self.sigma(f, a, ep2)
        XGK20 = sigma + ((dl**2)/2)*N*np.sin(f)*np.cos(f) * ( 1 + ((dl**2)/12)*(np.cos(f))**2 * ( 5 - (t**2)+9*(ni**2) + 4*(ni**4)     )  + ((dl**4)/360)*(np.cos(f)**4) * (61-58*(t**2)+(t**4) + 270*(ni**2) - 330*(ni**2)*(t**2))  )
        YGK20 = (dl*N* np.cos(f)) * (1+(((dl)**2/6)*(np.cos(f))**2) *(1-(t**2)+(ni**2))+((dl**4)/120)*(np.cos(f)**4)*(5-18*(t**2)+(t**4)+14*(ni**2)-58*(ni**2)*(t**2)) )
        X2000 = XGK20 * m 
        Y2000 = YGK20 * m + strefa*1000000 + 500000
        return(X2000, Y2000)
