# -*- coding: utf-8 -*-
"""
Created on Sun Apr 23 21:11:54 2023

@author: USER
"""
import math
import numpy as np 
from math import *
import argparse

o = object()

class Transformacje:
    def __init__(self, model):
        if model == "WGS84":
            self.a = 6378137.000
            self.b = 6356752.31424518
        elif model == "GRS80":
            self.a = 6378137.0
            self.b = 6356752.31414036
        elif model == "Krasowski":
            self.a = 6378245.0
            self.b = 6356752.314245
        else:
            raise NotImplementedError("nieobsługiwana elipsoida")
        self.flat = (self.a - self.b)/self.a
        self.e2 = sqrt(2 * self.flat - self.flat ** 2)
        self.ep2 = (2 * self.flat - self.flat ** 2)
        
    """
    Tranformacja współrzędnych geocentrycznych XYZ na współrzędne elipsoidalne fi, lambda, h
    """
    
    def XYZ2flh(self, X, Y, Z, output = "dec_degree"):
        FLH = []
        
        """Zastosowano algorytm Hirvonena, transformujący współrzędne prostokątne na współrzędne elipsoidalne. W procesie iteracyjnym, uzyskujemy dokładne wyniki"""
        
        r   = np.sqrt(X**2 + Y**2)           # promień
        Fi = atan(Z / (r * (1 - self.ep2)))    # pierwsze przybliilizenie
        fi = 0
        while abs(Fi - fi) > 0.000001/206265:    
            Fi = fi
            N = self.a / sqrt(1 - self.ep2 * sin(Fi)**2)
            h = r / cos(Fi) - N
            lat = atan((Z/r) * (((1 - self.ep2 * N/(N + h))**(-1))))
        lam = atan(Y/X)
        N = self.a / sqrt(1 - self.ep2 * (sin(fi))**2);
        h = r / cos(fi) - N  
        FLH.append([degrees(Fi), degrees(lam), h])
        if output == "dec_degree":
            return FLH
        elif output == "dms":
            lat = self.deg2dms(degrees(lat))
            lon = self.deg2dms(degrees(lon))
            return f"{lat[0]:02d}:{lat[1]:02d}:{lat[2]:.2f}", f"{lon[0]:02d}:{lon[1]:02d}:{lon[2]:.2f}", f"{h:.3f}"
        else:
            raise NotImplementedError ("nieobsługiwana elipsoida")
    
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
            XYZ = []
            while True:
                N=self.Np(f)
                X=(N+h)*np.cos(f)*np.cos(l)
                Y=(N+h)*np.cos(f)*np.sin(l)
                Z=(N*(1-self.ep2)+h)*np.sin(f)
                if abs(Xp-X)<(0.000001/206265):
                    break
                XYZ.append([X, Y, Z])
            return XYZ
        
        """
        Tranformacja współrzędnych geocentryczny do współrzędnych topocentrycznych
        """
        def Rneu(f,l):
            R = np.array([[-np.sin(f)*np.cos(l), -np.sin(l), np.cos(f)*np.cos(l)],
                          [-np.sin(f)*np.sin(l), np.cos(l), np.cos(f)*np.sin(l)],
                          [np.cos(f), 0., np.sin(f)]])
            return(R)
        
    
        def XYZ2NEU(self, X, Y, Z, X0, Y0, Z0):
            neu = []
            fi, lam, h = self.XYZ2flh(X, Y, Z, output = "degrees")
            
            Rneu = self.Rneu(f,l)
            for X, Y, Z in zip(X, Y, Z):
                X_sr = [X-X0, Y-Y0, Z-Z0] 
                X_rneu = R_neu.T@X_sr
                neu.append(X_rneu.T)
                
            return neu

        """
        Tranformacja współrzędnych fi, lambda do układu 2000
        """
        
        def sigma(self, f, a, ep2):
            A0 = 1 - ep2/4 - 3 * ep2**2/64 - 5 * ep2**3/256
            A2 = (3/8) * (ep2 + ep2**2/4 + 15 * ep2**3/128)
            A4 = (15/256) * (ep2**2 + 3 * ep2**3/4)
            A6 = 35 * ep2**3/3072
            sig = a * (A0 * f - A2 * np.sin(2 * f) + A4 * np.sin(4 * f) - A6 * np.sin(6 * f))
            return(sig)
        
        def GK2000(self, f, l, a, ep2):
            XY2000 = []
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
    
        b2 = (a**2) * (1-ep2)   #krotsza polos
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
        XY2000.append([X2000, Y2000])
        return(X2000, Y2000)
    
    """
    Tranformacja współrzędnych fi, lambda do układu 1992
    """
    
    def GK1992(self, f, l, a , ep2):
        XY92 = []
        lam0 = (19*np.pi)/180
        m = 0.9993
        b2 = (a**2) * (1-ep2)   #krotsza polos
        e2p = ( a**2 - b2 ) / b2   #drugi mimosrod elipsy
        dlam = l - lam0
        t = np.tan(f)
        ni = np.sqrt(e2p * (np.cos(f))**2)
        N = self.Np(f)

        sigma = self.sigma(f, a, ep2)

        xgk = sigma + ((dlam**2)/2)*N*np.sin(f)*np.cos(f) * ( 1+ ((dlam**2)/12)*(np.cos(f))**2 * ( 5 - (t**2)+9*(ni**2) + 4*(ni**4)     )  + ((dlam**4)/360)*(np.cos(f)**4) * (61-58*(t**2)+(t**4) + 270*(ni**2) - 330*(ni**2)*(t**2))  )
        ygk = (dlam*N* np.cos(f)) * (1+(((dlam)**2/6)*(np.cos(f))**2) *(1-(t**2)+(ni**2))+((dlam**4)/120)*(np.cos(f)**4)*(5-18*(t**2)+(t**4)+14*(ni**2)-58*(ni**2)*(t**2)) )
        
        x92 = xgk*m - 5300000
        y92 = ygk*m + 500000
        
        XY92.append([x92, y92])
        
        return XYZ92 
    
"""   def plik(self, plik_wynikowy, transf: type = str):
        dane = np.genfromtxt(plik_wynikowy,delimiter = " ")
        if transf == 'XYZ2BLH':
            FLH = self.XYZ2flh(dane[:,0], dane[:,1], dane[:,2])
            np.savetxt(f"plik_wynikowy_{transformacja}_{args.el}.txt", FLH, delimiter=' ', fmt='%7.10f %7.10f %7.3f')
        elif transf == "flh2XYZ":
            XYZ = self.flh2XYZ(np.deg2rad(dane[:,0], np.deg2rad(dane[:,1]), dane[:,2]))
            np.savetxt(f"plik_wynikowy_{transformacja}_{args.el}.txt", XYZ, delimiter =' ', fmt ='%7.3f %7.3f %7.3f' )
        elif transf == "XYZ2NEU":
            NEU = self.XYZ2NEU(dane[])"""

                               
"""if __name__ == "__main__":
    # Tworzenie parsera argumentów
    
    parser = argparse.ArgumentParser(description="Przykładowy program z użyciem argparse")
    parser.add_argument("--model", type=str, default="WGS84", choices=["WGS84", "GRS80", "Krasowski"],
                        help="Model elipsoidy (domyślnie: WGS84)")
    parser.add_argument("--Plik", type=str, help="Wartość X")
    parser.add_argument("--Elipsoida", type=str, help="Wartość Y")
    parser.add_argument("--Transformacja", type=str, help="Wartość Z")
    parser.add_argument("--output", type=str, default="dec_degree", choices=["dec_degree", "dms"],
                        help="Format wyników (domyślnie: dec_degree)")
    args = parser.parse_args()
    elipsoidy = {'WGS84':[6378137.000, 0.00669438002290], 'GRS80':[6378137.000, 0.00669438002290], 'KRASOWSKI':[6378245.000, 0.00669342162296]}
    transformacje = {'XYZ2BLH': 'XYZ2BLH','BLH2XYZ': 'BLH2XYZ','PL2000':'PL2000','PL1992':'PL1992', 'XYZ2NEUP':'XYZ2NEUP'}
    
    wybor = input(str("Jezeli chcesz wykonac kolejna transformacje wpisz TAK jesli chcesz zakonczyc ENTER: ")).upper()
    args.Elipsoida = None
    args.Plik= None
    args.Transformacja= None"""