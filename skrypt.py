# -*- coding: utf-8 -*-
"""
Created on Sun Apr 23 21:11:54 2023

@author: USER
"""
import numpy as np 
from math import *
from argparse import ArgumentParser


class Transformacje:
    def __init__(self, elip):
        print('class', elip)
        self.a = elip[0]
        self.ep2 = elip[1]
    
    """
    Definicja Np
    """
    
    def Np(self,f):
        N=self.a/np.sqrt(1-self.ep2*np.sin(f)**2)
        return N
    
    """
    Tranformacja współrzędnych geocentrycznych XYZ na współrzędne elipsoidalne fi, lambda, h
    """
    
    def XYZ2FLH(self, X, Y, Z):
        transformed = []
        
        """Zastosowano algorytm Hirvonena, transformujący współrzędne prostokątne na współrzędne elipsoidalne. W procesie iteracyjnym, uzyskujemy dokładne wyniki"""
        for X, Y, Z in zip(X, Y, Z):
            r   = np.sqrt(X**2 + Y**2)           # promień
            Fi = atan(Z / (r * (1 - self.ep2)))    # pierwsze przybliilizenie
            fi = 0
            while abs(Fi - fi) > 0.000001/206265:    
                Fi = fi
                N = self.Np(Fi)
                h = r / cos(Fi) - N
                lat = atan((Z/r) * (((1 - self.ep2 * N/(N + h))**(-1))))
                break
            lam = np.arctan(Y/X)
            N = self.a / sqrt(1 - self.ep2 * (sin(Fi))**2);
            h = r / cos(Fi) - N  
            transformed.append([np.rad2deg(Fi), np.rad2deg(lam), h])
        return transformed
        

        """
        Transformacja współrzędnych elipsoidalnych fi, lambda, h na współrzędne XYZ
        """
        def FLH2XYZ(self,Fi,lam,h):
            transformed = []
            while True:
                N=self.Np(Fi)
                Xk=(N+h)*np.cos(Fi)*np.cos(lam)
                Yk=(N+h)*np.cos(Fi)*np.sin(lam)
                Zk=(N*(1-self.ep2)+h)*np.sin(Fi)
                Xp = 0
                if abs(Xp-Xk)<(0.000001/206265):
                    break
                transformed.append([Xk, Yk, Zk])
            return transformed
        
        """
        Tranformacja współrzędnych geocentryczny do współrzędnych topocentrycznych
        """
        def Rneu(f,l):
            R = np.array([[-np.sin(f)*np.cos(l), -np.sin(l), np.cos(f)*np.cos(l)],
                          [-np.sin(f)*np.sin(l), np.cos(l), np.cos(f)*np.sin(l)],
                          [np.cos(f), 0., np.sin(f)]])
            return(R)
        
    
        def XYZ2NEU(self, X, Y, Z, X0, Y0, Z0):
            transformed = []
            r   = np.sqrt(X**2 + Y**2)           # promień
            Fi = atan(Z / (r * (1 - self.ep2)))    # pierwsze przybliilizenie
            fi = 0
            while abs(Fi - fi) > 0.000001/206265:    
                Fi = fi
                N = self.Np(f)
                h = r / np.cos(Fi) - N
                lat = atan((Z/r) * (((1 - self.ep2 * N/(N + h))**(-1))))
                break
            lam = np.arctan(Y/X)
            N = self.a / sqrt(1 - self.ep2 * (sin(fi))**2);
            h = r / cos(fi) - N  
            
            R_neu = self.Rneu(f,l)
            for X, Y, Z in zip(X, Y, Z):
                X_sr = [X-X0, Y-Y0, Z-Z0] 
                X_rneu = R_neu.T@X_sr
                transformed.append(X_rneu.T)
                
            return transformed

        """
        Tranformacja współrzędnych fi, lambda do układu 2000
        """
        
        def sigma(self, f):
            A0 = 1 - self.ep2/4 - 3 * self.ep2**2/64 - 5 * self.ep2**3/256
            A2 = (3/8) * (self.ep2 + self.ep2**2/4 + 15 * self.ep2**3/128)
            A4 = (15/256) * (self.ep2**2 + 3 * self.ep2**3/4)
            A6 = 35 * self.ep2**3/3072
            sig = self.a * (A0 * f - A2 * np.sin(2 * f) + A4 * np.sin(4 * f) - A6 * np.sin(6 * f))
            return(sig)
        
        def GK2000(self, f, l, m=0.999923):
            transformed = []
            for f, l in zip(f,l):
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
        
            b2 = (slef.a**2) * (1-self.ep2)   #krotsza polos
            e2p = ( self.a**2 - b2 ) / b2   #drugi mimosrod elipsy
            dl = l - l0
            t = np.tan(f)
            ni = np.sqrt(e2p * (np.cos(f))**2)
            N = self.Np(f)
            sigma = self.sigma(f)
            XGK20 = sigma + ((dl**2)/2)*N*np.sin(f)*np.cos(f) * ( 1 + ((dl**2)/12)*(np.cos(f))**2 * ( 5 - (t**2)+9*(ni**2) + 4*(ni**4)     )  + ((dl**4)/360)*(np.cos(f)**4) * (61-58*(t**2)+(t**4) + 270*(ni**2) - 330*(ni**2)*(t**2))  )
            YGK20 = (dl*N* np.cos(f)) * (1+(((dl)**2/6)*(np.cos(f))**2) *(1-(t**2)+(ni**2))+((dl**4)/120)*(np.cos(f)**4)*(5-18*(t**2)+(t**4)+14*(ni**2)-58*(ni**2)*(t**2)) )
            X2000 = XGK20 * m 
            Y2000 = YGK20 * m + strefa*1000000 + 500000
            transformed.append([X2000, Y2000])
            
        return transformed
    
    """
    Tranformacja współrzędnych fi, lambda do układu 1992
    """
    
    def GK1992(self, f, l, m = 0.9993):
        transformed = []
        lam0 = (19*np.pi)/180
        for f, l in zip(f,l):
            b2 = (self.a**2) * (1-self.ep2)   #krotsza polos
            e2p = ( self.a**2 - b2 ) / b2   #drugi mimosrod elipsy
            dlam = l - lam0
            t = np.tan(f)
            ni = np.sqrt(e2p * (np.cos(f))**2)
            N = self.Np(f)

            sigma = self.sigma(f)

            xgk = sigma + ((dlam**2)/2)*N*np.sin(f)*np.cos(f) * ( 1+ ((dlam**2)/12)*(np.cos(f))**2 * ( 5 - (t**2)+9*(ni**2) + 4*(ni**4)     )  + ((dlam**4)/360)*(np.cos(f)**4) * (61-58*(t**2)+(t**4) + 270*(ni**2) - 330*(ni**2)*(t**2))  )
            ygk = (dlam*N* np.cos(f)) * (1+(((dlam)**2/6)*(np.cos(f))**2) *(1-(t**2)+(ni**2))+((dlam**4)/120)*(np.cos(f)**4)*(5-18*(t**2)+(t**4)+14*(ni**2)-58*(ni**2)*(t**2)) )
            
            x92 = xgk*m - 5300000
            y92 = ygk*m + 500000
            
            transformed.append([x92, y92])
        
        return transformed 

    def plik(self, pliczek, transf):
        dane = np.genfromtxt(pliczek,delimiter = ' ')
        print('tu zapisuje')
        if transf == 'XYZ2FLH':
            print('tu zapisuje')
            transformed  = self.XYZ2FLH(dane[:,0], dane[:,1], dane[:,2])
            np.savetxt(f"dane_{transf}_{args.elip}.txt", transformed, delimiter=' ', fmt='%0.10f %0.10f %0.3f')
        elif transf == 'FLH2XYZ':
            transformed  = self.FLH2XYZ(np.deg2rad((dane[:,0])), np.deg2rad(dane[:,1]), dane[:,2])
            np.savetxt(f"plik_wynikowy_{transf}_{args.elip}.txt", transformed, delimiter =' ', fmt ='%0.3f %0.3f %0.3f' )
        elif transf == 'XYZ2NEU':
            transformed  = self.XYZ2NEU(dane[1:,0], dane[1:,1], dane[1:,2], dane[0,0], dane[0,1], dane[0,2])
            np.savetxt(f"plik_wynikowy_{transf}_{args.elip}.txt", transformed, delimiter =' ', fmt ='%0.3f %0.3f %0.3f' )
        elif transf == 'GK2000':
            transformed  = self.GK2000(np.deg2rad(dane[:,0]), np.deg2rad(dane[:,1]))
            np.savetxt(f"plik_wynikowy_{transf}_{args.elip}.txt", transformed, delimiter=' ', fmt='%0.3f %0.3f')
        elif transf == 'GK1992':
            transformed  = self.GK1992(np.deg2rad(dane[:,0]), np.deg2rad(dane[:,1]))
            np.savetxt(f"plik_wynikowy_{transf}_{args.elip}.txt", transformed, delimiter=' ', fmt='%0.3f %0.3f')

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('-pliczek', type=str, help='Wpisz sciezke do pliku z danymi wejsciowymi')
    parser.add_argument('-elip', type=str, help='Wybierz elipsoide sposrod dostepnych: WRS84, GRS80, KRASOWSKI')
    parser.add_argument('-transf', type=str, help='Wybierz transformacje, z ktorej chcesz skorzystac, sposrod dostepnych: XYZ2flh, flh2XYZ, saz2neu, GK2000, GK1992, XYZ2NEU')
    args = parser.parse_args()
    elip = {'WGS84': [6378137.000, 0.00669438002290], 'GRS80': [6378137.000, 0.00669438002290], 'KRASOWSKI': [6378245.000, 0.00669342162296]}
    transf = {'XYZ2FLH': 'XYZ2FLH', 'FLH2XYZ': 'FLH2XYZ','XYZ2NEU': 'XYZ2NEU', 'GK2000': 'GK2000', 'GK1992': 'GK1992'}
    # print(args.elip, elip[args.elip])
    # print(elip['WGS84'][1])
    # print(args.pliczek)
    # wsp = Transformacje(elip[args.elip.upper()])
    # wczyt = wsp.plik(args.pliczek, transf[args.transf.upper()])
    
    try:
        wsp = Transformacje(elip[args.elip])
        wczyt = wsp.plik(args.pliczek, transf[args.transf.upper()])
        print("Utworzono plik ze wspolrzednymi.")
    except AttributeError as e:
        print("Error:", e)    
    except FileNotFoundError:
        print("Nie znaleziono podanego pliku")
    except KeyError:
        print("Niepoprawna nazwa Elipsoidy lub Transformacji")
    except IndexError:
        print("Dane w podanym pliku sa w nieodpowiednim formacie")
    except ValueError:
        print("Dane w podanym pliku sa w nieodpowiednim formacie")
    finally:
        print("Mamy nadzieję, że nasz program był dla Ciebie użyteczny")