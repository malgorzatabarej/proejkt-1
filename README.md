
## PROJEKT 1
# Transformacje

Dzięki temu programowi użytkownik, możne przetransformować własne współrzędne w zależnosci od potrzeb.

### WYMAGANIA DO OBSŁUGI PROGRAMU:
```
Projekt tworzony był za pomocą programu Python w wersji 3.10.11
Zaimportowana została biblioteka math, numpy oraz argparse.
```

### OPCJE PROGRAMU:
**Dostępne są następujące elipsoidy:**
```
GRS80
WGRS84
KRASOWSKIEGO
```

**Dostępne są następujące transformacje:**
```
XYZ do BLH
BLH do XYZ
XYZ do NEU
BL do PL2000
BL do PL1992
```

### OPIS DZIAŁANIA PROGRAMU:
Program wymaga podania przez użytkownika danych, na podstawie następujących flag:
```
-dane wymaga podania nazwy oraz rozszerzenia pliku, z ktorego użytkownik chce pobrac wspołrzędne
-elip wymaga podania elipsoidy, z której uzytkownik chce skorzystac np. WGS84 (WIELKIMI LITERAMI) 
-transf wymaga podania nazwy transformacji, z której użytkownik chce skorzystać np. XYZ2FLH
```
*obsługiwane przez program transformacje oraz elipsoidy podane są powyżej*

### PRZYKŁADOWE WYWOŁANIE ZA POMOCĄ WIERSZU POLECEŃ
```
1. Otworzenie folderu, w którym znajduje się plik, za pomocą polecenia 'cd'
np. cd C:\Users\asus\Desktop\program

2. Uruchomienie programu:
np. python skrypt.py -dane wsp.txt -elip WGS84 -transf XYZ2FLH
```

Jeżeli wszystko zostało podane prawidłowo utworzy się plik tekstowy z wynikami oraz wyswietli się  następujący komunikat:

```
Utworzono plik ze wspołrzędnymi.
```

### PRZYKŁADOWE PLIKI ZE WSPÓŁRZĘDNYMI
Do działania programu ważny jest poprawny zapis danych w pliku wejściowym, dlatego też na początu tej strony zamieszczone zostały przez nas przykładowe pliki txt dla poszczególnych transformacji. 











