
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
XYZ do FLH
FLH do XYZ
XYZ do NEU
BL do PL2000 (oznakowane w programie jako GK2000)
BL do PL1992 (oznakowane w programie jako GK1992)
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
np. 
cd C:\Users\asus\Desktop\program

2. Uruchomienie programu:
np. 
python skrypt.py -dane wsp.txt -elip WGS84 -transf XYZ2FLH
```

Jeżeli wszystko zostało podane prawidłowo utworzy się plik tekstowy z wynikami oraz wyswietli się  następujący komunikat:

```
Utworzono plik ze wspołrzędnymi.
```

### PRZYKŁADOWE PLIKI ZE WSPÓŁRZĘDNYMI
Do działania programu ważny jest poprawny zapis danych w pliku wejściowym, dlatego też na początu tej strony zamieszczone zostały przez nas przykładowe pliki txt dla poszczególnych transformacji. 

### OMÓWIENIE DANYCH WEJŚCIOWYCH I WYJŚCIOWYCH

**XYZ do FLH**
```
DANE WEJSCIOWE 'DANE_XYZ.txt' I WYJSCIOWE
Plik wejściowy dla tej transformacji zawiera trzy kolumny, kolejno odpowiadające współrzędnym: X, Y, Z podanych w metrach. 
Po zastosowaniu tej transformacji użytkownik otrzymuje plik, zawierający trzy kolumny: kolejno szerokość geograficzną (fi), długość 
geograficzną (lambda) i wysokość punktu(h). Fi i lambda podane są w stopniach, zaś współrzedne h w metrach. 

```
**FLH DO XYZ**
```
DANE WEJŚCIOWE 'DANE_FLH.txt' I WYJSCIOWE
Plik wejściowy zawiera trzy kolumny, odpowiadające kolejno szerokości geograficznej podanej w stopniach, długości geograficznej podanej 
w stopniach oraz wysokości punktu podanej w metrach. W wyniku transformacji użytownik otrzymuje plik z trzema kolumnami i tą samą liczbą wierszy 
co w pliku wejściowym. Kolumny w pliku wyjściowym odpowiadają wartością: współrzędna X, współrzędna Y, współrzedna Z - każda z nich podana jest w metrach. 

```
**XYZ do NEU**
```
DANE WEJSCIOWE 'DANE_NEU.txt' I WYJSCIOWE
Plik wejściowy zawiera trzy kolumny, odpowiadające współrzędnym X, Y, Z podanych w metrach. W wyniku transformacji utworzony zostaje plik zawierający trzy kolumny, 
które odpowiadają współrzędnym X, Y, Z (w metrach) w układzie topocentrycznym. Plik wyjściowy zawiera o jeden wiersz mniej, niż plik wyjściowym. Jest to spowodowane tym, 
że program odczytuje pierwszy wiersz pliku wejściowego jako punt początkowy układu NEU, więc nie bierze go pod uwagę w przeliczaniu współrzędnych.
```

**BL do PL2000/PL1992**
*Szczegóły transformacji BL do PL2000
Przy tej transformacji należy wziąć po uwagę zakresy poszczególnych stref (wartości długości geodezyjnej). 
I STREFA:
<13.5,16.5>
II STREFA:
(16.5,19.5>
III STREFA:
(19.5,22.5>
IV STREFA:
(22.5,25.5>
Jeśli wartość długości geodezyjnej jest mniejsza niż 13,5 stopnia lub większa niż 22,5 stopnia to wynik może być poza zakresem układu PL2000. W takim przypadku wartości X i Y będą nieprawidłowe lub nieprzydatne.*

```
DANE WEJSCIOWE 'DANE_FL_92_20.txt'
Dane wejściowe i wyjściowe dla transformacji GK2000 i GK1992
Plik z danymi zawiera dwie kolumny, pierwsza z nich prezentuje wartości dla szerokości geograficznej podanej w stopniach, zaś druga prezentuje wartości dla 
długości geograficznej podanej w stopniach. Po transformacji użytkownik otrzymuje plik, zawierający dwie kolumny: współrzędne X i Y podane w metrach. 
```
*WSZYSTKIE DANE MAJĄ TYP INTEGER*

### ZNANE BŁĘDY PROGRAMU

```
Transformacja Krasowski do PL2000 daje błędne wyniki, mimo że jest ona dostępna w programie. Nie powinna być ona stosowana przez użytkownika.
```










