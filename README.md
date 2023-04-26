## PROJEKT 1
# Transformacje

Dzięki temu programowi, można transformować współrzędne, w zależnosci od potrzeb użytkownika.

**Program wymaga:**
```
Phyton w jakiejs wersji
Zaimprtowania biblioteki math, numpy oraz argparse
```

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

## OPIS DZIAŁANIA PROGRAMU
Program wymaga podania przez użytkownika danych, przypisanych do następujących flag:
```
-Plik wymaga podania sciezki do pliku, z ktorego chcemy pobrac dane
-Elipsoida wymaga podania konkretnej elipsoidy
-Transformacja wymaga podania nazwy transformacji, z której użytkownik chce skorzystać
```
*obsługiwane przez program transformacje oraz elipsoidy podane są powyżej*

```
Przykładowe polecenie wykonane w CMD:
sciezka -Plik costam -Elipsoida GRS80 -Transformacja XYZ2flh
```

Jeżeli wszystko zostało podane prawidłowo utworzy się plik tekstowy z wynikami oraz wyswietli się  następujący komunikat:

**Utworzono plik ze wspołrzędnymi.**












