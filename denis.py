def objekt_falt(sekund):
    tyngdekraft = 9.81
    distanse = tyngdekraft * sekund 
    return distanse

antall_sekund = int(input("Hvor lenge har objekt 1 falt?\n"))
distanse1 = objekt_falt(antall_sekund) # Kaller på funksjonen

antall_sekund = int(input("Hvor lenge har objekt 2 falt?\n"))
distanse2 = objekt_falt(antall_sekund) # Kaller på funksjonen

antall_sekund = int(input("Hvor lenge har objekt 3 falt?\n"))
distanse3 = objekt_falt(antall_sekund) # Kaller på funksjonen

print("Objektene har falt %, %, og % meter", distanse1, distanse2, distanse3)


def sjekk_dato(dato):
    if (dato >= 1 and dato <= 31):
        return True
    else:
        return False



# Koden som bruker funksjonen  år
input_dato = int(input("Oppgi leveringsdato: "))

if (sjekk_dato(input_dato) == False):
    print("Du skrev inn en ugyldig dato!")
else:
    print("Pakken blir levert %", input_dato)




