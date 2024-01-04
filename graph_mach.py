
from math import *
import random 
import numpy as np
from matplotlib import pyplot as plt #für die graphiken
import time

# x Werte und eine Liste, die zu jeden x-Werte eine Liste der y-Werte enthällt
def graph_off(x, y_werte, vergleich = False):
    y_listen = []
    Summe = []
    for i in range(len(y_werte[0])):
        y_listen += [[]]
        for j in y_werte:
            y_listen[i] += [j[i]]
            if i == 0:
                Summe += [sum(j)]
    #print(x,y_listen)

    #Das Bild
    plt.figure(figsize=(12,8))

    farben = [ 'dimgray' , 'silver', 'lightcoral', 'firebrick', 'darkred' ,'red', 'orangered' , 'darkorange', 'gold', 'yellow' , 'yellowgreen' , 'limegreen', 'forestgreen' , 'darkgreen' , 'teal' , 'royalblue' , 'navy' , 'indigo']

    for y in range(len(y_listen)):
        plt.plot(x,y_listen[y], label = "Bereich " + str(y+1) , color = farben[y])

    if vergleich:

        berg_summe = []
        x = np.linspace(0,2*pi,200)
        for i in range(len(y_listen[3])):
            berg_summe.append(y_listen[2][i] + y_listen[3][i])
        plt.plot(x, berg_summe, color= 'green')

        x = np.linspace(0,pi,100)
        vergleich_aufwärts = []
        vergleich_abwärts = []
        vergleich_hügel = []
        vergleich_ungleicher_hügel = []

        for i in x:
            a = acos(cos(i)**2)
            vergleich_aufwärts.append(100_000/(4*pi)* (2*i - a)  ) # ( a + 2*(pi/2-i))) # = pi + a - 2*i 
            vergleich_abwärts.append(100_000/(4*pi)* (2*pi -2*i -a))
            vergleich_hügel.append(100_000/(4*pi)*a)
            vergleich_ungleicher_hügel.append(100_000/(4*pi)*(1/2*a + 2/(pi**2)*sin(2*i)))

        plt.plot(x, vergleich_aufwärts , label = 'Vergleich', color = 'black') 
        plt.plot(x, vergleich_abwärts , color = 'black')  
        plt.plot(x, vergleich_hügel , color = 'black') 
        plt.plot(x, vergleich_ungleicher_hügel , color = 'red')

           

         
    
    #plt.plot(x, Summe , label = "Summe" , color = 'black')
    plt.xlabel('Winkel')
    plt.ylabel('Punkte')
    #plt.ylim(10,-10)
    plt.legend()
    plt.savefig('el_20240104_04_01.png')
    plt.show()

def graph_vergleich(x, y_werte):
    y_listen = []
    Summe = []
    for i in range(len(y_werte[0])):
        y_listen += [[]]
        for j in y_werte:
            y_listen[i] += [j[i]]
            if i == 0:
                Summe += [sum(j)]
    #print(x,y_listen)

    #Das Bild
    plt.figure(figsize=(12,8))
    farben = [ 'dimgray' , 'silver', 'lightcoral', 'firebrick', 'darkred' ,'red', 'orangered' , 'darkorange', 'gold', 'yellow' , 'yellowgreen' , 'limegreen', 'forestgreen' , 'darkgreen' , 'teal' , 'royalblue' , 'navy' , 'indigo']

    x = np.linspace(0,pi,100)

    fehler = []
    sinus = []
    sin_arccos = []
    mittelwert = []
    random_idee = []
    for i in range(len(y_listen[1])):
        a = acos(cos(x[i])**2)
        fehler.append(y_listen[1][i] - 100_000/(4*pi)*a/2)
        sinus.append(100_000/(4*pi)*((2/pi)**2)*asin(sin(x[i])*cos(x[i])))
        random_idee.append(100_000/(1/2*pi**2)*(1/(4*pi))*(sin(acos(sin(x[i]-pi/2))*2 + pi)))
        
    x_2 = np.linspace(0,pi/2,100)    
    for i in range(len(x_2)):    
        sin_arccos.append(100_000/(1/2*pi**2) *(1/(4*pi))* (sin(acos(x_2[i]*2/pi)*2)))
        mittelwert.append(100_000/(1/2*pi**2) *(1/(4*pi))* (1/3* sin(acos(x_2[i]*2/pi)*2)) + 2/3* 100_000/(4*pi)*((2/pi)**2)*asin(sin(x_2[i])*cos(x_2[i])))


    plt.plot(x, fehler, color = farben[7])
    plt.plot(x, sinus, color = 'black')
    plt.plot(x_2,sin_arccos, color = 'red')
    #plt.plot(x_2, mittelwert, color = 'indigo')
    plt.plot(x, random_idee )

    plt.xlabel('Winkel')
    plt.ylabel('Punkte')
    #plt.ylim(10,-10)
    plt.legend()
    plt.savefig('el_20231227_02_04.png')
    plt.show()

    
Daten = [[132, 0, 126, 0, 124, 0, 143, 0, 0, 111, 0, 120, 0, 134, 0, 110], [94, 10, 118, 9, 114, 5, 130, 9, 12, 115, 22, 111, 12, 119, 12, 108], [117, 22, 99, 23, 98, 31, 101, 26, 26, 98, 25, 100, 15, 104, 26, 89], [83, 38, 82, 35, 66, 54, 89, 40, 48, 98, 34, 76, 45, 90, 44, 78], [73, 56, 66, 45, 70, 42, 80, 58, 51, 80, 57, 69, 66, 59, 54, 74], [60, 72, 62, 76, 64, 67, 46, 52, 79, 58, 74, 53, 52, 60, 61, 64], [52, 81, 45, 81, 51, 71, 46, 78, 91, 41, 80, 35, 81, 43, 78, 46], [33, 93, 40, 81, 37, 94, 39, 84, 84, 41, 97, 26, 97, 26, 90, 38], [21, 93, 17, 106, 19, 121, 20, 117, 105, 19, 103, 11, 96, 25, 111, 16], [9, 130, 4, 120, 10, 113, 3, 120, 122, 7, 122, 5, 114, 7, 107, 7], [10, 127, 2, 120, 5, 123, 15, 118, 118, 4, 111, 5, 114, 3, 121, 4], [25, 100, 23, 98, 11, 98, 26, 109, 106, 26, 97, 21, 111, 20, 108, 21], [33, 80, 43, 94, 33, 89, 25, 98, 90, 39, 87, 37, 91, 27, 101, 33], [47, 63, 61, 81, 42, 88, 44, 88, 96, 48, 66, 44, 77, 36, 75, 44], [50, 54, 37, 59, 69, 79, 76, 69, 59, 52, 68, 49, 69, 72, 65, 73], [66, 56, 57, 50, 83, 42, 73, 57, 61, 85, 50, 85, 36, 72, 63, 64], [93, 39, 85, 34, 73, 39, 91, 38, 37, 70, 45, 81, 44, 88, 43, 100], [104, 29, 103, 23, 120, 20, 91, 18, 26, 84, 32, 88, 23, 89, 35, 115], [108, 7, 111, 13, 111, 12, 131, 12, 10, 112, 9, 118, 15, 113, 17, 101], [147, 0, 106, 0, 123, 0, 111, 0, 0, 128, 0, 136, 0, 135, 0, 114]]  
#Daten = np.array(Daten).reshape(int(len(Daten)/16),16).tolist()

Winkel = np.linspace(0,2*pi,20)
graph_off(Winkel, Daten)