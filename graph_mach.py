
from math import *
import random 
import numpy as np
from matplotlib import pyplot as plt #für die graphiken
import time

# x Werte und eine Liste, die zu jeden x-Werte eine Liste der y-Werte enthällt
def graph_off(x, y_werte):
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

    
    #plt.plot(x, Summe , label = "Summe" , color = 'black')
    plt.xlabel('Winkel')
    plt.ylabel('Punkte')
    #plt.ylim(10,-10)
    plt.legend()
    plt.show()
    plt.savefig('el_20231228_02_01.png')
    
Daten =   
Winkel = np.linspace(0,2*pi,100)
graph_off(Winkel, Daten)