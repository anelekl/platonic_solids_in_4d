
from math import *
import random 
import numpy as np
from matplotlib import pyplot as plt
import time

dim = 4

# nimmt ein paar Vektoren und macht sie orthogonal (rechtwinklig) und normal (länge 1)
def orthonormalisierung(ebene, dim = 4):
    M = np.array(ebene, dtype=np.float64).reshape(len(ebene),dim)    # macht aus der liste der 2 vekotren nen array (Datentyp für Matrixen / Tensoren/ Vektoren)
    M[0] = M[0] / sqrt(M[0]@M[0]) # Vektor durch Länge => länge 1

    # Gram- Schmidt verfahren zum orthogonalisieren    
    for k in range(1,len(ebene)): 
        summ =0
        for j in range(0,k):
            summ= np.dot(M[j].T,M[k])* M[j] 

        M[k] = M[k] - summ
        M[k] = M[k] / sqrt(M[k]@M[k])
     
    return M.reshape(len(ebene),dim) # Ausgabe als (Zeilen-)Vektor 

# nimmt 2 Vektoren (aus R4) und gibt eine dazu orthogonale (rechtwinklige) Ebene aus
def orthogonal_ebene(ebene, dim = 4):
    v1 = np.array(ebene[0],dtype=float).reshape(dim,1) # (Spalten-)Vektoren aus Liste
    v2 = np.array(ebene[1],dtype=float).reshape(dim,1)

    V = [ebene[0], ebene[1]]
    E= np.eye(dim) #[(1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1)]
    neu = 0
    for i in range(len(v1)):
        if (-10**(-4) > v1[i]) or (v1[i] > 10**(-4)) : #Falls != 0, aber halt mit rundungsfehler
            neu = i 
            """for j in range(dim):
                if j !=i:
                    V += [E[j].tolist()]"""
            #print(V)
            break # ich brauche V gar nicht

    # Basiswechselmatrix (multiplikation mit Vektor gibt an, wie oft jeder Basisvektor für den Vektor gebraucht wird) 
    
    B = np.eye(dim) #Einheitsmatrix -> wird Basiswechselmatrix
    for i in range(dim):
        if neu == i:
            B[i,i] = np.array([1/v1[i]])[0,0] 
        else:
            B[i,neu] = -v1[i] 

    v2 = np.dot(B,v2) # np.dot ist Matrixmultiplikation -> macht den Basiswechsel
    #neu2 = 0
    for i in range(0,dim):
        if (-10**(-4) > v2[i]) or (v2[i] > 10**(-4)) and i != neu : # um zu gucken, welcher Vektor der Einheitsmatrix linear unabhängig von den  der Ebene ist
            #neu2 = i
            for j in range(dim):
                if j != i and j != neu:
                    V += [E[j].tolist()] # Macht eine Basis der R4, die die Ebenenvektoren enthällt
            break
    #print(V)
    V = np.array(V).reshape(dim,dim)

    return  orthonormalisierung(V)[2:dim,:] #gibt die neuen Vektoren aus, die jetzt auch orthogonal zu dem rest sind

# Ebene in der gedreht wird
def rotation(winkel,ebene,ecken,dim = 4):

        # Aus Ebenen-Vektoren Orthonormale Vektoren machen
        vektoren = orthonormalisierung(ebene)

        g1 = vektoren[0].reshape(dim,1)
        g2 = vektoren[1].reshape(dim,1)

        # Aus zwei Orthonormalen Vektoren die Rotationsmatrix erstellen
        V = np.dot(g1,g1.T) + np.dot(g2,g2.T)
        W = np.dot(g1,g2.T) - np.dot(g2,g1.T)
        I = np.eye(dim)
        R = I + (cos(winkel)-1)*V + sin(winkel)*W #Rotationsmatrix
        
        ecken = np.array(ecken).reshape(len(ecken),dim)
        return np.dot(R,ecken.T)  # Ecken mal Rot. Matrix -> neue Ecken ausgegeben

def graph_off(x, y_werte):
    y_listen = []
    Summe = []
    for i in range(len(y_werte[0])):
        y_listen += [[]]
        for j in y_werte:
            y_listen[i] += [j[i]]
            if i == 0:
                Summe += [sum(j)]
    print(x,y_listen)

    #Das Bild
    plt.figure(figsize=(8,6))

    farben = [ 'dimgray' , 'silver', 'lightcoral', 'firebrick', 'darkred' ,'red', 'orangered' , 'darkorange', 'gold', 'yellow' , 'yellowgreen' , 'limegreen', 'forestgreen' , 'darkgreen' , 'teal' , 'royalblue' , 'navy' , 'indigo']

    for y in range(len(y_listen)):
        plt.plot(x,y_listen[y], label = "Bereich " + str(y+1) , color = farben[y])

    
    #plt.plot(x, Summe , label = "Summe" , color = 'black')
    plt.xlabel('x')
    plt.ylabel('y')
    #plt.ylim(10,-10)
    plt.legend()
    plt.show()

#Räume
class raum():
    
    def __init__(self, normalenvektor, dim = 4) -> None:
        self.normal = normalenvektor
        self.dim = dim
        pass

    # Ebene, die fest bleibt!!!    
    def rotation(self,winkel,ebene):
        E = orthonormalisierung(ebene)
        E = orthogonal_ebene(E)
        
        epsilon = rotation(winkel,E,self.normal).reshape(1,self.dim)
        return epsilon

    pass

n = 4 #Anzahl der Räume, meist 4, <= 4 !!!
winkel = pi/2

#Raum durch Normalenvekotor
Raums = []
for i in range(n):
    Raums += [raum([(0,0,0,1)])]

#Normal = np.array()    
#print("Normal", Normal)

Bereiche = []
D_vergl = [[1,1,1,1],[1,1,1,-1],[1,1,-1,1],[1,1,-1,-1],[1,-1,1,1],[1,-1,1,-1],[1,-1,-1,1],[1,-1,-1,-1],[-1,1,1,1],[-1,1,1,-1],[-1,1,-1,1],[-1,1,-1,-1],[-1,-1,1,1],[-1,-1,1,-1],[-1,-1,-1,1],[-1,-1,-1,-1]]
for i in range(2**len(Raums)):
        Bereiche += [0] 

# ---Hier gehts richtig los--------------------------------------------------------


anfang = time.time()
r = 500 #Radius
d = []
Winkel = np.linspace(0,2*pi,100)
Daten = []

for x in range(len(Winkel)):
    print("still running")
    for b in range(len(Bereiche)):
        Bereiche[b] = 0

    punkte = 0
    Normal = np.eye(dim)[-1].reshape(1,dim).tolist()   #[(0,0,0,1)]
    E = np.eye(dim)
    for i in range(0,len(Raums)-1):
        plus = (Raums[i].rotation(Winkel[x], [E[i],E[(i+1)%3]])).reshape(1,dim).tolist()
        Normal += plus     #  np.array(plus).reshape(i+2,4)
    #print(Normal) 
    #print("vorher", Daten)
    while punkte < 10_000:
        d = [0 for i in range(n)]
        #unnotige_punkte += 1
        p = np.array([[random.randint(-r,r)/r for i in range(dim)]]).reshape(dim,1)
        #print(p)
        if np.dot(p.T,p) <= 1:
            punkte += 1
            #print(punkte)
            for i in range(n):
                if np.dot(p.T,np.array(Normal[i]).reshape(dim,1)) < 0:
                    d[i] = -1
                else: d[i] = 1 
            for b in range(len(D_vergl)):
                if d == D_vergl[b]:
                    Bereiche[b] += 1
    #print(Bereiche)    
    #print(Daten)
    Daten += Bereiche
    #print(Daten)
print(anfang, time.time())
print(time.time()-anfang)
print("fertig")

with open('el_20231228_01.txt', "a") as speicher:
    print(punkte, ";", Winkel, ";", np.array(Daten).reshape(len(Winkel),16).tolist(), "\n \n", file=speicher)

