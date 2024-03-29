
from math import *
import random 
import numpy as np
#from matplotlib import pyplot as plt
import time

# nimmt ein paar Vektoren und macht sie orthogonal (rechtwinklig) und normal (länge 1)
def orthonormalisierung(ebene, dim=4):
    M = np.array(ebene, dtype=np.float64).reshape(len(ebene),dim)    # macht aus der liste der 2 vekotren nen array (Datentyp für Matrixen / Tensoren/ Vektoren)
    M[0] = M[0] / sqrt(np.dot(M[0],M[0])) # Vektor durch Länge => länge 1

    # Gram- Schmidt verfahren zum orthogonalisieren    
    for k in range(1,len(ebene)): 
        summ =0
        for j in range(0,k):
            summ += np.dot(M[j].T,M[k])* M[j] 

        M[k] = M[k] - summ
        M[k] = M[k] / sqrt(np.dot(M[k],M[k]))
     
    return M.reshape(len(ebene),dim) # Ausgabe als (Zeilen-)Vektor 

# nimmt 2 Vektoren (aus R4) und gibt eine dazu orthogonale (rechtwinklige) Ebene aus
def orthogonal_ebene(ebene, dim = 4):
    v1 = np.array(ebene[0],dtype=float).reshape(dim,1) # (Spalten-)Vektoren aus Liste
    
    V = [ebene[0], ebene[1]]
    E= np.eye(dim) # [(1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1)]
    neu = 0
    for i in range(len(v1)):
        if (-10**(-4) > v1[i]) or (v1[i] > 10**(-4)) : #Falls != 0, aber halt mit rundungsfehler
            neu = i 
            """for j in range(4):
                if j !=i:
                    V += [E[j]]"""
            #print(V)
            break # ich brauche V gar nicht

    # Basiswechselmatrix (multiplikation mit Vektor gibt an, wie oft jeder Basisvektor für den Vektor gebraucht wird) 
    if dim == 4:
        v2 = np.array(ebene[1],dtype=float).reshape(dim,1)
        B = np.eye(4) #Einheitsmatrix -> wird Basiswechselmatrix
        for i in range(4):
            if neu == i:
                B[i,i] = np.array([1/v1[i]])[0,0] 
            else:
                B[i,neu] = -v1[i] 

        v2 = np.dot(B,v2) # np.dot ist Matrixmultiplikation -> macht den Basiswechsel
        #neu2 = 0
        for i in range(0,4):
            if ((-10**(-4) > v2[i]) or (v2[i] > 10**(-4))) and i != neu : # um zu gucken, welcher Vektor der Einheitsmatrix linear unabhängig von den  der Ebene ist
                #neu2 = i
                for j in range(4):
                    if j != i and j != neu:
                        V += [E[j]] # Macht eine Basis der R4, die die Ebenenvektoren enthällt
                break
    elif dim == 3:
        for j in range(dim):
            if all(np.array([j,j,j]) != neu):
                V += [E[j].tolist()]

    V = np.array(V).reshape(dim,dim)
    V = orthonormalisierung(V,dim)

    for i in range(len(V)):
        for i in range(len(V)):
            if i != j and not np.allclose(V[i]@V[j], 0):
                print("FEEEHHHLLEEEER in der Ebene")
            if i==j and not np.allclose(V[i]@V[j], 1):
                print("FEEEHLER ebene nicht normal")

    E = orthonormalisierung(ebene)
    for i in range(2):
        if any( E[i] != V[i]):
            print("FEEEHLER ebene verschoben")            

    return  orthonormalisierung(V,dim)[2:dim,:] #gibt die neuen Vektoren aus, die jetzt auch orthogonal zu dem rest sind

# Ebene in der gedreht wird
def rotation(winkel,ebene,ecken):
        #print(winkel, ebene, ecken)
        # Aus Ebenen-Vektoren Orthonormale Vektoren machen
        vektoren = orthonormalisierung(ebene)

        if np.dot(vektoren[0], vektoren[1].T) !=0  :
            print("FEHLER Ebene falsch bei rotation" , 1)

        if not np.allclose(vektoren[0]@ vektoren[0].T , 1) or not np.allclose(np.dot(vektoren[1], vektoren[1].T), 1):   
            print("FEHLER Ebene falsch bei rotation" , 2)   
            vektoren[0] = vektoren[0] /sqrt(np.dot(vektoren[0], vektoren[0].T))
            vektoren[1] = vektoren[1] /sqrt(np.dot(vektoren[1], vektoren[1].T))

        if not np.allclose(np.dot(vektoren[0], vektoren[0].T), 1): 
            print("dir ist nicht mehr zu helfen!!! 1")    
            print(vektoren[0])

        if not np.allclose(np.dot(vektoren[1], vektoren[1].T), 1):    
            print("dir ist nicht mehr zu helfen!!! 2") 
            print(vektoren[1])

        g1 = vektoren[0].reshape(4,1)
        g2 = vektoren[1].reshape(4,1)

        # Aus zwei Orthonormalen Vektoren die Rotationsmatrix erstellen
        V = np.dot(g1,g1.T) + np.dot(g2,g2.T)
        W = np.dot(g1,g2.T) - np.dot(g2,g1.T)
        I = np.eye(4)
        R = I + (cos(winkel)-1)*V + sin(winkel)*W #Rotationsmatrix

        #Basiswechselmatrix
        #B = np.array(vektoren.tolist() + orthogonal_ebene(vektoren).tolist()).reshape(4,4)
        
        #Rot = np.array([[cos(winkel),-sin(winkel),0,0],[sin(winkel),cos(winkel), 0,0],[0,0,1,0],[0,0,0,1]])

        #R = B * Rot * B.T
        #print(ecken)
        ecken = np.array(ecken).reshape(len(ecken),4)
        # (np.exp( winkel * (np.dot(g1.T, ecken.T) * g2 - np.dot(g2.T, ecken.T)* g1) )).T # np.dot(R,ecken.T).T
                    
        return np.dot(R,ecken.T).T  # Ecken mal Rot. Matrix -> neue Ecken ausgegeben

def spiegelung(normal, punkte):
    punkte = np.array(punkte).reshape(len(punkte),len(punkte[0])).T
    normal = np.array(normal).reshape(len(normal),1)
    H = np.eye(len(normal)) - 2 / np.dot(normal.T,normal) * np.dot(normal, normal.T)
    return np.dot(H,punkte).T

def drei_zu_vier(punkte):
    punkte = np.array(punkte).tolist()
    for i in range(len(punkte)):
        punkte[i] = list(punkte[i]) + [0]
    return punkte    

def vier_zu_drei(punkte):
    for i in range(len(punkte)):
        punkte[i] = list(punkte[i])[0:len(punkte[i]-1)]
    return punkte    

def bereich_test(punkte, normalenvektoren):
    Bereiche = []
    D_vergl = []
    n = len(normalenvektoren)

    for i in range(2**len(normalenvektoren)):
        Bereiche += [0]
        a = i
        D_vergl += [[]]
        for j in range(n):
            if a% 2 == 1:
                D_vergl[i] += [1]
                a = a - 1
            else: 
                D_vergl[i] += [0]
            a = a/2
        D_vergl[i] = list(reversed(D_vergl[i]))

    d = [-1 for i in range(n)]

    for p in punkte:
        for i in range(n):
            if not np.dot(p,np.array(Normal[i]).reshape(4,1)) < 0:
                d[i] = 1 
            else: d[i] = 0    
        
        for b in range(len(D_vergl)):
            if d == D_vergl[b]:
                Bereiche[b] += 1    
    return Bereiche

def winkel(vektor1, vektor2):
    vektor1 = vektor1 / sqrt(vektor1@vektor1)
    vektor2 = vektor2 / sqrt(vektor2@vektor2)
    return acos(np.dot(vektor1, vektor2.T))

class raum():
    
    def __init__(self, normalenvektor) -> None:
        self.normal = normalenvektor
        pass

    # Ebene, die fest bleibt!!!    
    def rotation(self,winkel,ebene):
        E = orthonormalisierung(ebene)
        E = orthogonal_ebene(E)
        
        epsilon = rotation(winkel,E,self.normal).reshape(1,4)
        return epsilon

    pass

class koerper():

    def __init__(self, ecken) -> None:
        self.ecken = ecken
        pass
    # Ebene, die fest beleibt!!!    
    def rotation(self, winkel, ebene):
        E = orthogonal_ebene(orthonormalisierung(ebene))
        return rotation(winkel, E, self.ecken)
    
    pass

körper_drehung = False

if körper_drehung:
    # Funktioniert erstmal nur für Eckenfigur <3,3>

    #3 Ecken: die erste, die "mittlere", die an der alles aufeinander trifft: NICHT MIT ANGEBEN!!!
    a = - (sqrt(5)+1)/4
    c = (1 - (sqrt(5) + 1)/2)/2
    b = 1/2

    dreid_korper = [(1,0,0),(0,1,0),(0,0,1)]  #[(1, 0, 0), (0.5, sqrt(3)/2, 0) , (0.5, sqrt(3)/6, sqrt(2/3))] #  [(a,b,c), (b,c,a), (c,a,b)]  # [(1,0,0),(0,1,0),(0,0,1)] 

    #fester_korper = koerper(drei_zu_vier(dreid_korper))

    alle_korper = []
    for i in range(3):
        normal = np.cross(dreid_korper[(i+2)%3],dreid_korper[(i+1)%3])
        gespiegelt = spiegelung(normal, [dreid_korper[(i)%3]] )
        langer = drei_zu_vier(gespiegelt)
        alle_korper.append(koerper([langer]))

        
    for i in alle_korper:
        print("-------------" , i.ecken)

    fester_korper = koerper(drei_zu_vier(dreid_korper))
    print(fester_korper.ecken)

    Winkel = np.linspace( 0,pi,1001)
    distanz = [100 for i in range(3) ]
    wirkliche_distanz = [100 for i in range(3)]
    wirklicher_winkel = [10 for i in range(3)]
    kl_winkel = [10 for i in range(3) ]


    for w in Winkel:
        #print()
        gedrehter_korper = []
        for i in range(2):
            gedrehter_korper.append(alle_korper[i].rotation(w,[fester_korper.ecken[(i+2)%3],fester_korper.ecken[(i+1)%3]]))
        gedrehter_korper.append(alle_korper[2].rotation(-w,[fester_korper.ecken[(2+2)%3],fester_korper.ecken[(2+1)%3]]))
            
        """for i in range(3):
            print("gedereht ", i , gedrehter_korper[i])  """  
        
        
        if np.allclose(gedrehter_korper[0], gedrehter_korper[1]):
            print (w , "hier!!!!")


        for i in range(3):
            #print(gedrehter_korper[i])

            b = np.max(gedrehter_korper[i] - gedrehter_korper[(i+1)%3])
            if b < distanz[i]:
                distanz[i] = b
                kl_winkel[i] = w
            
            a = (gedrehter_korper[i] - gedrehter_korper[(i+1)%3])@(gedrehter_korper[i] - gedrehter_korper[(i+1)%3]).T
            if  a < wirkliche_distanz[i]:
                wirkliche_distanz[i] = a
                wirklicher_winkel[i] = w

    print(distanz)
    print(kl_winkel)
    print(wirkliche_distanz , wirklicher_winkel)



### --- Variablen --- ###
raumdrehung =True

anzahl_punkte = 1_000_000
anzahl_winkel = 100
Winkel = [1.82347889] #np.linspace(0,2*pi,anzahl_winkel)
r = 500 # Radius
dim = 4 # Dimension
a = 3 #Anzahl der drehebenen
Daten = [] #da werden die Daten gespeichert
Winkel_Daten = []

Raum_anfang = raum([(0,0,0,1)])
print(Raum_anfang.normal)
Ebenen_Matrix = np.array([[(1, 0, 0, 0), (0.5, sqrt(3)/2, 0, 0) , (0.5, sqrt(3)/6, sqrt(2/3), 0)]]).reshape(a,dim)

randomm = True

### --- --- ###


punkte_zaehler = 0
punkte = []
anfangs_zeit = time.time()

if raumdrehung:

    if not randomm:
        zeit = time.time()
        for i in  range(-r,r):
            for j in range(-r,r):
                for k in range(-r,r):
                    for l in range(-r,r):
                        p = np.array([i/r,j/r,k/r,l/r]).reshape(1,4)
                        if np.dot(p,p.T) <= 1:
                            punkte += [p]
        print(time.time()-zeit)             
        print(len(punkte))


    if randomm:
        while punkte_zaehler < anzahl_punkte:
            p = np.array([[random.randint(-r,r)/r for i in range(4)]]).reshape(1,4)
            
            if np.dot(p,p.T) <= 1: #wenn in Späre
                punkte_zaehler += 1
                punkte += [p]


    for w in Winkel:

        zeit = time.time()
        print("still running" , w)
        Normal = [Raum_anfang.normal]
        Winkel_Neu = []

        for i in range(a):
            if i != 2:
                neu = (Raum_anfang.rotation(w, [Ebenen_Matrix[i],Ebenen_Matrix[(i+1)%3]])).reshape(1,4).tolist()

            else: neu = (Raum_anfang.rotation(-w, [Ebenen_Matrix[i],Ebenen_Matrix[(i+1)%3]])).reshape(1,4).tolist()
            neuer_winkel = winkel(np.array(Raum_anfang.normal) , np.array(neu))

            Normal.append(neu)
            Winkel_Neu.append(neuer_winkel)

        Daten.append(bereich_test(punkte, Normal))
        print(time.time()-zeit)


    #with open('el_20240106_07.txt', "a") as speicher:
     #   print(punkte_zaehler, "; \n", np.array(Daten).reshape(len(Winkel),16).tolist(), " \n ;", Winkel,  "\n \n", file=speicher)

    print(time.time()- anfangs_zeit)    
    print("fretig")
    print(Daten)