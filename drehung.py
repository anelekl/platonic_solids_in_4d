
from math import *
import random 
import numpy as np
#from matplotlib import pyplot as plt
import time

# nimmt ein paar Vektoren und macht sie orthogonal (rechtwinklig) und normal (länge 1)
# Dimension kann angegeben werden, default 4
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
#durch orthogonal_raum eigentlich überflüssig geworden
def orthogonal_ebene(ebene, dim = 4):
    v1 = np.array(ebene[0],dtype=float).reshape(dim,1) # (Spalten-)Vektoren aus Liste
    
    V = [ebene[0], ebene[1]]
    E= np.eye(dim) # Einheitsmatrix [(1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1)]
    neu = 0

    #sucht Einheitsvektor, der v1 ersetzen kann
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
        B = np.eye(dim) #Einheitsmatrix -> wird Basiswechselmatrix

        for i in range(dim):
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

    # fehlersuche
    for i in range(len(V)):
        for i in range(len(V)):
            if i != j and not np.allclose(V[i]@V[j], 0):
                print("FEEEHHHLLEEEER in der Ebene")
            if i==j and not np.allclose(V[i]@V[j], 1):
                print("FEEEHLER ebene nicht normal")

    #fehlersuche
    E = orthonormalisierung(ebene)
    for i in range(2):
        if any( E[i] != V[i]):
            print("FEEEHLER ebene verschoben")            

    return  orthonormalisierung(V,dim)[2:dim,:] #gibt die neuen Vektoren aus, die jetzt auch orthogonal zu dem rest sind

# wie orthogonale Ebene, nur für alle dimensionen verwendbar (dim von objekt und vor ortho objekt definierbar)
# Vektoren, zu denen etwas orthogonal sein soll, dim von raum in dem das ganze ist
def orthogonal_raum(vektoren, dim = 4):
    for i in range(len(vektoren)):
        while len(vektoren[i]) < dim:
            vektoren[i] = drei_zu_vier([vektoren[i]])
        if len(vektoren[i]) > dim:
            print("dimensionsfehler ")
            return

    V = [np.array(vektoren[i],dtype=float).reshape(dim,1) for i in range(len(vektoren))]
    E = np.eye(dim) 
    neu = []
    B = np.eye(dim) #Einheitsmatrix -> wird Basiswechselmatrix

    for j in range(len(V)): 
        v = V[j]
        
        # schafft Basiswechselmatrix
        if j!= 0:
            for i in range(dim):
                if neu[j-1] == i:
                    B[i,i] = np.array([1/V[j-1][i]])[0,0] 
                else:
                    B[i,neu[j-1]] = -V[j-1][i] 

        v = np.dot(B,v) # np.dot ist Matrixmultiplikation -> macht den Basiswechsel
        

        #sucht Einheitsvektor, der LU ist
        for i in range(len(v)):
            if (-10**(-4) > v[i]) or (v[i] > 10**(-4)) and all( np.full((len(neu)), i) != neu): #Falls != 0, aber halt mit rundungsfehler
                neu += [i] 
                break
    #print(neu)
    for j in range(dim):
        print(np.full((len(neu)),j))
        print( all(np.full((len(neu)), j) != neu))
        if not any(np.full((len(neu)), j) == neu):
            V += [np.array(E[j]).reshape(dim,1)] # Macht eine Basis des R^dim, die die Ebenenvektoren enthällt
            #print(j, V)
    #print(V)
    V = np.array(V).reshape(dim,dim)
    return  orthonormalisierung(V,dim)[len(vektoren):dim,:]              



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

# Spiegelung von den Punkten (punkte) entlang einer Ebene, dessen NOrmalenvektor angegeben ist (normal) [glaube ich]
def spiegelung(normal, punkte):
    punkte = np.array(punkte).reshape(len(punkte),len(punkte[0])).T
    normal = np.array(normal).reshape(len(normal),1)
    H = np.eye(len(normal)) - 2 / np.dot(normal.T,normal) * np.dot(normal, normal.T)
    return np.dot(H,punkte).T

# 3-dim zu 4-dim Vektoren durch hinten 0 hinzufügen
def drei_zu_vier(punkte):
    l = len(punkte[0])
    print(punkte)
    punkte = np.array(punkte).tolist()
    for i in range(len(punkte)):
        punkte[i] = list(punkte[i]) + [0]
    return np.array(punkte).reshape(l+1)    

# 4-dim zu 3-dim vektoren durch hinten ein entfernen
def vier_zu_drei(punkte):
    for i in range(len(punkte)):
        punkte[i] = list(punkte[i])[0:len(punkte[i]-1)]
    return punkte    

# testet in welchem Bereich zwischen den Ebenen, die zu den Normalenvektoren gehören, die Punkte liegen und gibt Liste mit anzahl der Pnkten in jedem Bereich an
# wäre super, wenn dass für meherere Dimensionen funktionieren würde, nicht nur 4
def bereich_test(punkte, normalenvektoren):
    Bereiche = []
    D_vergl = [] #wird Liste der Breichsnamen im binärformat, um dann mit den Punkten verglichen werden zu können
    n = len(normalenvektoren) # Anzahl Ebene

    # für jeden Bereich, erstellt bereichsnamen in D_vergl
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

    d = [-1 for i in range(n)] #wird zu den richtungsangaben für die Punkte (links bzw rechts von den Ebenen --> -1 bzw +1)
    
    # bestimmt Ort für jeden punkt
    for p in punkte:
        for i in range(n):
            #Skalarprodukt < 0 ==> auf der einen Seite 
            if not np.dot(p,np.array(Normal[i]).reshape(4,1)) < 0:
                d[i] = 1 
            # Skalarprodukt >= 0 ==> auf der anderen Seite
            else: d[i] = 0    
        
        # ordnet Punkt einem Bereich zu und zählt ihn
        for b in range(len(D_vergl)):
            if d == D_vergl[b]:
                Bereiche[b] += 1    
    return Bereiche

#berechnet Winkel zwischen 2 vektoren
def winkel(vektor1, vektor2):
    vektor1 = vektor1 / sqrt(vektor1@vektor1)
    vektor2 = vektor2 / sqrt(vektor2@vektor2)
    return acos(np.dot(vektor1, vektor2.T))

# Raum ist durch Normalenvektor bestimmt. Also immer 1 dim weniger als der Raum, in dem er definiert ist
# also in dem Fall immer 3d räume im 4d raum. wäre super, wenn 4d zu allgmeiner Dimension n geändert werden könnte
class raum():
    
    #def durch Vekotr der im rechten Winkel dadrauf steht (normalenvektor)
    def __init__(self, normalenvektor) -> None:
        self.normal = normalenvektor
        pass

    # Ebene, die fest bleibt (bei 4d rotation bleibt immer eine 2d Ebene unverändert und in einer dazu rechtwinkligen wird gedreht)  
    # Ebene durch 2 vektoren angeben, Winkel um den gedreht wird (winkel)
    def rotation(self,winkel,ebene):
        E = orthonormalisierung(ebene) 
        E = orthogonal_raum(E,4)
        
        epsilon = rotation(winkel,E,self.normal).reshape(1,4)
        return epsilon

    pass

#Raum wird eigentlich nicht gebraucht, bzw hat bis jetzt noch die gleichen Eigenschaften, wie Körper, aber Körper wird durch die Ecken des Körpers definiert (ecken)
# die liegen offensichtlich in einem Raum und die Rotation ist genau wie bei Räumen
class koerper():

    def __init__(self, ecken) -> None:
        self.ecken = ecken
        pass
    
    # Ebene, die fest beleibt    
    def rotation(self, winkel, ebene):
        E = orthogonal_raum(orthonormalisierung(ebene),4)
        return rotation(winkel, E, self.ecken)
    
    pass

körper_drehung = False ## körper_drehung macht das Winkel finden

if körper_drehung:
    # Funktioniert erstmal nur für Eckenfigur <3,3>

    #3 Ecken: die erste, die "mittlere", die an der alles aufeinander trifft: NICHT MIT ANGEBEN -- weil die in 0 gesetzt wurde
    # ich glaube, die sind alle unnötig, weil sie nicht mehr gebraucht werden
    a = - (sqrt(5)+1)/4
    c = (1 - (sqrt(5) + 1)/2)/2
    b = 1/2

    # hier stehen ausgeklammert die ganzen Körper, die ich so bruache. wäre super, wenn die in einer Variablen (dictionary) wären und mann sich aussuchen kann, was man hier dann abruft
    # Körper, der nicht gedreht wird
    dreid_korper = [(1,0,0),(0,1,0),(0,0,1)]  #[(1, 0, 0), (0.5, sqrt(3)/2, 0) , (0.5, sqrt(3)/6, sqrt(2/3))] #  [(a,b,c), (b,c,a), (c,a,b)]  # [(1,0,0),(0,1,0),(0,0,1)] 

    #fester_korper = koerper(drei_zu_vier(dreid_korper))

    #setzt die Körper aneinander, wie das bei ner Eckenfigur von <3,3> zu erwarten ist
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

    Winkel = np.linspace( 0,pi,1001) # winkel, die ausprobiert werden, wäre super, wenn man das durch Variablen irgenwo eingeben könnte
    distanz = [100 for i in range(3) ] # variable, in die dann der Abstand zwischen den Ecken geschrieben wird
    wirkliche_distanz = [100 for i in range(3)] #irgendwas, weil dinge nicht funktioniert haben 
    wirklicher_winkel = [10 for i in range(3)] 
    kl_winkel = [10 for i in range(3) ]


    for w in Winkel:
        #print()
        gedrehter_korper = []
        # dreht die 3 körper und fühgt sie gedrehter_korper hinzu
        for i in range(2):
            gedrehter_korper.append(alle_korper[i].rotation(w,[fester_korper.ecken[(i+2)%3],fester_korper.ecken[(i+1)%3]]))
        gedrehter_korper.append(alle_korper[2].rotation(-w,[fester_korper.ecken[(2+2)%3],fester_korper.ecken[(2+1)%3]]))
            
        """for i in range(3):
            print("gedereht ", i , gedrehter_korper[i])  """  
        
        # spruckt richtigen Winkel aus, wenn er getrofen wird
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
raumdrehung =False

anzahl_punkte = 1_000_000 
anzahl_winkel = 100 
Winkel = [1.82347889] #np.linspace(0,2*pi,anzahl_winkel) #Winkel um den/die gedreht wird
r = 500 # Radius für regelmäßige Anordnung
dim = 4 # Dimension
a = 3 #Anzahl der drehebenen
Daten = [] #da werden die Daten gespeichert
Winkel_Daten = []

Raum_anfang = raum([(0,0,0,1)])
#print(Raum_anfang.normal)
Ebenen_Matrix = np.array([[(1, 0, 0, 0), (0.5, sqrt(3)/2, 0, 0) , (0.5, sqrt(3)/6, sqrt(2/3), 0)]]).reshape(a,dim)

randomm = True

### --- --- ###


punkte_zaehler = 0
punkte = []
anfangs_zeit = time.time()


# Macht daten, wie datan_mach_programm.py nur werden hier für jeden Winkel die gleichen Punkte genommen. Das ist deutlich schneller, will ich aber eigentlich nicht
# Also ich hätte gerne die möglichkeit, mit einem boolean hin und her schalten zu können
if raumdrehung:

    # regelmäßige Anordnung der Punkte, ist aber scheiße, sollte trotzdem als Möglichkeit vorhanden sein
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

    # random anordnung der Punkte (bruacht super viel Zeit, dass muss sehr viel schneller werden)
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

        #Drehung der Ebenen mit Fallunterscheidung
        for i in range(a):
            if i != 2:
                neu = (Raum_anfang.rotation(w, [Ebenen_Matrix[i],Ebenen_Matrix[(i+1)%3]])).reshape(1,4).tolist()

            else: neu = (Raum_anfang.rotation(-w, [Ebenen_Matrix[i],Ebenen_Matrix[(i+1)%3]])).reshape(1,4).tolist()
            neuer_winkel = winkel(np.array(Raum_anfang.normal) , np.array(neu))

            Normal.append(neu)
            Winkel_Neu.append(neuer_winkel)

        Daten.append(bereich_test(punkte, Normal)) #testet in welchem Bereich die Punkte sind und gibt eine Liste mit Anzahl für alle bereiche aus (siehe oben)
        print(time.time()-zeit) #misst die Zeit


    #with open('el_20240106_07.txt', "a") as speicher:
     #   print(punkte_zaehler, "; \n", np.array(Daten).reshape(len(Winkel),16).tolist(), " \n ;", Winkel,  "\n \n", file=speicher)

    print(time.time()- anfangs_zeit)    
    print("fretig")
    print(Daten) #das kann gut in ein Dokument geschrieben werden 
    # Bezeichungsformat "name von Programm, dass die Daten erstellt hat"_"Datum, sodass es eindeutig ist, also mit Minuten und sekunden oder so"
    #alternativ Datum nur bis tag und ne automatische durchnummerierung

print(orthogonal_raum([[1,0,0,3],[5,9,27,15]],4))
print(orthogonal_ebene([[1,0,0,3],[5,9,27,15]],4))