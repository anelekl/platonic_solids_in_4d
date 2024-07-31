from math import *
import random
import numpy as np
#from matplotlib import pyplot as plt
import time

#nimmt eine Liste von Zeilenvektoren, gibt eine Orthonormalbasis des Spanns der
#Vektoren durch Gram-Schmidt (als Zeilenvektoren einer Matrix)
def orthonormalisierung(vektoren) -> np.ndarray:
    M = np.array(vektoren, dtype=np.float64)
    anzahl = len(vektoren)
    k = 0
    while k<anzahl:
        M[k] = M[k] - sum([M[j]@M[k] * M[j] for j in range(k)])
        if(np.linalg.norm(M[k]) < 1e-10):
            anzahl -= 1
            for i in range(k,anzahl):
                M[i]=M[i+1]
            continue
        M[k] = M[k] / np.linalg.norm(M[k])
        k += 1
    return M[:anzahl]

#nimmt eine Liste von Vektoren, gibt eine Orthonormalbasis des zum Spann der
#Vektoren orthogonalen Raums (als Zeilenvektoren einer Matrix)
def orthogonal_raum(vektoren, dim:int = 4) -> np.ndarray:
    ortho_mat = orthonormalisierung(vektoren)
    k = ortho_mat.shape[0]
    return orthonormalisierung(np.concatenate((ortho_mat,np.eye(dim))))[k:]

def orthogonal_ebene(vektoren,dim:int=4) -> np.ndarray:
    assert(len(vektoren)+2 == dim)
    return orthogonal_raum(vektoren,dim)

#Gibt Matrix zurück, die eine Rotation um den Winkel in der aus den ersten bei-
#den Basisvektoren gebildeten Ebene beschreibt
def rot_matrix(winkel:float,basis:np.ndarray=np.eye(4)) -> np.ndarray:
    dim = len(basis)
    assert(dim >= 2)
    mat2d = np.array([[np.cos(winkel),np.sin(winkel)],
                    [-np.sin(winkel),np.cos(winkel)]])
    return basis@np.concatenate(
           (np.concatenate((mat2d,np.zeros((2,dim-2))),axis=1),np.eye(dim)[2:])
           )@np.linalg.inv(basis)

def rotation(winkel:float,ebene,ecken) -> np.ndarray:
    assert(len(ebene[0]) == 4 and len(ecken[0]) == 4)
    vektoren = orthonormalisierung(ebene)
    assert(len(vektoren) == 2)
    basis = np.concatenate(vektoren,orthogonal_raum(vektoren))
    R = rot_matrix(winkel,basis)
    return ecken@R.T

def spiegelung(normale, punkte):
    normale/=np.linalg.norm(normale)
    return np.array([punkt-2*normale*normale@punkt for punkt in punkte])

def drei_zu_vier(punkte):
    return np.concatenate((np.array(punkte),np.zeros((len(punkte),1))),axis=1)

def vier_zu_drei(punkte):
    return np.array(punkte)[:,:3]

#gibt die Orientierungen der Punkte zu den Ebenen als Dictionary zurück
def bereich_test(punkte, normalen):
    bereiche = {}
    for punkt in punkte:
        bereich=[0 if punkt@normale<0 else 1 for normale in normalen]
        if bereich in bereiche:
            bereiche[bereich]+=1
        bereiche[bereich]=1
    return bereiche

#berechnet Winkel zwischen 2 vektoren
def winkel(vektor1, vektor2):
    vektor1 = vektor1 / np.linalg.norm(vektor1)
    vektor2 = vektor2 / np.linalg.norm(vektor2)
    return acos(vektor1@vektor2)

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

#print(orthogonal_raum([[1,0,0,3],[5,9,27,15]],4))
#print(orthogonal_ebene([[1,0,0,3],[5,9,27,15]],4))
print(rot_matrix(0))
