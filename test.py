
from math import *
import random 
import numpy as np

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
def basisergenzung(ebene, dim = 4):
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
        #print(V)    
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

    return  orthonormalisierung(V,dim) 

def orthogonal_ebene(ebene, dim = 4):
    return basisergenzung(ebene, dim)[2:dim,:] #gibt die neuen Vektoren aus, die jetzt auch orthogonal zu dem rest sind


def rotation(winkel,ebene,ecken):
        print("Winkel : ", winkel)
        print("Ebene: ", ebene)
        print("Ecken (Anfang) : ", ecken)
        # Aus Ebenen-Vektoren Orthonormale Vektoren machen
        vektoren = orthonormalisierung(ebene)

        print("Ebene ortho. : ", vektoren)

        if not np.allclose(np.dot(vektoren[0], vektoren[1].T), 0)  :
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

        g1 = vektoren[0].reshape((4,1))
        g2 = vektoren[1].reshape((4,1))

        # Aus zwei Orthonormalen Vektoren die Rotationsmatrix erstellen
        V = np.dot(g1,g1.T) + np.dot(g2,g2.T)
        print("V : " , V)
        W = np.dot(g1,g2.T) - np.dot(g2,g1.T)
        print("W : ", W)
        I = np.eye(4)
        print("I : ", I)
        print("cos : ", cos(winkel))
        print("sin : ", sin(winkel))
        print("Ding mal cos : ", np.dot((cos(winkel)-1),V))
        print("Ding mal sin : ", np.dot(sin(winkel),W))

        R = I + np.dot((cos(winkel)-1),V) + np.dot(sin(winkel),W) #Rotationsmatrix
        print("R : ", R)

        #Basiswechselmatrix
        B = basisergenzung(vektoren).reshape(4,4)
        
        for i in range(4):
            for j in range(4):
                if not np.allclose(B[i]@B[j].T , I[i,j]):
                    print("scheiße")     
        
        Rot = np.array([[cos(winkel),-sin(winkel),0,0],[sin(winkel),cos(winkel), 0,0],[0,0,1,0],[0,0,0,1]])

        R_2 = np.dot(B,np.dot(Rot,B.T)) #B * Rot * B.T 

        ecken = np.array(ecken).reshape(len(ecken),4)
        gedreht = np.dot(R,ecken.T).T # (np.exp( winkel * (np.dot(g1.T, ecken.T) * g2 - np.dot(g2.T, ecken.T)* g1) )).T # np.dot(R,ecken.T).T
        gedreht_2 = np.dot(R_2,ecken.T).T

        for i in range(len(ecken)):
            w = np.dot(ecken[i], gedreht[i].T)* 1 / (sqrt(np.dot(gedreht[i],gedreht[i].T))* sqrt(np.dot(ecken[i], ecken[i].T)))
            w_2 = np.dot(ecken[i], gedreht_2[i].T)* 1 / (sqrt(np.dot(gedreht_2[i],gedreht_2[i].T))* sqrt(np.dot(ecken[i], ecken[i].T)))
            
            if w < 1 and w_2 < 1:
                print("Diff : ", acos(w) - acos(w_2))
            #print()
            #if w <= 1:
               # print(acos(w), winkel, acos(w)-winkel)
            if not np.allclose(w, cos(winkel)):
                print("Falsche drehnung")

            if not np.allclose(sqrt(np.dot(gedreht[i],gedreht[i].T)), sqrt(np.dot(ecken[i], ecken[i].T))):
                print("Längeverschiebung!!!!!!-----------")  

        print("Ecken (Ende) : ", ecken)
        print("gedrehte Ecken : " , np.dot(R,ecken.T).T )
        print("Winkel (echt) : ", acos( np.dot(ecken[i], gedreht[i].T) * 1 / (sqrt(np.dot(gedreht[i],gedreht[i].T))* sqrt(np.dot(ecken[i], ecken[i].T))) ))
        print(" Skalarprodukt : ", np.dot(ecken[i], gedreht[i].T))
        print("Länge anfang : ", np.dot(ecken[i], ecken[i].T))
        print("Länge danach : ", sqrt(np.dot(gedreht[i],gedreht[i].T)))

        return np.dot(R,ecken.T).T  # Ecken mal Rot. Matrix -> neue Ecken ausgegeben

fehler = 0
einser_fehler = 0
ortho_fehler = 0
winkel_fehler = 0

rotation(pi, [(1,1,0,0),(2,0,1,0)], [(1,2,1,1)])

"""for x in (pi, 0): #  np.linspace(0,2*pi,20):
    V = np.array([[random.randint(-500,500) for i in range(4)] for j in range(4)])

    E = np.eye(4)
    #print(" HIer " ,V[0] + V[1] + orthogonal_ebene([V[0],V[1]]).tolist())
    Liste = []
    Liste +=  [V[0].tolist()] + [V[1].tolist()] + orthogonal_ebene([V[0],V[1]]).tolist()
    #print("Liste ", Liste)
    O = np.array(Liste).reshape(4,4)
    #print(O)

    for i in range(4):
        for j in range(4):
            if not np.allclose(O[i]@O[j], E[i,j]) and not (i<2 and j<2):
                ortho_fehler += 1

    V = orthonormalisierung(V)
    
    

    for i in range(4):
        for j in range(4):
            if not np.allclose(np.dot(V[i],V[j]), E[i,j]):
                if i == j:
                    einser_fehler += 1
                fehler += 1

    Rot = rotation(x, [V[0],V[1]], [V[2]]).reshape(1,4)

    if (V[2]@Rot.T / (sqrt(V[2]@V[2].T)*sqrt(Rot@Rot.T))) > 1:
        print((V[2]@Rot.T / (sqrt(V[2]@V[2].T)*sqrt(Rot@Rot.T))))
    elif not np.allclose(acos((V[2]@Rot.T / (sqrt(V[2]@V[2].T)*sqrt(Rot@Rot.T)))), x) :
        winkel_fehler += 1
        print(acos((V[2]@Rot.T / (sqrt(V[2]@V[2].T)*sqrt(Rot@Rot.T))))-x)
            
print(fehler)
print(einser_fehler)
print(fehler-einser_fehler)
print(ortho_fehler)
print(winkel_fehler)"""


a = 234532
b = str(a)
c = list(b)
print(a ,b ,c)