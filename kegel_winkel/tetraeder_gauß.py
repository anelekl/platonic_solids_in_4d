
from math import *
import numpy as np
import random

def make_tetraeder(winkel):
    a = cos(winkel)
    b_2 = sqrt(1 - a**2)
    b_3_4 = (a - a**2)/b_2
    c_3 = sqrt(1 - a**2 - b_3_4**2)
    c_4 = (a - a**2 - b_3_4**2) / c_3
    d_4 = sqrt(1 - a**2 - b_3_4**2 - c_4**2 )

    # Normalenvektoren der 4 3-Ebenen, die den Tetraeder bilden
    return [np.array([1,0,0,0]), np.array([a, b_2, 0,0]), np.array([a, b_3_4, c_3,0 ]), np.array([a, b_3_4, c_4, d_4])]

def raumwinkel_tetraeder(tetraeder, sample_size, dim):
    sample = list(filter(lambda x:np.linalg.norm(x)<1,2*np.random.random((sample_size,dim))-1))
    
    # Raumwinkel durch Normalenvektoren
    return len(list(filter(lambda x: all([ np.dot(x,v)< 0 for v in tetraeder ]) ,sample)))/len(sample)

print(raumwinkel_tetraeder(make_tetraeder(0.1),10_000,4))

def orthogonal_raum(vektoren, dim:int = 4) -> np.ndarray:
    if(len(vektoren)>0):
        dim=len(vektoren[0])
    else:
        return np.eye(dim)
    return diff(np.eye(dim),vektoren)

def schnitt(raum1,raum2):
    return orthogonal_raum(np.concatenate((raum1,raum2)).T)[:,:len(raum1)]@raum1

def projektion(punkt, ebene):
    ortho = orthogonal_raum(ebene)
    for i in range(len(ortho)):
        punkt = punkt - ortho[i]*np.dot(punkt,ortho[i])
    return punkt

def dreieck_winkel(tetraeder, ecken_index, andere_ecken_index):
    j = ecken_index
    ecke = tetraeder[j] 
    ortho = orthogonal_raum(ecke)
    raum = [tetraeder[ecken_index], tetraeder[andere_ecken_index[0]], tetraeder[andere_ecken_index[1]]]
    ebene = schnitt(ortho,raum)
    ecken_neu = []

    for i in j + andere_ecken_index:
        ecken_neu += projektion(ebene, tetraeder[i])
    
    # np.dot(1-0 , 2-0 )
    return acos( np.dot(ecken_neu[1]-ecken_neu[0], ecken_neu[2]-ecken_neu[0]) / (np.linalg.norm(ecken_neu[1]-ecken_neu[0]) * np.linalg.norm(ecken_neu[2]-ecken_neu[0])) ) 

tetra = [[1,0,0,0],[0,1,0,0],[0,0,1,0],[sqrt(1/2), sqrt(1/2),0,0]]
print( dreieck_winkel(tetra, 0, [1,2]))