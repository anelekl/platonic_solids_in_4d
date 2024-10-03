
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

    return [np.array([1,0,0,0]), np.array([a, b_2, 0,0]), np.array([a, b_3_4, c_3,0 ]), np.array([a, b_3_4, c_4, d_4])]

def raumwinkel_tetraeder(tetraeder, sample_size, dim):
    sample = list(filter(lambda x:np.linalg.norm(x)<1,2*np.random.random((sample_size,dim))-1))
    
    return len(list(filter(lambda x: all([ np.dot(x,v)< 0 for v in tetraeder ]) ,sample)))/len(sample)

print(raumwinkel_tetraeder(make_tetraeder(0.1),10_000,4))