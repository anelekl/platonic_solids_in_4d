
from math import *
import numpy as np
from matplotlib import pyplot as plt #für die graphiken

#hol
def load_dict_from_file():
    f = open('saved_winkel_raumwinkel.txt','r')
    data=f.read()
    f.close()
    return eval(data)

winkel_liste = np.linspace(0,pi,200)
raumwinkel = load_dict_from_file()

gauss_winkel = {}
key_list = list(raumwinkel.keys())


for i in range(int(len(raumwinkel)/2)):
    gauss_winkel[pi/2-winkel_liste[i]] = raumwinkel[winkel_liste[i]]
    raumwinkel.pop(key_list[-i-1]) 

#print(gauss_winkel)
#print(raumwinkel)

#plt.plot(list(raumwinkel.keys()), [ (w) for w in list(raumwinkel.values())])
# 3,3,3  3,4,3  4,3,3  5,3,3  3,3,4  3,3,5
#plt.scatter([1.394403364, 0.8644751331, 1.047197551, 0.3400114323, 1.212200941, 0.5003842446 ], [5,24, 16, 600, 8, 120])
plt.plot(raumwinkel.keys(), raumwinkel.values())
plt.plot(gauss_winkel.keys(), raumwinkel.values())
plt.plot(raumwinkel.keys(), [(1-cos(w))/2 for w in list(raumwinkel.keys()) ], color="green" )
#plt.savefig("2024_19_08_pidurchwinkel.png")
plt.show()
