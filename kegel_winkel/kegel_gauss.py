
from math import *
import numpy as np
from matplotlib import pyplot as plt #für die graphiken


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
    raumwinkel.pop(key_list[-i])

print(gauss_winkel)
print(raumwinkel)

plt.plot(raumwinkel.keys(), raumwinkel.values())
plt.plot(gauss_winkel.keys(), raumwinkel.values())
plt.show()