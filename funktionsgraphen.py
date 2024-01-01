
from math import *
import numpy as np
from matplotlib import pyplot as plt #für die graphiken

plt.figure(figsize=(12,8))

farben = [ 'dimgray' , 'silver', 'lightcoral', 'firebrick', 'darkred' ,'red', 'orangered' , 'darkorange', 'gold', 'yellow' , 'yellowgreen' , 'limegreen', 'forestgreen' , 'darkgreen' , 'teal' , 'royalblue' , 'navy' , 'indigo']
x = np.linspace(0,pi,200)

farbe = 0
for w in np.linspace(0.0001,pi-0.0001,10):
    farbe += 1
    y_werte_plus = []
    y_werte_minus = []
    for i in x:
        y_werte_plus.append(acos((sin(i)**2)*cos(w) + cos(i)**2))
        y_werte_minus.append(acos(-(sin(i)**2)*cos(w) - cos(i)**2))

    plt.plot(x, y_werte_plus , color = farben[farbe])
    plt.plot(x, y_werte_minus , color = farben[farbe] , linestyle = 'dashed')

plt.xlabel('Winkel')
plt.ylabel('Fläche')
plt.show()