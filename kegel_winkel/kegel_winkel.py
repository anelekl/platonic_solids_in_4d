
from math import *
import random 
import numpy as np
from matplotlib import pyplot as plt #für die graphiken

def kegel_verhältnis(dim,sample_size,winkel):
    sample=list(filter(lambda x:np.linalg.norm(x)<1,2*np.random.random((sample_size,dim))-1))
    return len(list(filter(lambda x:acos(x[0]/np.linalg.norm(x))<winkel,sample)))/len(sample)

for dim in range(4,5):
    #Öffnungswnkel des Kegels
    winkel_liste = [0.5003842446] #np.linspace(0,pi/2,100)
    verhältniss = []

    r = 10_000
    anzahl_punkte = 1_000_000
    punkte_zaehler = 0
    punkte = []
    
    for winkel in winkel_liste:
        punkte_zaehler = 0
        punkte = []
    
        innen = 0
        aussen = 0
        
        while punkte_zaehler < anzahl_punkte:
                    p = np.array([random.randint(-r,r)/r]+[ random.randint(-r,r)/r for i in range(1,dim)]).reshape(1,dim)
                    
                    if np.dot(p,p.T) <= 1: #wenn in Späre
                        punkte_zaehler += 1
                        punkte += [p]
        #print(len(punkte), anzahl_punkte)
        for p in punkte:
            #print(p)
            #print([p[0][i]**2 for i in range(p.shape[1])])
            #print( p[0][0] , " / " , sum([p[0][i]**2 for i in range(p.shape[1])]) )
            #print(p[0][0]/ sum([p[0][i]**2 for i in range(p.shape[1])]) )
            if acos(p[0][0]/ sqrt(np.dot(p,p.T)) ) > winkel:
                aussen += 1
                #print("a")
            else: 
                innen += 1
                #print("i")

        verhältniss += [innen/anzahl_punkte/2]
        print(innen/anzahl_punkte/2)
        if innen + aussen != anzahl_punkte:
             print("Verkackt", aussen, innen)
        #print(winkel , "; " , innen/anzahl_punkte)
    #print(verhältniss)
    plt.plot(winkel_liste,verhältniss, label=dim)

#plt.plot(winkel_liste, [(1-cos(winkel_liste[i])) for i in range(len(winkel_liste))], color="black")
#plt.plot(winkel_liste, [(1-cos(winkel_liste[i]))**2 for i in range(len(winkel_liste))], color="lightgreen")
#plt.xlabel("Vergleich: 1-cos schwarz, (1-cos)^2 grün")
plt.legend()
plt.savefig("2024_07_18_09_kegel.png")
plt.show() 

#winkel-raumwinkel dict

winkel_raumwinkel = dict(zip(winkel_liste, verhältniss))

def save_dict_to_file(dic):
    f = open('saved_winkel_raumwinkel_03.txt','w')
    f.write(str(dic))
    f.close()

def load_dict_from_file():
    f = open('saved_winkel_raumwinkel.txt','r')
    data=f.read()
    f.close()
    return eval(data)

#save_dict_to_file(winkel_raumwinkel)