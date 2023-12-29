
from math import *
import numpy as np

def side_len(n):
    return pi -((2*pi)/n) 

def winkel_kugel(a, b, c):
    return acos((cos(a) - cos(b) * cos(c)) / (sin(b) * sin(c)))

def area(innenwinkel):
    return sum(innenwinkel) - (len(innenwinkel) - 2)*pi

def raumwinkel_seiten(seiten):
    winkel = [winkel_kugel(seiten[i],seiten[(i+1)%3],seiten[(i+2)%3]) for i in range(3)]
    return area(winkel)
    
def raumwinkel(args):
    (x,y,z)=args
    seiten = [side_len(i) for i in [x,y,z]]
    return raumwinkel_seiten(seiten)

archimedian = [(3, 3, 3), (4, 4, 4), (5, 5, 5), (6, 6, 3), (8, 8, 3), (10, 10, 3), (6, 6, 5), (6, 6, 4), (4, 6, 8), (4, 6, 10)]

# 3333: 0.4326938
# 3344: 0.7836531040612146
# 3355: 1.1692084617955367
# 3444: 1.1081734479693928
# 3454: 1.4153040904516636
#33333: 0.8386023640061511

number_list = []
for i in archimedian:
    print(i, ':', raumwinkel(i) / pi, "pi")
    number_list.append(raumwinkel(i) / pi)

for i in range(3,11):
    print('(', i, ', 4 , 4 ) :' , raumwinkel((i,4,4))/ pi , "pi")

Winkel = np.linspace(0,2*pi,200)
Daten  = []
for i in Winkel:
    if acos(-1/sqrt(2+tan(i)**2)) == 0:
        Daten.append(0)
    else:    
        Daten.append(raumwinkel_seiten((pi/2,acos(-1/sqrt(2+tan(i)**2)),acos(-1/sqrt(2+tan(i)**2)))))

with open('el_20231229_07.txt', "a") as speicher:
    print("; \n", np.array(Daten).reshape(len(Winkel),1).tolist(), "\n ;", Winkel,  "\n \n", file=speicher)

arch_group = {(3, 3, 3) : 0.17547965609182192,
               (4, 4, 4) : 0.5,
               (5, 5, 5) : 0.9427508529512999,
               (6, 6, 6) : 0.608173447969393, 
               (8, 8, 3) : 0.8918265520306073,
               (10, 10, 3) : 1.2322795271987699,
               (6, 6, 5) : 1.3524163823495672,
               (6, 6, 4) : 1.0000000000000007, 
               (4, 6, 8) : 1.2500000000000002,
               (4, 6, 10) : 1.5,
               (3, 3, 3, 3) : 0.4326938,
               (3, 3, 4, 4) : 0.7836531040612146,
               (3, 3, 5, 5) : 1.1692084617955367,
               (3, 4, 4, 4) : 1.1081734479693928,
               (3, 4, 5, 4) : 1.4153040904516636,
               (3, 3, 3, 3, 3) : 0.8386023640061511}

"""
def find_combinations(numbers, target_sum, start_index=0, partial=[]):
    count=0
    current_sum = sum(partial)

    if round(current_sum, 4) == target_sum:
        count+=1
        print(partial)
        print([list(arch_group.keys())[list(arch_group.values()).index(n)] for n in partial])
        print()
    elif current_sum > target_sum:
        return count

    for i in range(start_index, len(numbers)):
        count += find_combinations(numbers, target_sum, i, partial + [numbers[i]])



count=find_combinations(list(arch_group.values()), 4)
print("Anzahl:",count)

"""
