import random
import numpy as np
#from matplotlib import pyplot as plt
import time
from copy import deepcopy
import sys

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

#Berechnet raum1\raum2(im Sinne der Mengendifferenz)
def diff(raum1,raum2):
    mat=orthonormalisierung(raum2)
    k=len(mat)
    return orthonormalisierung(np.concatenate((mat,raum1)))[k:]

#nimmt eine Liste von Vektoren, gibt eine Orthonormalbasis des zum Spann der
#Vektoren orthogonalen Raums (als Zeilenvektoren einer Matrix)
def orthogonal_raum(vektoren, dim:int = 4) -> np.ndarray:
    if(len(vektoren)>0):
        dim=len(vektoren[0])
    else:
        return np.eye(dim)
    return diff(np.eye(dim),vektoren)

def orthogonal_ebene(vektoren,dim:int=4) -> np.ndarray:
    if(len(vektoren)>0):
        dim=len(vektoren[0])
    assert(len(vektoren)+2 == dim)
    return orthogonal_raum(vektoren,dim)

#Gibt Matrix zurück, die eine Rotation um den Winkel in der aus den ersten bei-
#den Basisvektoren gebildeten Ebene beschreibt, so dass ein Winkel von 90 Grad
#bei einer Multiplikation von links den zweiten in den ersten Basisvektor über-
#führt.
def rot_matrix(winkel:float,basis:np.ndarray=np.eye(4)) -> np.ndarray:
    dim = len(basis)
    assert(dim >= 2)
    mat2d = np.array([[np.cos(winkel),np.sin(winkel)],
                    [-np.sin(winkel),np.cos(winkel)]])
    return basis.T@np.concatenate(
           (np.concatenate((mat2d,np.zeros((2,dim-2))),axis=1),np.eye(dim)[2:])
           )@basis

def rotation(winkel:float,ebene,ecken,offset=None) -> np.ndarray:
    if type(offset)==type(None):
        offset=np.zeros((1,len(ebene[0])))
    vektoren = orthonormalisierung(ebene)
    assert(len(vektoren) == 2)
    basis = np.concatenate((vektoren,orthogonal_raum(vektoren)))
    R = rot_matrix(winkel,basis)
    return offset+(ecken-offset)@R.T

def spiegelung(normale, punkte,offset=None):
    if type(offset)==type(None):
        offset=np.zeros((1,len(normale)))
    normale/=np.linalg.norm(normale)
    return punkte-(punkte-offset)@normale*normale*2
def spiegelung_gen(spiegel_raum,punkte,offset=None):
    assert(len(punkte)>0)
    if type(offset)==type(None):
        offset=np.zeros((1,len(punkte[0])))
    normal_raum=orthogonal_raum(spiegel_raum)
    return punkte-2*(punkte-offset)@normal_raum.T@normal_raum

def drei_zu_vier(punkte):
    return np.concatenate((punkte,np.zeros((len(punkte),1))),axis=1)

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
    return np.acos(vektor1@vektor2)

def schnitt(raum1,raum2):
    return orthogonal_raum(np.concatenate((raum1,raum2)).T)[:,:len(raum1)]@raum1

class plat_solid:
    def __init__(self,schäfli_list,higher_order=None,vec_ortho=None):
        if higher_order!=None:
            self.points=schäfli_list
            self.higher_order=higher_order
            self.vec_ortho=vec_ortho
            return
        if len(schäfli_list)==1:
            num_points=schäfli_list[0]
            self.angle=2*np.pi/num_points
            matrix_rot=np.array([[np.cos(self.angle),-np.sin(self.angle)],[np.sin(self.angle),
                                                                np.cos(self.angle)]])
            self.points=[np.array([0,0]),np.array([1,0])]
            self.vec_ortho=[np.array([0,1])]
            for _ in range(num_points-2):
                self.points.append(self.points[-1]+matrix_rot@(self.points[-1]-
                                                               self.points[-2]))
                self.vec_ortho.append(matrix_rot@self.vec_ortho[-1])
            self.vec_ortho.append(matrix_rot@self.vec_ortho[-1])
            self.higher_order=[[{i} for i in range(num_points)],
                               [{i,(i+1)%num_points} for i in range(num_points)],
                               [{i for i in range(num_points)}]]
            self.angle_in=np.pi-self.angle
            return
        dim=len(schäfli_list)+1
        a=plat_solid(schäfli_list[:-1]).copy(dim)
        parts=[a]
        self.vec_ortho=[np.eye(dim)[-1]]
        #winkel herausfinden (TODO)
        angle_rot=2*np.pi/schäfli_list[-1]
        boundary_m3=a.get_points(dim-3,0)
        mat_ortho=orthogonal_raum(boundary_m3[1:]-boundary_m3[0]) if len(boundary_m3)>1 else orthogonal_raum([],3)
        mat_ortho=diff(mat_ortho,[np.eye(dim)[-1]])
        b=a.copy()
        c=a.copy()
        c.points=list(rotation(angle_rot,mat_ortho,c.points,boundary_m3[0]))
        dir_b=sum(boundary_m3)/len(boundary_m3)-sum(b.points)/len(b.points)
        dir_c=sum(boundary_m3)/len(boundary_m3)-sum(c.points)/len(c.points)
        dir_b/=np.linalg.norm(dir_b)
        dir_c/=np.linalg.norm(dir_c)
        boundary_m2_b=min([i for (i,x) in enumerate(b.higher_order[dim-2]) if 0 in x])
        boundary_m2_c=max([i for (i,x) in enumerate(c.higher_order[dim-2]) if 0 in x])
        boundary_m3_indices=a.reduce_order(dim-3,0,0)
        point_b=b.points[min([x for x in b.reduce_order(dim-2,0,boundary_m2_b) if not x in boundary_m3_indices])]
        vec_b=dir_b@(boundary_m3[0]-point_b)*dir_b
        points_c=[x for (i,x) in enumerate(c.points) if i in c.reduce_order(dim-2,0,boundary_m2_c) and not i in boundary_m3_indices]
        for point_c in points_c:
            vec_c=dir_c@(boundary_m3[0]-point_c)*dir_c
            mat=np.stack((vec_b,vec_c,point_b-point_c),axis=1)
            ortho_basis=orthogonal_raum(mat)
            if len(ortho_basis)!=1:
                print("huh!")
                exit(1)
            if abs(ortho_basis[0][0]+ortho_basis[0][1])>1e-5:
                continue
            coefs=ortho_basis[0]/ortho_basis[0][2]
            angle_meet=np.acos(1+coefs[1])
            break
        b.rotate(dim-3,0,angle_meet,np.eye(dim)[-1])
        c.rotate(dim-3,0,angle_meet,np.eye(dim)[-1])
        boundary=b.get_points(dim-2,boundary_m2_b)
        vec_b=sum(b.points)/len(b.points)-sum(boundary)/len(boundary)
        vec_c=sum(c.points)/len(c.points)-sum(boundary)/len(boundary)
        vec_b/=np.linalg.norm(vec_b)
        vec_c/=np.linalg.norm(vec_c)
        self.angle_in=np.acos(vec_b@vec_c)
        self.angle=np.pi-self.angle_in
        
        #körper konstruieren
        self.points=a.points.copy()
        self.higher_order=deepcopy(a.higher_order)
        unmatched=self.get_unmatched(schäfli_list[-1])
        while len(unmatched)>0:
            to_be_completed=unmatched[0]
            to_be_mirrored=adjacent=0
            for i in range(len(self.higher_order[dim-2])):
                if not to_be_completed in self.higher_order[dim-2][i]:
                    continue
                adj_adjacent=[j for j in range(len(self.higher_order[dim-1]))
                                if i in self.higher_order[dim-1][j]]
                if len(adj_adjacent)==0:
                    sys.exit(1)
                if len(adj_adjacent)==2:
                    continue
                adjacent=i
                to_be_mirrored=adj_adjacent[0]
                break
            to_be_added=parts[to_be_mirrored].copy()
            points_target=self.get_points(dim-2,adjacent)
            index_adapted=0
            for i in range(len(to_be_added.higher_order[dim-2])):
                points=to_be_added.get_points(dim-2,i)
                if np.linalg.norm(sum(points)-sum(points_target))<1e-5:
                    index_adapted=i
                    break
            to_be_added.flip(dim-2,index_adapted)
            point_out=sum(to_be_added.points)/len(to_be_added.points)-sum(points_target)/len(points_target)
            point_out/=np.linalg.norm(point_out)
            to_be_added.points=list(rotation(self.angle,[self.vec_ortho[to_be_mirrored],point_out],to_be_added.points,points_target[0]))
            self.vec_ortho.append(rotation(self.angle,[self.vec_ortho[to_be_mirrored],point_out],[self.vec_ortho[to_be_mirrored]])[0])
            index_change=[[]]
            for i in range(len(to_be_added.points)):
                isIn=False
                for j in range(len(self.points)):
                    if np.linalg.norm(to_be_added.points[i]-self.points[j])>1e-5:
                        continue
                    index_change[0].append(j)
                    isIn=True
                    break
                if isIn:
                    continue
                index_change[0].append(len(self.points))
                self.points.append(to_be_added.points[i])
                self.higher_order[0].append({index_change[0][-1]})
                continue
            for (i,order) in enumerate(to_be_added.higher_order):
                if i==0:
                    continue
                index_change.append([])
                for (j,obj) in enumerate(order):
                    obj_changed={index_change[i-1][k] for k in obj}
                    if obj_changed in self.higher_order[i]:
                        index_change[i].append(self.higher_order[i].index(obj_changed))
                        continue
                    index_change[i].append(len(self.higher_order[i]))
                    self.higher_order[i].append(obj_changed)
            parts.append(to_be_added)
            unmatched=self.get_unmatched(schäfli_list[-1])
        self.higher_order.append([{i for i in range(len(self.higher_order[dim-1]))}])
        
    def copy(self,dim=None):
        if dim==None:
            dim=len(self.points[0])
        assert(dim>=len(self.points[0]))
        return plat_solid(list(np.concatenate((self.points,np.zeros((len(self.points),dim-len(self.points[0])))),axis=1)),deepcopy(self.higher_order),list(np.concatenate((self.vec_ortho,np.zeros((len(self.vec_ortho),dim-len(self.vec_ortho[0])))),axis=1)))
    def reduce_order(self,order,orderNew,index:int):
        objects_included={index}
        for order_curr in range(order-1,orderNew-1,-1):
            objects_included={i for i in range(len(self.higher_order[order_curr])) if any([i in self.higher_order[order_curr+1][j] for j in objects_included])}
        return objects_included
    #Gibt die Liste der Punkte zurück, welche im index-ten order-dimensionalen Teilkörper enthalten sind.
    def get_points(self,order,index):
        return [self.points[i] for i in self.reduce_order(order,0,index)]
    def flip(self,order,index):
        plane=self.get_points(order,index)
        plane_basis=plane[1:]-plane[0]
        self.points=list(spiegelung_gen(plane_basis,self.points,plane[0]))
        self.vec_ortho=list(spiegelung_gen(plane_basis,self.vec_ortho))
    def rotate(self,order,index,angle,dir):
        points_staying=self.get_points(order,index)
        vec_out=sum(self.points)/len(self.points)-sum(points_staying)/len(points_staying)
        vec_out/=np.linalg.norm(vec_out)
        self.points=list(rotation(angle,[dir,vec_out],self.points,sum(points_staying)/len(points_staying)))
        self.vec_ortho=list(rotation(angle,[dir,vec_out],self.points))
    def get_unmatched(self,target):
        dim=len(self.points[0])
        return [x for x in range(len(self.higher_order[dim-3])) if len([1 for z in self.higher_order[dim-1] if any([x in self.higher_order[dim-2][y] for y in z])])!=target]

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
if __name__=="__main__":
    körper_drehung = False ## körper_drehung macht das Winkel finden

    if körper_drehung:
        # Funktioniert erstmal nur für Eckenfigur <3,3>

        #3 Ecken: die erste, die "mittlere", die an der alles aufeinander trifft: NICHT MIT ANGEBEN -- weil die in 0 gesetzt wurde
        # ich glaube, die sind alle unnötig, weil sie nicht mehr gebraucht werden
        a = - (np.sqrt(5)+1)/4
        c = (1 - (np.sqrt(5) + 1)/2)/2
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

        Winkel = np.linspace( 0,np.pi,1001) # winkel, die ausprobiert werden, wäre super, wenn man das durch Variablen irgenwo eingeben könnte
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
    Ebenen_Matrix = np.array([[(1, 0, 0, 0), (0.5, np.sqrt(3)/2, 0, 0) , (0.5, np.sqrt(3)/6, np.sqrt(2/3), 0)]]).reshape(a,dim)

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
