from math import *
import numpy as np
from matplotlib import pyplot as plt
import random

def dreieck(p1,p2,p3):
    

def dreieck_alt(p1,p2,p3):
    p = [np.array(p1), np.array(p2), np.array(p3)]
    # Translation
    p = np.array([p[i] - np.array(p1) for i in [0,1,2]])
    # Rotation 1
    if p[1][1] != 0:
        alpha = atan(- p[1][2]/p[1][1]) # -c/b
    else:
        alpha = pi/2
    D = np.array([[1,0,0] ,[0, cos(alpha), -sin(alpha)],[0, sin(alpha), cos(alpha)]])
    p = p@D.T # (D@(p.T)).T
    if p[1][0] != 0:
        beta = atan(- p[1][1]/p[1][0] )
    else: 
        beta = pi/2
    D2 = np.array([[cos(beta), -sin(beta), 0], [sin(beta), cos(beta), 0], [0,0,1]])
    p = p@ D2.T
    # Rotation 2
    if p[2][1] != 0:
        gamma = atan(- p[2][2]/p[2][1]) 
    else: 
        gamma = pi/2
    D3 = np.array([[1,0,0] ,[0, cos(gamma), -sin(gamma)],[0, sin(gamma), cos(gamma)]])
    p = p@D3.T
    A = abs(p[1][0] * p[2][1] /2)
    return A

N = 30 # Punkte
def kurve(Q, epsi = 0.01, P = [np.array([0,i/N,0]) for i in range(N+1)], K = 100_000):
    N = len(P)-1
    ev = np.array([epsi, 0, 0])
    #P = [np.array([0,i/N,0]) for i in range(N)]
    A = 0

    for i in range(N):
        A += dreieck(P[i], P[i+1],Q)
    print(A)
    for j in range(K):
        i = random.randint(0,N-2)
        if dreieck(P[i], P[i+1],Q) + dreieck(P[i+1], P[i+2],Q) >= dreieck(P[i] , P[i+1] + ev,Q) + dreieck(P[i+1] + ev, P[i+2],Q):
            P[i+1] += ev
        elif dreieck(P[i], P[i+1],Q) + dreieck(P[i+1], P[i+2],Q) > dreieck(P[i] , P[i+1] - ev,Q) + dreieck(P[i+1] - ev, P[i+2],Q):
            P[i+1] -= ev

    A = 0
    for i in range(N):
        A += dreieck(P[i], P[i+1],Q)
    print(A)
    return [P,A]

def kurve_2(Q1, Q2, epsi,  P = [np.array([0,i/N,0]) for i in range(N+1)], K = 100_000):
    N = len(P)-1
    ev = np.array([epsi, 0, 0])
    A = 0

    for i in range(N):
        A += dreieck(P[i], P[i+1],Q1) + dreieck(P[i], P[i+1],Q2)
    print(A)

    for j in range(K):
        i = random.randint(0,N-2)
        if dreieck(P[i], P[i+1],Q1) + dreieck(P[i+1], P[i+2],Q1) + dreieck(P[i], P[i+1],Q2) + dreieck(P[i+1], P[i+2],Q2) >= dreieck(P[i] , P[i+1] + ev,Q1) + dreieck(P[i+1] + ev, P[i+2],Q1) + dreieck(P[i] , P[i+1] + ev,Q2) + dreieck(P[i+1] + ev, P[i+2],Q2):
            P[i+1] += ev
        elif dreieck(P[i], P[i+1],Q1) + dreieck(P[i+1], P[i+2],Q1) + dreieck(P[i], P[i+1],Q2) + dreieck(P[i+1], P[i+2],Q2) > dreieck(P[i] , P[i+1] - ev,Q1) + dreieck(P[i+1] - ev, P[i+2],Q1) + dreieck(P[i] , P[i+1] - ev,Q2) + dreieck(P[i+1] - ev, P[i+2],Q2):
            P[i+1] -= ev

    A = 0
    for i in range(N):
        A += dreieck(P[i], P[i+1],Q1) + dreieck(P[i], P[i+1],Q2)
    print(A)

    return [P,A]

def exp_verteilung(N_div2, Qx):
    N = 2* N_div2
    lim = (-log(1/N)*2)
    q = [lim*2*i/(N) for i in range(int(N/2))] 
    P = [ Qx*(1 - exp(-q[i])) for i in range(int(N/2))] + [Qx] + [Qx + (1-Qx)*exp(-q[-(i+1)]) for i in range(1,int(N/2))]
    P = [np.array([0,P[i],0]) for i in range(len(P))] 
    #print(len(P))
    return P


# --- START ---
einP = False
zweiP = True

if zweiP:
    Q1 = [1,0,0.2]
    Q_1 = [1,0,0]
    Q2 = [1,1,0.1]
    Q_2 = [1,1,0]

    epsi = 1/(N*50)
    K = 200

    P = kurve_2(Q_1, Q_2 ,epsi=epsi, K=30_000)[0]

    x = [P[i][0] for i in range(len(P))]
    y = [P[i][1] for i in range(len(P))]

    P2 = kurve_2(Q1, Q2 ,epsi=epsi,P=P, K=700_000)[0]

    x2 = [P2[i][0] for i in range(len(P2))]
    y2 = [P2[i][1] for i in range(len(P2))]

    plt.figure()
    plt.scatter(Q1[0], Q1[1])
    plt.scatter(Q2[0], Q2[1])
    plt.plot(x,y)
    plt.plot(x2,y2)
    plt.savefig("2dto3d_rand2_08.png")
    plt.show()


if einP:
    Q = [1,1/3,0]
    Q2 = [1,1/3,0.1]

    epsi = 1/(N*50)
    K = 200

    q = exp_verteilung(20,Q[1])

    P = kurve(Q,epsi, K=100_000)[0]
    #P = copy(P1)

    #P = kurve(Q,epsi,P1, 30_000)[0]

    datei = open('kurve_2d.txt','a')
    datei.write(str(P))

    x = [P[i][0] for i in range(len(P))]
    y = [P[i][1] for i in range(len(P))]


    #P3 = kurve(Q2,epsi=0.1, P=P , K = 4_000 )[0]
    P2 = kurve(Q2,epsi, P, K=100_000 )[0]

    # --- PLOT ---


    x2 = [P2[i][0] for i in range(len(P2))]
    y2 = [P2[i][1] for i in range(len(P2))]

    plt.figure()
    plt.scatter(Q[0], Q[1])
    plt.plot(x,y)
    plt.plot(x2,y2)
    plt.savefig("2dto3d_rand_19.png")
    plt.show()



    #print(P)