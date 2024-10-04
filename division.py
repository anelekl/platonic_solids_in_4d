import time
from math import *
import numpy as np

start = time.time()
x = pi
y = 3.141
for i in range(int(7e6)):
    x/y

print(time.time()-start)

print(1/6 * -(1) * np.linalg.det(np.array([[0.39,0.16,-0.77],[-0.85,0.92,0.202],[-0.34,-0.34,0.59]])))
print([1,2,3]*3)