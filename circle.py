import scipy
from scipy.optimize import fsolve
import numpy as np
import math
# x,y is the coordinate, r is the radius

# parameter,dimensional
r = 5e-1
stepPoints = 500
######################################################

pointsList = np.arange(-r, r+ 2.0 * r /
                       stepPoints * 0.99, 2.0 * r / stepPoints)
result = []
for i in pointsList:
    if(i>r):
        i=r
    result.append( math.sqrt(r * r - i * i) )
result = np.array(result)
result = np.abs(result)
print(pointsList,result)
np.savetxt('ansCirPositive.txt', np.array([pointsList, result[:]]).T)
np.savetxt('ansCirNegative.txt', np.array([pointsList, -result[:]]).T)