import scipy
from scipy.optimize import fsolve
import numpy as np
import math

# x,y is the coordinate, b is the half disctance between source and sink, ratio is Sigma/Vinf


def rankineShape(y, paras):
    x, ratio, b = paras
    return y - ratio / 2.0 / math.pi * math.atan((2.0 * b * y) / (x * x + y * y - b * b))
#    return y + ratio / 2.0 / math.pi * (math.atan2(y,x+b)-math.atan2(y,x-b))

flowRate = 20.0
velocityInf = 1.0
b = 6.0
stepPoints = 50
Ratio = flowRate / velocityInf
Range = math.sqrt(b * b + Ratio * b / math.pi)
pointsList = np.arange(-Range, Range + 2.0 * Range /
                       stepPoints - 1e-10, 2.0 * Range / stepPoints)
result = []
for i in pointsList:
    paras = [i, Ratio, b]
    result.append(scipy.optimize.fsolve(rankineShape, 7.0, args=paras))
result = np.array(result)
result = np.abs(result)

print(result)
np.savetxt('ans.dat', np.array([pointsList,result[:,0]]).T)
