import scipy
from scipy.optimize import fsolve
import numpy as np
import math

# x,y is the coordinate, b is the half disctance between source and sink, ratio is Sigma/Vinf


def rankineShape(y, paras):
    x, ratio, b = paras
    return y + ratio / 2.0 / math.pi * (math.atan2(y, x + a) - math.atan2(y, x - a))


# parameter,dimensional
flowRate = 2
a = 1
stepPoints = 500
######################################################

velocityInf = 1.0
Ratio = flowRate / velocityInf
Range = math.sqrt(a * a + Ratio * a / math.pi)
pointsList = np.arange(-Range, Range + 2.0 * Range /
                       stepPoints - 1e-10, 2.0 * Range / stepPoints)
result = []
for i in pointsList:
    paras = [i, Ratio, a]
    # a is the initial value in fsolve function
    result.append(scipy.optimize.fsolve(rankineShape, a, args=paras))
result = np.array(result)
result = np.abs(result)
np.savetxt('ansDinmensional.txt', np.array([pointsList, result[:, 0]]).T)
# NonDimensioanal
# See website xuzhaoyue1995.github.io
paras = [0.0, Ratio, a]
h = scipy.optimize.fsolve(rankineShape, a, args=paras)
a2 = a / (2.0 * h)
flowRate = a2 / a * flowRate
Ratio = flowRate / velocityInf
Range = math.sqrt(a2 * a2 + Ratio * a2 / math.pi)
pointsList = np.arange(-Range, Range + 2.0 * Range /
                       stepPoints - 1e-10, 2.0 * Range / stepPoints)
result2 = []
for i in pointsList:
    paras = [i, Ratio, a2]
    result2.append(scipy.optimize.fsolve(rankineShape, a2, args=paras))
result2 = np.array(result2)
result2 = np.abs(result2)
result2[0]=0.0
result2[-1]=0.0
halfMap2=np.array([pointsList, result2[:, 0], np.zeros(np.shape(pointsList))]).T
#totalMap=pointsList[-2:0:-1] maybe useful
np.savetxt('ansNonDinmensional.txt', halfMap2)
halfMap2[:,1]=-halfMap2[:,1]
np.savetxt('ansNonDinmensionalNegative.txt', halfMap2)
paras = [0.0, Ratio, a2]
h = scipy.optimize.fsolve(rankineShape, a, args=paras)
print('The nondimensional X-axis a is', a2, 'The nondimensional source strength is',
      flowRate, 'The long axis is ', Range, 'The short axis is ', h)
