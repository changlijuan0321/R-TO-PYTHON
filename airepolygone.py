import numpy as np
from subfunctions import *  
from scipy.spatial import ConvexHull
import math

def airepolygone(matrix, pointenKM):
    long = len(matrix[ :1])
    aire = np.mat(np.zeros(long, 1))
    for j in range(1, long-1):
        aire[j] = pointenKM[j, 1] * pointenKM[(j+1), 2] - pointenKM[(j+1), 1] * pointenKM[j, 2]
    aire[long] = pointenKM[long,1] * pointenKM[1, 2] - pointenKM[1, 1] * pointenKM[long, 2]

    #计算多边形的面积
    nbelement = len(scipy.spatial.ConvexHull(matrix))

    if nbelement == long:
        AIRE = (1/2) * sum(aire)
    else:
        AIRE = -1

    matCx = np.mat(np.zeros(long, 1))
    matCy = np.mat(np.zeros(long, 1))

    for h in range(1, long-1):
        matCx[h] = (pointenKM[h, 1] + pointenKM[(h+1), 1]) * (pointenKM[h, 1] * pointenKM[(h+1), 2] - pointenKM[(h+1), 1] * pointenKM[h, 2])
        matCy[h] = (pointenKM[h, 2] + pointenKM[(h+1), 2]) * (pointenKM[h, 1] * pointenKM[(h+1), 2] - pointenKM[(h+1), 1] * pointenKM[h, 2])

    matCx[long] = (pointenKM[long, 1]+pointenKM[1, 1]) * ((pointenKM[long, 1] * pointenKM[1, 2])  -pointenKM[1, 1] * pointenKM[long, 2])
    matCy[long] = (pointenKM[long, 2]+pointenKM[1, 2]) * ((pointenKM[long, 1] * pointenKM[1,2]) - pointenKM[1, 1] * pointenKM[long,2])

    #计算体心坐标
    Cx = sum(matCx)
    Cy = sum(matCy)
    Cx = (1/(6*AIRE)*Cx)
    Cy = (1/(6*AIRE)*Cy)

    #convert Cx and Cy to degree
    dlon = pointenKM[(long + 1), 1]
    dlat = pointenKM[(long + 1), 1]
    Cx = Cx/dlon
    Cy = Cy/dlat

    R = np.sqrt(abs(AIRE)/math.pi)
    rslt = [Cx, Cy, R, abs(AIRE)]
    return rslt


