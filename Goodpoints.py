# Newcrosscircle.R 找到了所有圆两两之间的交集，Goodpoints.R 负责找到出现在所有圆中的点
# This function seeks the "good points" among the crossing points
# found before. These good points represent the intersection
# region of all circles. Good points are any points that are
# located in all circles. These good points allow us to find the
# area where we expect that the target host is located. These
# good points will be the vertices of the polygon. Note that the
# polygon is the approximation of this intersection region. The
# longitude and the latitude of good points are saved in a file
# whose name is the IP address of the given target host and the
# file is located in the folder "Bonpoint". It should be notedd
# that if we have two good points, the intersection region will
# be the circle whose diameter is the geographic distance of
# these two points.
import sys
import os
import numpy as np
from subfunctions import *  

def goodpoints(target, delaysNodes, lonlatNnodes):
    nlandmark = len(delaysNodes)
    nomfic = str.join("Crossingpoint/crossingpoint-", target)
    ficrslt = str.join("Bonpoint/zonecross-", target)
    
    if os.path.exists(ficrslt):
        os.remove(ficrslt)

    with open(nomfic, 'r') as f:
        npoints = np.arange(1, nlandmark+1)
        # seeks the smallest radius.                 
        mindelaysNodes = min(delaysNodes)
        # the least distance between any landmark toward a target.
        pos = delaysNodes.index(mindelaysNodes)
        bondist = delaysNodes[pos]     
        ptcrossing = np.asmatrix(f.read())
        long = len(ptcrossing[:, 0])
        test = 0
        for j in range(0, long):
            compteur = 0
            D1 = calculdist(deg2rad(lonlatNnodes[pos, 0]), deg2rad(lonlatNnodes[pos, 1]), deg2rad(ptcrossing[j, 2]), deg2rad(ptcrossing[j, 3]))
            if D1 <= bondist:
                 # seek if the give b point is inside all circles
                for z in range(0, npoints):
                    D2 = calculdist(deg2rad(lonlatNnodes[z, 0]), deg2rad(lonlatNnodes[z, 1]), deg2rad(ptcrossing[j, 2]), deg2rad(ptcrossing[j, 3]))
                    if D2 <= delaysNodes[z]:
                        compteur += 1
            #如果点j到所有的landmark之间的距离都小于landmark到target的距离，那么该点就属于goodpoint，存入ficrslt，用于最后的绘图
            if compteur == len(npoints):
                with open(ficrslt, 'a') as f:
                    f.write(ptcrossing[j, :])
                    test += 1
        if test == 0:
            with open(ficrslt, 'a') as f:
                    f.write('')

if __name__ == '__main__':
    target = "123.123.123.123"
    delaysNodes = np.array([10, 20, 30])
    lonlatNnodes = np.array([
        [30.01, 50.21],
        [30.01, 50.21],
        [30.01, 50.21]
    ])
    goodpoints(target, delaysNodes, lonlatNnodes)