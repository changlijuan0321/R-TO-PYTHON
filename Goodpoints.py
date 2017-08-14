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
    
    with open('nomfic.txt') as f:
        npoints = np.array(1,nlandmark)  
        #处理可能有问题
        mindelaysNodes = min(delaysNodes)
        pos = delaysNodes.index(mindelaysNodes)
        bondist = delaysNodes[pos]     
        ptcrossing = np.asmatrix(f.read())
        #处理可能有问题
        long = len(ptcrossing[:,1])
        test = 0
        for j in range(1, long):
            compteur = 0
            D1 = calculdist(deg2rad(lonlatNnodes[pos, 1]), deg2rad(lonlatNnodes[pos, 2]), deg2rad(ptcrossing[j, 1]), deg2rad(ptcrossing[j, 2]))

            if D1 <= bondist:
                for z in range(1,npoints):
                    D2 = calculdist(deg2rad(lonlatNnodes[z, 1]), deg2rad(lonlatNnodes[z, 2]), deg2rad(ptcrossing[j, 1]), deg2rad(ptcrossing[j, 2]))
                    if D2 <= delaysNodes[z]:
                        compteur =  compteur +1
            #如果点j到所有的landmark之间的距离都小于landmark到target的距离，那么该点就属于goodpoint，存入ficrslt，用于最后的绘图
            if compteur == len(npoints):
                with open('ficrslt.txt') as f:
                    f.write(ptcrossing[j, :])
                    test = test + 1
        if test == 0:
            with open('ficrslt.txt') as f:
                    f.write('')

