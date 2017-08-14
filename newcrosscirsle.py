import numpy as np
import math
import os
import sys
from subfunctions import *  

def newcrosscircle(target, delaysNnodes, lonlatNnodes):
    nlandmarks = len(delaysNnodes)
    veclines = plotcircle(lonlatNnodes[1, 1], lonlatNnodes[1, 2], delaysNnodes[1])
    for j in range(2,nlandmarks):
        matline = plotcircle(lonlatNnodes[j, 1], lonlatNnodes[j, 2], delaysNnodes[j])
        veclines = np.vstack(veclines, matline)

    L = len(veclines)
    lineA = np.mat(np.zeros((101, 2)))
    lineB = np.mat(np.zeros((101, 2)))

    nomfic = str.join("Crossingpoint/crossingpoint",target)
    
    with open('nomfic.txt') as f:
        f.read()

    kl = 0
    indice = 1
    bon = 0
    #"bon"表示是否有交叉点

    while indice <= L:
        kl = kl+1
        h = kl
        compteur = indice
        compteur = compteur+2
        while compteur <= L and h < nlandmarks:
            R1 = delaysNnodes[kl]
            lon1 = lonlatNnodes[kl, 1]
            lat1 = lonlatNnodes[kl, 2]

            h = h+1
            R2 = delaysNnodes[h]
            lon2 = lonlatNnodes[h, 1]
            lat2 = lonlatNnodes[h, 2]

            lineA = np.vstack[veclines, indice:(indice+1)]
            lineB = np.vstack[veclines, compteur:(compteur+1)]
            #分配两个matrix存两个圆

            compteur = compteur+2

            #测试这两个圆是否相交
            if cross(lon1, lat1, R1, lon2, lat2, R2) == 1:
                ptsegments1 = findsegments(lineA, R2, lon2, lat2)
                ptsegments2 = findsegments(lineB, R1, lon1, lat1)

            for comp in range(1, ptsegments1):
                if ptsegments1[comp,1] != 0:
                    #计算lineB与lineA的相交一条线段eq1
                    eq1 = eqdroite(ptsegments1[comp, 1], ptsegments1[comp, 2], ptsegments1[comp, 3], ptsegments1[comp, 4])
                    for ind in range(1,ptsegments2):
                        #找到lineB与lineA的一条线段eq2
                        if ptsegments2[ind, 1]!=0:
                            eq2 = eqdroite(ptsegments2[ind, 1], ptsegments2[ind, 2], ptsegments2[ind, 3], ptsegments2[ind, 4])
                            #计算两个直线的交点ptcross
                            ptcross = pointinter(eq1, eq2, ptsegments2[comp, 1], ptsegments2[comp, 3], ptsegments2[ind, 1], ptsegments2[ind, 3])
                            #两个线段不一定在圆内相交，因为eq1、eq2只是两个圆相交的线段，但不表示这两者就是对应的，所以如果相交就存入之前新建的nomfic文件，不然就继续
                            if ptcross.isdigit:
                                D1 = calculdist(deg2rad(ptcross[1]), deg2rad(ptcross[2]), deg2rad(lon1),deg2rad(lat1))
                                D2 = calculdist(deg2rad(ptcross[1]), deg2rad(ptcross[2]), deg2rad(lon2),deg2rad(lat2))
                                if(D1 <= R1 and D2 <= R2):
                                    if(ptcross[1]>=-180 and ptcross[1]<=180) and (ptcross[2]>=-90 and ptcross[2]<=90):
                                        f.write(ptcross)
                                        bon = bon+1
    indice = indice+2
    return bon






