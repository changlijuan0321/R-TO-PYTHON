#以所有地标为圆心，以测得的半径画圆，在两两圆之间判断是否相交，并判断交集的点
########################################################################
    ## The function Newcrosscircle seeks the crossing points between
    ## the circles that have as center the geographic position of the
    ## landmarks and as radius the geographic distance from each
    ## landmark toward a given target host. For each landmark, it
    ## draws the circle whose radius is the geographic distance
    ## between this landmark and the target host, and centered at the
    ## position of the landmark. This function returns a variable. If
    ## the value of this variable is equal to zero we have no crossing
    ## points so we can not localize this target host. If its value is
    ## different to zero, the longitude and the latitude of each
    ## crossing point are saved in a file whose name is the IP address
    ## of the target host and located in the folder "Crossingpoint"
######################################################################
    ## delaysNnodes represents the distance, in km, between a given
    ## landmark and the target. It can be seen as the radius of the
    ## circle where the center is the geographic position
    ## (longitude,latitude), in degrees, of the landmarks. Lonlat
    ## contains the longitude et latitude of our set of landmarks. It
    ## will be use as the center of our circle.
import numpy as np
import math
import os
import sys
from subfunctions import *  

def newcrosscircle(target, delaysNnodes, lonlatNnodes):
    nlandmarks = len(delaysNnodes)
    veclines = list(plotcircle(lonlatNnodes[0, 0], lonlatNnodes[0, 1], delaysNnodes[0]))
    for j in range(0, nlandmarks):
        matline = list(plotcircle(lonlatNnodes[j, 0], lonlatNnodes[j, 1], delaysNnodes[j]))
        veclines = np.c_[veclines, matline]

    L = len(veclines)
    lineA = np.mat(np.zeros((101, 2)))
    lineB = np.mat(np.zeros((101, 2))) 
    #nomfic用来收集所有圆的交集
    nomfic = "Crossingpoint/crossingpoint"+target
    f = open(nomfic, 'w')
    kl = 0
    indice = 0
    bon = 0
    #"bon"表示是否有交叉点

    while indice <= L:
        kl += 1
        h = kl
        compteur = indice
        compteur += 2
        while compteur <= L and h < nlandmarks:
            #第一个圆的交点和半径
            R1 = delaysNnodes[kl-1]
            lon1 = lonlatNnodes[kl-1, 0]
            lat1 = lonlatNnodes[kl-1, 1]
            h += 1
            #第二个圆的交点和半径
            R2 = delaysNnodes[h-1]
            lon2 = lonlatNnodes[h-1, 0]
            lat2 = lonlatNnodes[h-1, 1]
            #形成第一个圆第二个圆的交点
            lineA = veclines[:, indice:(indice + 1)]
            lineB = veclines[:, compteur:(compteur + 1)]
            compteur = compteur+2
            #测试这两个圆是否相交
            ptsegments1 = []
            ptsegments2 = []
            if cross(lon1, lat1, R1, lon2, lat2, R2) == 1:
                # find segment's first circle that cross the second circle 
                ptsegments1 = findsegments(lineA, R2, lon2, lat2)
                #find segment's second circle that cross first circle
                ptsegments2 = findsegments(lineB, R1, lon1, lat1)

            for comp in range(0, len(ptsegments1)):
                if ptsegments1[comp,0] != 0:
                    #计算lineB与lineA的相交一条线段eq1
                    eq1 = eqdroite(ptsegments1[comp, 0], ptsegments1[comp, 1], ptsegments1[comp, 2], ptsegments1[comp, 3])
                    for ind in range(0,ptsegments2):
                        #找到lineB与lineA的一条线段eq2
                        if ptsegments2[ind, 0]!=0:
                            eq2 = eqdroite(ptsegments2[ind, 0], ptsegments2[ind, 1], ptsegments2[ind, 2], ptsegments2[ind, 3])
                            #计算两个直线的交点ptcross
                            ptcross = pointinter(eq1, eq2, ptsegments2[comp, 1], ptsegments2[comp, 3], ptsegments2[ind, 1], ptsegments2[ind, 3])
                            #两个线段不一定在圆内相交，因为eq1、eq2只是两个圆相交的线段，但不表示这两者就是对应的，所以如果相交就存入之前新建的nomfic文件，不然就继续
                            if int(ptcross):
                                D1 = calculdist(deg2rad(ptcross[1]), deg2rad(ptcross[2]), deg2rad(lon1), deg2rad(lat1))
                                D2 = calculdist(deg2rad(ptcross[1]), deg2rad(ptcross[2]), deg2rad(lon2), deg2rad(lat2))
                                if(D1 <= R1 and D2 <= R2):
                                    if(ptcross[1] >= -180 and ptcross[1] <= 180) and (ptcross[2] >= -90 and ptcross[2] <= 90):
                                        f.write(str(ptcross[1: -1]))
                                        bon += 1
        indice += 2
    return bon

if __name__ == '__main__':
    target = '60.247.13.2'
    delaysNodes = np.array([10, 20, 30])
    lonlatNnodes = np.array([
        [30.01, 50.21],
        [30.01, 50.21],
        [30.01, 50.21]
    ])
    newcrosscircle(target, delaysNodes, lonlatNnodes)




