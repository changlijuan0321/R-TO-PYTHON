## This function provide the localisation of a target without
## plotting the points in a figure. "delaysNnodes" contains the
## geographical distance constraint between landmarks and target
## hosts. "lonlatNnodes" contains the position of all succesfull
## landmarks. "target" the names of target hosts

from pandas import read_csv
import sys
import os
import math
import numpy as np
import lib.tools
from subfunctions import *  

def runangle(target, delaysNnodes, lonlatNnodes, rtt):
    nblandmark = len(delaysNnodes)
    rttNnodes = rtt[goodlandRTT]
    posminrtt = goodlandRTT.index(rttNnodes)
    figure = "Figure-R/location_estimate_of_" + target + ".png"

    #boolpolygone is a variable to test whether it has a polygon or not
    boolpolygone = 0
    mindelaysNodes = min(delaysNodes)
    pos = delaysNodes.index(mindelaysNodes)
    ficslt = "Bonpoint/zonecross-" +  target

    if os.path.exists(ficslt):
        with open('ficrslt.txt') as f:
            if len(count.fields(ficslt) == 0):
                #count.fields 计算每行的字段数?
                with open('Rslt_localisation/estimate-loc.dat') as f:
                    f.write([target, lonlatNnodes[pos, 0], lonlatNnodes[pos, 1]])
                LONGITUDE = lonlatNnodes[pos, 0]
                LATITUDE = lonlatNnodes[pos, 1]
                with open('Rslt_localisation/rayon-centroid-v2.dat') as f:
                    f.write([target, lonlatNnodes[pos]])
                AIRE = math.pi*(delaysNnodes[pos]**2)
                with open('Rslt_localisation/aire-polygone-v2.dat') as f:
                    f.write([target, lonlatNnodes[pos]])
                rayon = delaysNnodes[pos]
                boolpolygone = 2
            else:
                pointcross = np.asmatrix(f.read())   
                long = len(pointcross[:, 1])
                angles = np.mat(np.zeros(long, 1))
                
                #计算起点，距离target最近的点（x0,y0），到其他所有交点的角度
                posmin = pointcross[:, 1].index()
                xzero = pointcross[posmin, 0]
                yzero = pointcross[posmin, 1]
                for j in range(0,long):
                    angles[j] = rad2deg(calculangle(xzero, yzero, pointcross[j, 1], pointcross[j,2]))
                #排序角度
                angles = sorted(angles)
                #以距离target最近的点为起点，按角度排序排列其他的点的集合pointordonne
                pointordonne = np.mat(np.zeros(2,long))
                if len(angles) >= 3:
                    compteur = 1
                    while compteur <= long:
                        for k in range(1,long):
                            ang = rad2deg(calculangle(xzero, yzero, pointcross[k, 1], pointcross[k, 2]))
                            if ang == angles:
                                pointordonne[compteur, 0] = pointcross[k, 0]
                                pointordonne[compteur, 1] = pointcross[k, 1]
                        compteur = compteur + 1
       
                    pointenKM = degre2km(pointordonne)
                    pointendegre = pointordonne
                    CxCyRAIRE = airepolygone(pointendegre,pointenKM)
                    #如果面积不为0
                    if CxCyRAIRE[3] != 0 and CxCyRAIRE[3] != 1:
                        with open('Rslt_localisation/estimate-loc.dat') as f:
                            f.write([target, formatC(CxCyRAIRE[0], digits = 2, format = "f"), formatC(CxCyRAIRE[1], digits = 2, format = "f")])

                        LONGITUDE = CxCyRAIRE[0]
                        LATITUDE = CxCyRAIRE[1]
                        with open('Rslt_localisation/rayon-centroid-v2.dat') as f:
                            f.write([target, CxCyRAIRE[2]])
                        with open('Rslt_localisation/aire-polygone-v2.dat') as f:
                            f.write([target, CxCyRAIRE[3]])
                        boolpolygone = 1

                    elif CxCyRAIRE[3] == 1:
                        aire = area.polygon(pointcross)
                        rayon = sqrt(aire/math.pi)
                        LONGITUDE = centroid.polygon(pointcross)[0]
                        LATITUDE  = centroid.polygon(pointcross)[1]
                        
                        with open('Rslt_localisation/estimate-loc.dat') as f:
                            f.write([target, LONGITUDE, LATITUDE])
                        with open('Rslt_localisation/rayon-centroid-v2.dat') as f:
                            f.write([target, rayon])   
                        with open('Rslt_localisation/aire-polygone-v2.dat') as f:
                            f.write([target, aire])

                    else:
                        for zt in range(1,long):
                            D = calculdist(deg2rad(pointcross[0, 0]), deg2rad(pointcross[0, 1]), deg2rad(pointcross[zt,0]), deg2rad(pointcross[zt,1]))
                            if D != 0:
                                rayon = D/2
                                CxCyRAIRE[2] = rayon
                                CxCyRAIRE[3] = math.pi*(rayon**2)

                        lonbar = 0
                        latbar = 0

                        for jk in range(0,long):
                            lonbar = lonbar + pointcross[jk,0]
                            latbar = latbar + pointcross[jk,1]

                        lonbar = lonbar/long
                        latbar = latbar/long
                        CxCyRAIRE[0] = lonbar
                        CxCyRAIRE[1] = latbar

                        with open('Rslt_localisation/estimate-loc.dat') as f:
                            f.write([target, formatC(CxCyRAIRE[0], digits = 2, format = "f"), formatC(CxCyRAIRE[1], digits = 2, format = "f")])
                    
                        LONGITUDE = CxCyRAIRE[0]
                        LATITUDE = CxCyRAIRE[1] 
                        rayon = CxCyRAIRE[2]
                        boolpolygone = 2

                        with open('Rslt_localisation/rayon-centroid-v2.dat') as f:
                            f.write([target, rayon])
                        with open('Rslt_localisation/aire-polygone-v2.dat') as f:
                            f.write([target, aire])      
                #只有两个交点，不够构成一个多边形	
                else:
                    #计算错误区域，我们有两个属于所有圆的两个点 
                    #calculate error area,we have two points that belong to all circles
                    if len(angles) == 2:
                        diametre = calculdist(deg2rad(pointcross[0, 0]), deg2rad(pointcross[0, 1]), deg2rad(pointcross[1, 2]), deg2rad(pointcross[1, 1]))
                        rayon = diametre/2
                        A2 = math.pi*(rayon**2)

                        with open('Rslt_localisation/aire-polygone-v2.dat') as f:
                            f.write([target, A2])
                        with open('Rslt_localisation/rayon-centroid-v2.dat') as f:
                            f.write([target, rayon])
                        lonbar = 0
                        latbar = 0

                        for j1 in range(0, long):
                            lonbar = lonbar + pointcross[j1, 0]
                            latbar = latbar + pointcross[j1, 1]

                        lonbar = lonbar/long
                        latbar = latbar/long

                        with open('Rslt_localisation/estimate-loc.dat') as f:
                            f.write([target,formatC(lonbar, digits = 2, format = "f"), formatC(latbar, digits = 2, format = "f")])
                        
                        LONGITUDE = lonbar
                        LATITUDE = latbar
                        boolpolygone = 2
                    else:
                        # we have a point belonging to all the circles
                        rayon = delaysNnodes[pos]
                        with open('Rslt_localisation/rayon-centroid-v2.dat') as f:
                            f.write([target, rayon])

                        A2 = math.pi*(rayon**2)

                        with open('Rslt_localisation/aire-polygone-v2.dat') as f:
                            f.write([target, A2])
                        with open('Rslt_localisation/estimate-loc.dat') as f:
                            f.write([target, lonlatNnodes[pos, 1], lonlatNnodes[pos, 2]])
                    
                        LONGITUDE = lonlatNnodes[pos, 0]
                        LATITUDE = lonlatNnodes[pos, 1]
                        boolpolygone = 2
    else:
        with open('Rslt_localisation/rayon-centroid-v2.dat') as f:
            f.write([target, 0])  
        with open('Rslt_localisation/aire-polygone-v2.dat') as f:
            f.write([target, 0])
        with open('Rslt_localisation/estimate-loc.dat') as f:
            f.write([target, 0])
    
        



 







                
