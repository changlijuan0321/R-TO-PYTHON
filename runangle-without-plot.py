from pandas import read_csv
import sys
import os
import math
import numpy as np
import lib.tools
from subfunctions import *  

def runangle(target, delayNnodes, lonlatNnodes, rtt):
    nblandmark = len(delayNnodes)
    rttNnodes = rtt[goodlandRTT]
    posminrtt = goodlandRTT.index(rttNnodes)
    figure = str.join("Figure-R/location_estimate_of_",target,".png")

    #boolpolygone is a variable to test whether it has a polygon or not
    boolpolygone = 0
    mindelaysNodes = min(delayNodes)
    pos = delaynNodes.index(mindelaysNodes)
    ficslt = str.join("Bonpoint/zonecross-", target)

    if os.path.exists(ficslt):
        with open('ficrslt.txt') as f:
            if len(count.fields(ficslt) == 0):
                #count.fields 计算每行的字段数?
                with open('Rslt_localisation/estimate-loc.dat') as f:
                    f.write([target, lonlatNnodes[pos, 1], lonlatNnodes[pos, 2]])
                LONGITUDE = lonlatNnodes[pos,1]
                LATITUDE = lonlatNnodes[pos,2]
                with open('Rslt_localisation/rayon-centroid-v2.dat') as f:
                    f.write([target, lonlatNnodes[pos]])
                AIRE = math.pi*(delaysNnodes[pos]**2)
                with open('Rslt_localisation/aire-polygone-v2.dat') as f:
                    f.write([target, lonlatNnodes[pos]])
                rayon = delayNnodes[pos]
                boolpolygone = 2
            else:
                pointcross = np.asmatrix(f.read())   
                long = len(pointcross[:, 1])
                angles = np.mat(np.zeros(long,1))
                
                #计算起点，距离target最近的点（x0,y0），到其他所有交点的角度
                posmin = pointcross[:,1].index()
                xzero = pointcross[posmin,1]
                yzero = pointcross[posmin,2]
                for j in range(1,long):
                    angles[j] = rad2deg(calculangle(xzero,yzero,pointcross[j,1],pointcross[j,2]))
                #排序角度
                angles = sorted(angles)
                #以距离target最近的点为起点，按角度排序排列其他的点的集合pointordonne
                pointordonne = np.mat(np.zeros(2,long))
                if len(angles) >= 3:
                    compteur = 1
                    while compteur <= long:
                        for k in range(1,long):
                            ang = rad2deg(calculangle(xzero,yzero,pointcross[k,1],pointcross[k,2]))
                            if ang == angles:
                                #处理可能有问题、、R语言的identical是要完全一样
                                pointordonne[compteur, 1] = pointcross[k, 1]
                                pointordonne[compteur, 2] = pointcross[k, 2]
                        compteur = compteur + 1

                    pointenKM = degre2km(pointordonne)
                    pointendegre = pointordonne
                    CxCyRAIRE = airepolygone(pointendegre,pointenKM)
                    #如果面积不为0
                    if CxCyRAIRE[4] != 0 and CxCyRAIRE[4] != 1:
                        with open('Rslt_localisation/estimate-loc.dat') as f:
                            f.write([target, formatC(CxCyRAIRE[1], digits = 2, format = "f"), formatC(CxCyRAIRE[2], digits = 2, format = "f")])

                        LONGITUDE = CxCyRAIRE[1]
                        LATITUDE = CxCyRAIRE[2]
                        with open('Rslt_localisation/rayon-centroid-v2.dat') as f:
                            f.write([target, CxCyRAIRE[3]])
                        with open('Rslt_localisation/aire-polygone-v2.dat') as f:
                            f.write([target, CxCyRAIRE[4]])
                        boolpolygone = 1

                    elif CxCyRAIRE[4] == 1:
                        aire = area.polygon(pointcross)
                        rayon = sqrt(aire/math.pi)
                        LONGITUDE = centroid.polygon(pointcross)[1]
                        LATITUDE  = centroid.polygon(pointcross)[2]
                        
                        with open('Rslt_localisation/estimate-loc.dat') as f:
                            f.write([target, LONGITUDE, LATITUDE])
                        with open('Rslt_localisation/rayon-centroid-v2.dat') as f:
                            f.write([target, rayon])   
                        with open('Rslt_localisation/aire-polygone-v2.dat') as f:
                            f.write([target, aire])

                    else:
                        for zt in range(2,long):
                            D = calculdist(deg2rad(pointcross[1,1]), deg2rad(pointcross[1,2]), deg2rad(pointcross[zt,1]), deg2rad(pointcross[zt,2]))
                            if D != 0:
                                rayon = D/2
                                CxCyRAIRE[3] = rayon
                                CxCyRAIRE[4] = math.pi*(rayon**2)

                        lonbar = 0
                        latbar = 0

                        for jk in range(1,long):
                            onbar = lonbar + pointcross[jk,1]
                            latbar = latbar + pointcross[jk,2]

                        lonbar = lonbar/long
                        latbar = latbar/long
                        CxCyRAIRE[1] = lonbar
                        CxCyRAIRE[2] = latbar

                        with open('Rslt_localisation/estimate-loc.dat') as f:
                            f.write([target, formatC(CxCyRAIRE[1], digits = 2, format = "f"), formatC(CxCyRAIRE[2], digits = 2, format = "f")])
                    
                        LONGITUDE = CxCyRAIRE[1]
                        LATITUDE = CxCyRAIRE[2] 
                        rayon = CxCyRAIRE[3]
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
                        diametre = calculdist(deg2rad(pointcross[1,1]), deg2rad(pointcross[1,2]), deg2rad(pointcross[2,1]), deg2rad(pointcross[2,2]))
                        rayon = diametre/2
                        A2 = math.pi*(rayon**2)

                        with open('Rslt_localisation/aire-polygone-v2.dat') as f:
                            f.write([target, A2])
                        with open('Rslt_localisation/rayon-centroid-v2.dat') as f:
                            f.write([target, rayon])
                        lonbar = 0
                        latbar = 0

                        for j1 in range(1,long):
                            lonbar = lonbar + pointcross[j1, 1]
                            latbar = latbar + pointcross[j1, 2]

                        lonbar = lonbar/long
                        latbar = latbar/long

                        with open('Rslt_localisation/estimate-loc.dat') as f:
                            f.write([target,formatC(lonbar, digits = 2, format = "f"), formatC(latbar, digits = 2, format = "f")])
                        
                        LONGITUDE = lonbar
                        LATITUDE = latbar
                        boolpolygone = 2
                    else:
                        # we have a point belonging to all the circles
                        rayon = delayNnodes[pos]
                        with open('Rslt_localisation/rayon-centroid-v2.dat') as f:
                            f.write([target, rayon])

                        A2 = math.pi*(rayon**2)

                        with open('Rslt_localisation/aire-polygone-v2.dat') as f:
                            f.write([target, A2])
                        with open('Rslt_localisation/estimate-loc.dat') as f:
                            f.write([target, lonlatNnodes[pos, 1], lonlatNnodes[pos, 2]])
                    
                        LONGITUDE = lonlatNnodes[pos, 1]
                        LATITUDE = lonlatNnodes[pos, 2]
                        boolpolygone = 2
    else:
        with open('Rslt_localisation/rayon-centroid-v2.dat') as f:
            f.write([target, 0])  
        with open('Rslt_localisation/aire-polygone-v2.dat') as f:
            f.write([target, 0])
        with open('Rslt_localisation/estimate-loc.dat') as f:
            f.write([target, 0])
    
        



 







                
