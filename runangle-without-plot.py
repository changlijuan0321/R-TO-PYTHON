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
#import lib.tool
from subfunctions import * 
from io import StringIO 
import scipy

def runangle(target, delaysNnodes, lonlatNnodes, rtt):
    nblandmark = len(delaysNnodes)
    # rtt contains the min RTT between target host and landmarks
    #TODO
    goodlandRTT = [1, nblandmark]
    rttNnodes = [goodlandRTT]
    # posminrtt = goodlandRTT.index(rttNnodes)
    figure = "Figure-R/location_estimate_of_" + target + ".png"
    # we output the result in png format, and create the file .png
    # boolpolygone is a variable to test whether it has a polygon or not
    boolpolygone = 0
    mindelaysNodes = min(delaysNodes)
    pos = delaysNodes.index(mindelaysNodes)
    ficrslt = "Bonpoint/zonecross-" +  target
    # we test if the file zonecross exits which expressed the targeted host is localizable
    if os.path.exists(ficrslt):
        with open(ficrslt, 'r') as f:
            # we test if we have no crossing point
            if boolpolygone != 0:
                # TODO 
                with open('Rslt_localisation/estimate-loc.dat', 'w') as f:
                    f.write(str([target, lonlatNnodes[pos, 0], lonlatNnodes[pos, 1]]))
                LONGITUDE = lonlatNnodes[pos, 0]
                LATITUDE = lonlatNnodes[pos, 1]
                with open('Rslt_localisation/rayon-centroid-v2.dat', 'w') as f:
                    f.write(str([target, lonlatNnodes[pos]]))
                # calculation surface of the circle
                AIRE = math.pi*(delaysNnodes[pos]**2)
                with open('Rslt_localisation/aire-polygone-v2.dat', 'w') as f:
                    f.write(str([target, lonlatNnodes[pos]]))
                rayon = delaysNnodes[pos]
                boolpolygone = 2
            # we have crossing points
            else:
                c = StringIO(f.read())     
                pointcross = np.loadtxt(c, delimiter=',')
                long = len(pointcross[:, 1])
                angles = np.zeros((long, 1))
                
                #计算起点，距离target最近的点（x0,y0），到其他所有交点的角度
                #posmin = pointcross.index(min(pointcross[:, 0]))
                #TODO
                posmin = 1
                xzero = pointcross[posmin, 0]
                yzero = pointcross[posmin, 1]
                for j in range(0, long):
                    angles[j] = rad2deg(calculangle(xzero, yzero, pointcross[j, 0], pointcross[j,1]))
                #排序角度
                angles = sorted(angles)
                #以距离target最近的点为起点，按角度排序排列其他的点的集合pointordonne
                pointordonne = np.zeros((2, long))
                if len(angles) >= 3:
                    compteur = 1
                    while compteur <= long:
                        for k in range(1, long):
                            ang = rad2deg(calculangle(xzero, yzero, pointcross[k, 0], pointcross[k, 1]))
                            if ang == angles:
                                pointordonne[compteur, 0] = pointcross[k, 0]
                                pointordonne[compteur, 1] = pointcross[k, 1]
                        compteur = compteur + 1
                    # calculate the area and center of the polygon(Cx,Cy), R=sqrt(Aire/pi)
                    # convert degrees to km (lon,lat)->(dlon*lon,dlat*lat)
                    # dlon and dlat represent the conversion coefficients for each polygon
                    pointenKM = degre2km(pointordonne)
                    pointendegre = pointordonne
                    # calculate of the centroid of the polygon and the radius of the circle of this polygon, airepolygone
                    # is the function that returns the centroid, radius and area of the polygon, and traces the vertices of the polygon
                    CxCyRAIRE = airepolygone(pointendegre,pointenKM)

                    #如果面积不为0
                    if CxCyRAIRE[3] != 0 and CxCyRAIRE[3] != 1:
                        # we trace the centroid of the polygon (estimate of the location of the point)
                        with open('Rslt_localisation/estimate-loc.dat', 'w') as f:
                            f.write(str([target, roun(CxCyRAIRE[0], 2), round(CxCyRAIRE[1], 2)]))

                        LONGITUDE = CxCyRAIRE[0]
                        LATITUDE = CxCyRAIRE[1]
                        with open('Rslt_localisation/rayon-centroid-v2.dat', 'w') as f:
                            f.write(str([target, CxCyRAIRE[2]]))
                        with open('Rslt_localisation/aire-polygone-v2.dat', 'w') as f:
                            f.write(str([target, CxCyRAIRE[3]]))
                        # we have a polygon we put the variable a 1
                        boolpolygone = 1

                    elif CxCyRAIRE[3] == 1:
                        aire = area.polygon(pointcross)
                        rayon = sqrt(aire/math.pi)
                        LONGITUDE = scipy.cluster.hierarchy.centroid(pointcross[0])
                        LATITUDE  = scipy.cluster.hierarchy.centroid(pointcross[1])
                        #TODO 
                        with open('Rslt_localisation/estimate-loc.dat', 'w') as f:
                            f.write(str([target, LONGITUDE, LATITUDE]))
                        with open('Rslt_localisation/rayon-centroid-v2.dat', 'w') as f:
                            f.write(str([target, rayon]))   
                        with open('Rslt_localisation/aire-polygone-v2.dat', 'w') as f:
                            f.write(str([target, aire]))
                    # if AIRE == 0 we have identical points in our matrix and in fact we have only 2 points or 1 point we
                    # look for the barycentre
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

                        with open('Rslt_localisation/estimate-loc.dat', 'w') as f:
                            f.write(str([target, formatC(CxCyRAIRE[0], digits = 2, format = "f"), formatC(CxCyRAIRE[1], digits = 2, format = "f")]))
                    
                        LONGITUDE = CxCyRAIRE[0]
                        LATITUDE = CxCyRAIRE[1] 
                        rayon = CxCyRAIRE[2]
                        boolpolygone = 2

                        with open('Rslt_localisation/rayon-centroid-v2.dat', 'w') as f:
                            f.write(str([target, rayon]))
                        with open('Rslt_localisation/aire-polygone-v2.dat', 'w') as f:
                            f.write(str([target, aire]))      
                #只有两个交点，不够构成一个多边形	
                else:
                    #计算错误区域，我们有两个属于所有圆的两个点 
                    #calculate error area,we have two points that belong to all circles
                    if len(angles) == 2:
                        diametre = calculdist(deg2rad(pointcross[0, 0]), deg2rad(pointcross[0, 1]), deg2rad(pointcross[1, 0]), deg2rad(pointcross[1, 1]))
                        rayon = diametre/2
                        A2 = math.pi*(rayon**2)

                        with open('Rslt_localisation/aire-polygone-v2.dat', 'w') as f:
                            f.write(str([target, A2]))
                        with open('Rslt_localisation/rayon-centroid-v2.dat', 'w') as f:
                            f.write(str([target, rayon]))
                        lonbar = 0
                        latbar = 0

                        for j1 in range(0, long):
                            lonbar = lonbar + pointcross[j1, 0]
                            latbar = latbar + pointcross[j1, 1]

                        lonbar = lonbar/long
                        latbar = latbar/long

                        with open('Rslt_localisation/estimate-loc.dat', 'w') as f:
                            f.write(str([target, round(lonbar,  2), round(latbar,  2)]))
                        LONGITUDE = lonbar
                        LATITUDE = latbar
                        boolpolygone = 2
                    else:
                        # we have a point belonging to all the circles
                        rayon = delaysNnodes[pos]
                        with open('Rslt_localisation/rayon-centroid-v2.dat', 'w') as f:
                            f.write(str([target, rayon]))

                        A2 = math.pi*(rayon**2)

                        with open('Rslt_localisation/aire-polygone-v2.dat', 'w') as f:
                            f.write(str([target, A2]))
                        with open('Rslt_localisation/estimate-loc.dat', 'w') as f:
                            f.write(str([target, lonlatNnodes[pos, 0], lonlatNnodes[pos, 1]]))
                    
                        LONGITUDE = lonlatNnodes[pos, 0]
                        LATITUDE = lonlatNnodes[pos, 1]
                        boolpolygone = 2
    # ficslt 文件不存在，即无法定位
    else:
        with open('Rslt_localisation/rayon-centroid-v2.dat', 'w') as f:
            f.write(str([target, 0]))  
        with open('Rslt_localisation/aire-polygone-v2.dat', 'w') as f:
            f.write(str([target, 0]))
        with open('Rslt_localisation/estimate-loc.dat', 'w') as f:
            f.write(str([target, 0]))



if __name__ == '__main__':
    target = "123.123.123.123"
    delaysNodes = [10, 20, 30]
    lonlatNnodes = np.array([
        [30.01, 50.21],
        [40.01, 60.21],
        [10.01, 20.21]
    ])
    rtt = 12
    runangle(target, delaysNodes, lonlatNnodes, rtt)
        



 







                
