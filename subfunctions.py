import math 
import numpy as np

def ptlonlat(ln1, lt1, dist, tcin):
    d = dist/6378.16
    lon1 = math.pi*(ln1/180)
    lat1 = math.pi*(lt1/180)
    tc = math.pi*(tcin/180)
    lat = math.asin(math.sin(lat1)*math.cos(d) + math.cos(lat1)*math.sin(d)*math.cos(tc))
    dlon = math.atan2(math.sin(tc)*math.sin(d)*math.cos(lat1), math.cos(d)-math.sin(lat1)*math.sin(lat))
    lon = ((lon1-dlon+math.pi)%(2*math.pi))-math.pi
    list1 = [180*lon/math.pi, 180*lat/math.pi]
    return list1 

def plotcircle(lon, lat, r):
    angles = [ 360.0*i/100 for i in range(102) ]  
    points = map(lambda  x: ptlonlat(lon, lat, r, x), angles) 
    return points

def plotcircleV2(lon, lat, r):
    angles = [ 360.0*i/100 for i in range(102) ]  
    points = map(lambda  x: ptlonlat(lon, lat, r, x), angles) 
    #?
    #?

def cross(lon1, lat1, R1, lon2, lat2, R2):
    distC1C2 = calculdist(deg2rad(lon1), deg2rad(lat1), deg2rad(lon2), deg2rad(lat2))
    if abs(R1-R2) < distC1C2 and distC1C2 < (R1+R2):
        rep = 1
    else:
        rep = 0
    return rep

def findsegments(lineA, RayonB, lon2, lat2):
    pt = np.mat(np.zeros((3, 3)))
    extremite1 = calculdist(deg2rad(lineA[1, 1]), deg2rad(lineA[1, 2]), deg2rad(lon2), deg2rad(lat2))
    extremite2 = calculdist(deg2rad(lineA[2, 1]), deg2rad(lineA[2, 2]), deg2rad(lon2), deg2rad(lat2))
    if (extremite1 < RayonB and extremite2 > RayonB) or (extremite1 > RayonB and extremite2 < RayonB):
        pt = [np.mat(lineA[1, 1], lineA[1, 2], lineA[2, 1], lineA[2, 2])]
    for i in range(2, 100):
        extremite1= calculdist(deg2rad(lineA[i,1]), deg2rad(lineA[i,2]), deg2rad(lon2), deg2rad(lat2))
        extremite2 = calculdist(deg2rad(lineA[i+1,1]), deg2rad(lineA[i+1,2]), deg2rad(lon2), deg2rad(lat2))
        if (extremite1 < RayonB and extremite2 > RayonB) or (extremite1 > RayonB and extremite2 < RayonB):
            pointcross = mat(lineA[i, 1], lineA[i, 2], lineA[i+1, 1], lineA[i+1, 2])
            pt = [hstack((pt, pointcross))]
    return pt

def eqdroite(x1, y1, x2, y2):
    a = (y1-y2)/(x1-x2)
    b = y2-a*x2
    list2 = [a, b]
    return list2

def pointinter(eq1, eq2, lon1, lon11, lon2, lon22):
    x = (eq2[2]-eq1[2]) / (eq1[1]-eq2[1])
    y = (eq2[1]*x) + eq2[2]
    pt = ["init", "init"]
    matlon = np.mat((lon1, lon11, lon2, lon22), size=(1, 1))
    maxlon = max(matlon)
    minlon = min(matlon)
    if x >= minlon and x <= maxlon:
        pt = [x, y]
        return pt

def calculdist(lon1, lat1, lon2, lat2):
    d = 2*math.asin(math.sqrt((math.sin((lat1-lat2)/2))**2+math.cos(lat1)*math.cos(lat2)*(math.sin((lon2-lon1)/2))**2))
    d = d*6371
    return d

def deg2rad(x):
    y = (x*math.pi)/180
    return y

def rad2deg(angle):
    andegre = (180/math.pi)*angle
    return andegre
    
def degre2km(points):
    lon2km = np.mat(np.zeros(1, len(points[:1]) - 1))
    lat2km = np.mat(np.zeros(1, len(points[:1]) - 1))
    for i in range(2, len(points[:1])-1):
        D1 = calculdist(deg2rad(points[1,1]), deg2rad(points[1,2]), deg2rad(points[i,1]), deg2rad(points[i,2]))
        D2 = calculdist(deg2rad(points[1,1]+1), deg2rad(points[1,2]), deg2rad(points[i,1]), deg2rad(points[i,2]))
        lon2km[i] = abs(D2 - D1)
        D3 = calculdist(deg2rad(points[1,1]),deg2rad(points[1,2]+1),deg2rad(points[i,1]),deg2rad(points[i,2]))
        lat2km[i] = abs(D1 - D3)
    dlon = np.mean(lon2km)
    dlat = np.mean(lat2km)
    points[:1] = points[:1] * dlon
    points[:2] = points[:2] * dlat
    points = np.hstack(points, [dlon, dlat])

#Tracepolygone 函数未处理

def calculangle(lon1, lat1, lon2, lat2):
    lon1 = deg2rad(lon1)
    lat1 = deg2rad(lat1)
    lon2 = deg2rad(lon2)
    lat2 = deg2rad(lat2)
    tc1 = mod(atan2(sin(lon2-lon1)*cos(lat2),cos(lat1)*sin(lat2)-sin(lat1)*cos(lat2)*cos(lon2-lon1)), 2*pi)
    return tcl

#Tracebestline函数未处理

def fonctx(a, x, b):
    y = a*x + b
    return y

def mod(y, x):
    resultat = y - x*math.floor(y/x)
    return resultat










