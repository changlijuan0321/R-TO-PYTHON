import math 
import numpy as np
import airepolygone
import matplotlib.pyplot as plt

def ptlonlat(ln1, lt1, dist, tcin):
    # convert radius of the earth from from km to radians
    d = dist / 6378.16
    # convert angles from degrees to radians
    lon1 = np.pi * (ln1 / 180)
    lat1 = np.pi * (lt1 / 180)
    tc = np.pi * (tcin / 180)
    lat = np.arcsin(np.sin(lat1) * np.cos(d) +
                    np.cos(lat1) * np.sin(d) * np.cos(tc))
    dlon = np.arctan2(np.sin(tc) * np.sin(d) * np.cos(lat1),
                      np.cos(d) - np.sin(lat1) * np.sin(lat))
    lon = ((lon1 - dlon + np.pi) % (2 * np.pi)) - np.pi
    return [180 * lon / np.pi, 180 * lat / np.pi]

def plotcircle(lon, lat, r):
    # this function gives the coordinates of points which form a given circle
    angles = [360.0 * i / 100 for i in range(102)]
    points = map(lambda x: ptlonlat(lon, lat, r, x), angles) 
    return np.array(list(points))

def cross(lon1, lat1, R1, lon2, lat2, R2):
    distC1C2 = calculdist(lon1, lat1, lon2, lat2)
    if abs(R1 - R2) < distC1C2 and distC1C2 < (R1 + R2):
        return True
    else:
        return False

def findsegments(lineA, radiusB, lon2, lat2):
    pt = np.zeros((3, 3))
    point1 = calculdist(lineA[0, 0], lineA[0, 1], lon2, lat2)
    point2 = calculdist(lineA[1, 0], lineA[1, 1], lon2, lat2)
    if ((point1 < radiusB and point2 > radiusB) or
        (point1 > radiusB and point2 < radiusB)):
        pt = [lineA[0, 0], lineA[0, 1], lineA[1, 0], lineA[1, 1]]
    for i in range(1, 101):
        point1 = calculdist(lineA[i, 0], lineA[i, 1], lon2, lat2)
        point2 = calculdist(lineA[i + 1, 0], lineA[i + 1, 1], lon2, lat2)
        if ((point1 < radiusB and point2 > radiusB) or
            (point1 > radiusB and point2 < radiusB)):
            pointcross = [lineA[i, 0], lineA[i, 1],
                          lineA[i + 1, 0], lineA[i + 1, 1]]
            pt = np.c_[pt, pointcross]
    return pt

def eqdroite(x1, y1, x2, y2):
    a = (y1 - y2) / (x1 - x2)
    b = y2 - a * x2
    return [a, b]

def pointinter(eq1, eq2, lon1, lon11, lon2, lon22):
    x = (eq2[1] - eq1[1]) / (eq1[0] - eq2[0])
    y = (eq2[1] * x) + eq2[1]
    pt = ["init", "init"]
    matlon = [lon1, lon11, lon2, lon22]
    maxlon = max(matlon)
    minlon = min(matlon)
    if x >= minlon and x <= maxlon:
        pt = [x, y]
    return pt

def calculdist(lon1, lat1, lon2, lat2):
    lon1, lat1, lon2, lat2 = map(deg2rad, (lon1, lat1, lon2, lat2))
    d = 2 * np.arcsin(np.sqrt((np.sin((lat1 - lat2) / 2)) ** 2 + np.cos(lat1) * np.cos(lat2) * (np.sin((lon2 - lon1) / 2)) ** 2))
    d = d * 6371
    return d

def deg2rad(x): return x * np.pi / 180

def rad2deg(angle): return 180 / np.pi * angle

def degre2km(points):
    lon2km = np.mat(np.zeros(1, len(points[:1]) - 1))
    lat2km = np.mat(np.zeros(1, len(points[:1]) - 1))
    for i in range(2, len(points[:1])-1):
        D1 = calculdist(deg2rad(points[0, 0]), deg2rad(points[0, 1]), deg2rad(points[i, 0]), deg2rad(points[i, 1]))
        D2 = calculdist(deg2rad(points[0, 0]+1), deg2rad(points[0, 1]), deg2rad(points[i, 0]), deg2rad(points[i, 1]))
        lon2km[i] = abs(D2 - D1)
        D3 = calculdist(deg2rad(points[0, 0]), deg2rad(points[0, 1]+1), deg2rad(points[i,0]), deg2rad(points[i,1]))
        lat2km[i] = abs(D1 - D3)
    dlon = np.mean(lon2km)
    dlat = np.mean(lat2km)
    points[:1] = points[:1] * dlon
    points[:2] = points[:2] * dlat
    points2 = [points, [dlon, dlat]]
    return points2

def calculangle(lon1, lat1, lon2, lat2):
    lon1 = deg2rad(lon1)
    lat1 = deg2rad(lat1)
    lon2 = deg2rad(lon2)
    lat2 = deg2rad(lat2) 
    tc1 = mod(math.atan2(math.sin(lon2-lon1) * math.cos(lat2), math.cos(lat1)*math.sin(lat2) - math.sin(lat1) * math.cos(lat2) * math.cos(lon2-lon1)), 2*math.pi)
    return tc1

def fonctx(a, x, b):
    y = a*x + b
    return y

def mod(y, x):
    resultat = y - x*math.floor(y/x)
    return resultat
