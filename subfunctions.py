#!usr/bin/env python3
# -*- coding:utf-8 -*-

"""Util Functions

"""

import math
import numpy as np
# import airepolygone
# import matplotlib.pyplot as plt


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
    """判断圆心和半径分别为 ((lon1, lat1), R1) 和 ((lon2, lat2), R2) 的两个圆是否相
    交，若相交则返回 True

    """
    distC1C2 = calculdist(lon1, lat1, lon2, lat2)
    # 排除了相切的情况
    if abs(R1 - R2) < distC1C2 and distC1C2 < (R1 + R2):
        return True
    else:
        return False


def findsegments(lineA, radiusB, lon2, lat2):
    """找到 lineA 与圆 ((lon2, lat2), radiusB) 之间的相交线段

    """
    points = np.zeros((1, 4))
    # 判断线段 point1->point2 是否与圆 B 相交
    for i in range(np.size(lineA, 0) - 1):
        distance1 = calculdist(lineA[i, 0], lineA[i, 1], lon2, lat2)
        distance2 = calculdist(lineA[i + 1, 0], lineA[i + 1, 1], lon2, lat2)
        if ((distance1 < radiusB and distance2 > radiusB) or
                (distance1 > radiusB and distance2 < radiusB)):
            tmp = np.array([lineA[i, 0], lineA[i, 1],
                            lineA[i + 1, 0], lineA[i + 1, 1]])
            points = np.vstack((points, tmp))
    return points[1:]


def equation(x1, y1, x2, y2):
    a = (y1 - y2) / (x1 - x2)
    b = y2 - a * x2
    return [a, b]


def pointinter(eq1, eq2, lon1, lon11, lon2, lon22):
    """计算两个直线方程 eq1 与 eq2 在纬度范围为 [lon1, lon11] 之间的交点

    若存在交点则返回交点坐标，否则返回空

    """
    x = - (eq2[1] - eq1[1]) / (eq2[0] - eq1[0])
    y = (eq2[0] * x) + eq2[1]
    pt = []
    matlon = [lon1, lon11, lon2, lon22]
    maxlon = max(matlon)
    minlon = min(matlon)
    if x >= minlon and x <= maxlon:
        pt = [x, y]
    return pt


def calculdist(lon1, lat1, lon2, lat2, R=6371):
    """计算经纬度分别为 (lon1, lat1) 和 (lon2, lat2) 的两个点之间的距离

    输入单位为 degree，输出单位为 km

    Haversine formula:
        a = hav(lat1 - lat2) + cos(lat1) * cos(lat2) * hav(lon1 - lon2)

    hav() 的定义为:
        hav(x) = sin(x / 2) ** 2

    所以角距离 c 表示为:
        c = 2 * arctan(sqrt(a) / sqrt(1 - a))
        或者
        c = 2 * arcsin(min(1, sqrt(a)))

    所以两个坐标之间的距离:
        d = R * c

    For more details: http://www.movable-type.co.uk/scripts/latlong.html

    """
    lon1, lat1, lon2, lat2 = map(deg2rad, (lon1, lat1, lon2, lat2))
    havlat = np.sin((lat1 - lat2) / 2)
    havlon = np.sin((lon1 - lon2) / 2)
    ad = 2 * np.arcsin(
        np.sqrt(
            (havlat ** 2) + np.cos(lat1) * np.cos(lat2) * (havlon ** 2)))
    d = ad * R
    return d


def deg2rad(x):
    """单位转换: degree -> radian

    """
    return x * np.pi / 180


def rad2deg(angle): 
    """单位转换: radian -> degree

    """
    return 180 / np.pi * angle


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
