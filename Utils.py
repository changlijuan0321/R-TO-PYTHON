#!usr/bin/env python3
# -*- coding:utf-8 -*-

"""Util Functions

"""

import numpy as np


def ptlonlat(ln1, lt1, dist, tcin):
    """以 (ln1, lt1) 为起点，dist 为距离，tcin 为发射角的对应的线段终点坐标

    """
    # dist 转换为相对地球半径的弧度
    d = dist / 6378.16
    lon1, lat1 = map(deg2rad, [ln1, lt1])
    tc = np.pi * (tcin / 180)
    dlon = np.arctan2(np.sin(tc) * np.sin(d) * np.cos(lat1),
                      np.cos(d) - np.sin(lat1) * np.sin(lat1))

    lat = np.arcsin(np.sin(lat1) * np.cos(d) +
                    np.cos(lat1) * np.sin(d) * np.cos(tc))
    # 限定经度范围为 [-π, π]
    lon = ((lon1 - dlon + np.pi) % (2 * np.pi)) - np.pi
    return [rad2deg(lon), rad2deg(lat)]


def plotcircle(lon, lat, radius):
    """以 (lon, lat) 为圆心，以 radius 为半径画圆
    具体方法：
        1. 在圆上取 101 个点，首尾闭合
        2. 使用 plotlonlat() 计算每个点的经纬度坐标
        3. 返回二维数组，每一行为一个点的经纬度坐标
    """
    # this function gives the coordinates of points which form a given circle
    angles = [360.0 * i / 100 for i in range(102)]
    points = map(lambda x: ptlonlat(lon, lat, radius, x), angles)
    return np.array(list(points))


def cross(lon1, lat1, r1, lon2, lat2, r2):
    """判断圆心和半径分别为 ((lon1, lat1), r1) 和 ((lon2, lat2), r2) 的两个圆是否相
    交，若相交则返回 True

    """
    distc1c2 = calcdist(lon1, lat1, lon2, lat2)
    # 排除了相切的情况
    if abs(r1 - r2) < distc1c2 and distc1c2 < (r1 + r2):
        return True
    else:
        return False


def findsegments(line_a, radiusB, lon2, lat2):
    """找到 line_a 与圆 ((lon2, lat2), radiusB) 之间的相交线段

    """
    points = np.zeros((1, 4))
    # 判断线段 point1->point2 是否与圆 B 相交
    for i in range(np.size(line_a, 0) - 1):
        distance1 = calcdist(line_a[i, 0], line_a[i, 1], lon2, lat2)
        distance2 = calcdist(line_a[i + 1, 0], line_a[i + 1, 1], lon2, lat2)
        if (((distance1 < radiusB and distance2 > radiusB) or
             (distance1 > radiusB and distance2 < radiusB))):
            tmp = [line_a[i, 0], line_a[i, 1], line_a[i + 1, 0], line_a[i + 1, 1]]
            points = np.vstack((points, tmp))
    return points[1:]


def equation(x1, y1, x2, y2):
    """计算以 (x1, y1) 和 (x2, y2) 为端点的线段所在直线的方程
    """
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
    maxlon, minlon = max(matlon), min(matlon)
    if x >= minlon and x <= maxlon:
        pt = [x, y]
    return pt


def calcdist(lon1, lat1, lon2, lat2, R=6371):
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
    ad = 2 * np.arcsin(np.sqrt(
            (havlat ** 2) + np.cos(lat1) * np.cos(lat2) * (havlon ** 2)))
    d = ad * R
    return d


def deg2rad(deg):
    """单位转换: degree -> radian

    """
    return np.pi * deg / 180


def rad2deg(rad):
    """单位转换: radian -> degree

    """
    return 180 * rad / np.pi


def degre2km(points):
    lon2km = np.zeros(1)
    lat2km = np.zeros(1)
    xzero, yzero = points[0]
    for point in points[1:]:
        d1 = calcdist(xzero, yzero, point[0], point[1])
        d2 = calcdist(xzero + 1, yzero, point[0], point[1])
        d3 = calcdist(xzero, yzero + 1, point[0], point[1])
        lon2km = np.append(lon2km, abs(d1 - d2))
        lat2km = np.append(lat2km, abs(d1 - d3))
    dlon = np.mean(lon2km)
    dlat = np.mean(lat2km)
    points = points * [dlon, dlat]
    return np.vstack((points, [dlon, dlat]))


def calcangle(lon1, lat1, lon2, lat2):
    lon1, lat1, lon2, lat2 = map(deg2rad, (lon1, lat1, lon2, lat2))
    tc1 = np.mod(
        np.arctan2(np.sin(lon2 - lon1) * np.cos(lat2),
                   (np.cos(lat1) * np.sin(lat2) -
                    np.sin(lat1) * np.cos(lat2) * np.cos(lon2 - lon1))),
        2 * np.pi
    )
    return rad2deg(tc1)


# Tracepolygone <- function(mat){
#       long<-length(mat[,1])
#       nbelement<-length(chull(mat))
#       if( nbelement==long ){
#         ## on trace les sommets du polygone
#               for( i in 1:long ){
#                       points(mat[i,1],mat[i,2],pch=20,cex=0.5)
#               }
#         ## on joint les sommet du polygone
#               for( j in 1:(long-1) ){
#                       lines(c(mat[j,1], mat[(j+1),1]), c(mat[j,2], mat[(j+1), 2]),
#                               col="blue")
#               }
#         ## joindre le dernier somment du polygone au premier point
#               lines(c(mat[long,1], mat[1,1]), c(mat[long,2], mat[1,2]), col="blue")
#       }
#       title(sub="The area of the tagged polygon represents the confidence region",
#               col.sub="blue", cex.sub=1)
# }


# Tracebestline()
# Tracebestline <- function (z,dists,rawdelays) {
#   plot(dists[z,], rawdelays[z,], type="n",
#        xlab=paste("distance (km)"), ylab=paste("delay (ms)"))
#   points(dists[z,][-z], rawdelays[z,][-z])
#     ## text(dists[z,][-z],rawdelays[z,][-z],cex=0.7,pos=1)
#   x <- fonctx(bls[z,2], 15000, bls[z,1])
#   lines(c(0,15000), c(bls[z,1],x), col="red")
# }


# fonctx
# fonctx <- function (a,x,b)
# {
#   y <- a*x+b
#   y
# }
