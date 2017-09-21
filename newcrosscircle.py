#!usr/bin/env python3
# -*- coding:utf-8 -*-

"""The function Newcrosscircle seeks the crossing points between the circles
that have as center the geographic position of the landmarks and as radius the
geographic distance from each landmark toward a given target host. For each
landmark, it draws the circle whose radius is the geographic distance between
this landmark and the target host, and centered at the position of the
landmark. This function returns a variable. If the value of this variable is
equal to zero we have no crossing points so we can not localize this target
host. If its value is different to zero, the longitude and the latitude of each
crossing point are saved in a file whose name is the IP address of the target
host and located in the folder "Crossingpoint"
################################################################################
delaysNnodes represents the distance, in km, between a given landmark and the
target. It can be seen as the radius of the circle where the center is the
geographic position (longitude,latitude), in degrees, of the landmarks. Lonlat
contains the longitude et latitude of our set of landmarks. It will be use as
the center of our circle.

"""

import numpy as np
from subfunctions import *


def newcrosscircle(target, delaysNnodes, lonlatNnodes):
    """以所有地标为圆心，以测得的半径画圆，在两两圆之间判断是否相交，并判断交集的点

    """
    nlandmarks = len(delaysNnodes)  # number of landmarks
    veclines = plotcircle(lonlatNnodes[0, 0],
                          lonlatNnodes[0, 1],
                          delaysNnodes[0])
    for i in range(1, nlandmarks):
        matline = plotcircle(lonlatNnodes[i, 0],
                             lonlatNnodes[i, 1],
                             delaysNnodes[i])
        veclines = np.c_[veclines, matline]

    L = np.size(veclines, 1)
    # lineA = np.zeros((101, 2))
    # lineB = np.zeros((101, 2))

    # 在 crossingpoint 目录下新建一个文件 crossingpoint-<ip> ，收集所有圆的交集
    nomfic = "Crossingpoint/crossingpoint" + target
    f = open(nomfic, 'w', encoding="utf-8")

    current_land, current_circle, count = 0, 0, 0   # count 表示是否有交叉点
    while current_circle <= L:
        target_land = current_land + 1
        target_circle = current_circle + 2
        while target_circle <= L and target_land <= nlandmarks:
            # 两个圆的圆心和半径
            R1 = delaysNnodes[current_land]
            lon1 = lonlatNnodes[current_land, 0]
            lat1 = lonlatNnodes[current_land, 1]
            R2 = delaysNnodes[target_land]
            lon2 = lonlatNnodes[target_land, 0]
            lat2 = lonlatNnodes[target_land, 1]

            # 两个圆上的点
            lineA = veclines[:, current_circle:(current_circle + 2)]
            lineB = veclines[:, target_circle:(target_circle + 2)]

            # 判断这两个圆是否相交
            ptsegments1 = []
            ptsegments2 = []
            if cross(lon1, lat1, R1, lon2, lat2, R2):
                # 找到第一个圆中与第二个圆相交的线段
                ptsegments1 = findsegments(lineA, R2, lon2, lat2)
                # 找到第二个圆中与第一个圆相交的线段
                ptsegments2 = findsegments(lineB, R1, lon1, lat1)
            else:
                target_land += 1
                target_circle += 2
                continue

            for first in ptsegments1:
                # 计算 lineB 与 lineA 的相交一条线段 eq1
                eq1 = equation(first[0], first[1], first[2], first[3])
                for second in ptsegments2:
                    # 找到 lineB 与 lineA 的一条线段 eq2
                    eq2 = equation(second[0], second[1], second[2], second[3])
                    # 计算两个直线的交点 ptcross
                    ptcross = pointinter(eq1, eq2,
                                         first[0], first[2],
                                         second[0], second[2])
                    # 两个线段不一定在圆内相交，因为 eq1、eq2 只是两个圆相交的线
                    # 段，但不表示这两者就是对应的，所以如果返回的结果 ptcross
                    # 非空，则说明两个确实在给定的线段范围内相交，结果就存入之前
                    # 新建的 nomfic 文件，不然就说明两个线段并没有交点，继续下一
                    # 个点继续
                    if ptcross:
                        D1 = calculdist(ptcross[0], ptcross[1], lon1, lat1)
                        D2 = calculdist(ptcross[0], ptcross[1], lon2, lat2)
                        if D1 <= R1 and D2 <= R2:
                            if (ptcross[0] >= -180 and ptcross[0] <= 180 and
                                    ptcross[1] >= -90 and ptcross[1] <= 90):
                                f.write(str(ptcross)[1: -1])
                                count += 1
            target_land += 1
            target_circle += 2
        current_land += 1
        current_circle += 2
    return count

if __name__ == '__main__':
    target = '60.247.13.2'
    delaysNnodes = np.array([13.0, 10.0, 7.0, 6.0])
    lonlatNnodes = np.array([
        [116.2068, 39.9074],
        [116.3084, 39.9305],
        [116.3826, 39.8583],
        [116.4104, 39.9050]
    ])
    newcrosscircle(target, delaysNnodes, lonlatNnodes)
