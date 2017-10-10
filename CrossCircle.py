#!usr/bin/env python3
# -*- coding:utf-8 -*-

"""CrossCircle.py 计算以地标为圆心，画圆得到的交点的坐标集合，交给 GoodPoint.py
处理

The function CrossCircle seeks the crossing points between the circles that
have as center the geographic position of the landmarks and as radius the
geographic distance from each landmark toward a given target host. For each
landmark, it draws the circle whose radius is the geographic distance between
this landmark and the target host, and centered at the position of the
landmark. This function returns a variable. If the value of this variable is
equal to zero we have no crossing points so we can not localize this target
host. If its value is different to zero, the longitude and the latitude of each
crossing point are saved in a file whose name is the IP address of the target
host and located in the folder "Crossingpoint"
###############################################################################
delays represents the distance, in km, between a given landmark and the
target. It can be seen as the radius of the circle where the center is the
geographic position (longitude, latitude), in degrees, of the landmarks. Lonlat
contains the longitude and latitude of our set of landmarks. It will be use as
the center of our circle.

"""

import numpy as np
from Utils import plotcircle, findsegments, pointinter, \
                  calcdist, cross, equation


def crosscircle(target, delays, lonlats):
    """以所有地标为圆心，以测得的半径画圆，在两两圆之间判断是否相交，并判断交集的点

    """
    nlandmarks = len(delays)  # 地标个数
    veclines = plotcircle(lonlats[0, 0], lonlats[0, 1], delays[0])
    for i in range(1, nlandmarks):
        matline = plotcircle(lonlats[i, 0], lonlats[i, 1], delays[i])
        veclines = np.c_[veclines, matline]

    cols = np.size(veclines, 1)
    # line_a = np.zeros((101, 2))
    # line_b = np.zeros((101, 2))

    # 在 crossingpoint 目录下新建一个文件 crossingpoint-<ip> ，收集所有圆的交集
    crossingpoint = "Crossingpoint/crossingpoint-"
    nomfic = "".join([crossingpoint, target])
    cp = open(nomfic, 'w', encoding="utf-8")

    current_land, current_circle, count = 0, 0, 0   # count 表示是否有交叉点
    while current_circle < cols:
        target_land = current_land + 1
        target_circle = current_circle + 2
        while target_circle < cols and target_land < nlandmarks:
            # 两个圆的圆心和半径
            r1 = delays[current_land]
            lon1 = lonlats[current_land, 0]
            lat1 = lonlats[current_land, 1]
            r2 = delays[target_land]
            lon2 = lonlats[target_land, 0]
            lat2 = lonlats[target_land, 1]

            # 两个圆上的点
            line_a = veclines[:, current_circle:(current_circle + 2)]
            line_b = veclines[:, target_circle:(target_circle + 2)]

            # 判断这两个圆是否相交
            ptsegments1 = []
            ptsegments2 = []
            if cross(lon1, lat1, r1, lon2, lat2, r2):
                # 找到第一个圆中与第二个圆相交的线段
                ptsegments1 = findsegments(line_a, r2, lon2, lat2)
                # 找到第二个圆中与第一个圆相交的线段
                ptsegments2 = findsegments(line_b, r1, lon1, lat1)
            else:
                target_land += 1
                target_circle += 2
                continue

            for first in ptsegments1:
                # 计算 line_b 与 line_a 的相交一条线段 eq1
                eq1 = equation(first[0], first[1], first[2], first[3])
                for second in ptsegments2:
                    # 找到 line_b 与 line_a 的一条线段 eq2
                    eq2 = equation(second[0], second[1], second[2], second[3])
                    # 计算两个直线的交点 ptcross
                    ptcross = pointinter(eq1, eq2,
                                         first[0], first[2],
                                         second[0], second[2])
                    # 两个线段不一定在圆内相交，因为 eq1、eq2 只是两个圆相交的
                    # 线段，但不表示这两者就是对应的，所以如果返回的结果
                    # ptcross 非空，则说明两个确实在给定的线段范围内相交，结果
                    # 就存入之前新建的 nomfic 文件，不然就说明两个线段并没有交
                    # 点，继续下一个点继续
                    if ptcross:
                        d1 = calcdist(ptcross[0], ptcross[1], lon1, lat1)
                        d2 = calcdist(ptcross[0], ptcross[1], lon2, lat2)
                        if d1 <= r1 and d2 <= r2:
                            if ((ptcross[0] >= -180 and ptcross[0] <= 180 and
                                 ptcross[1] >= -90 and ptcross[1] <= 90)):
                                cp.write("".join([str(ptcross)[1:-1], "\n"]))
                                count += 1
            target_land += 1
            target_circle += 2
        current_land += 1
        current_circle += 2
    cp.close()
    return count

if __name__ == '__main__':
    TEST_TARGET = '60.247.13.2'
    TEST_DELAYS = np.array([13.0, 10.0, 7.0, 6.0])
    TEST_LONLATNODES = np.array([
        [116.2068, 39.9074],
        [116.3084, 39.9305],
        [116.3826, 39.8583],
        [116.4104, 39.9050]
    ])
    crosscircle(TEST_TARGET, TEST_DELAYS, TEST_LONLATNODES)
