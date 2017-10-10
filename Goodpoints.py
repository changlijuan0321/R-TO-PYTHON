#!usr/bin/env python3
# -*- coding:utf-8 -*-

"""GoodPoints.py 负责筛选出现在所有圆中的点，即得到目标区域的轮廓

This function seeks the "good points" among the crossing points found
before. These good points represent the intersection region of all
circles. Good points are any points that are located in all circles. These good
points allow us to find the area where we expect that the target host is
located. These good points will be the vertices of the polygon. Note that the
polygon is the approximation of this intersection region. The longitude and the
latitude of good points are saved in a file whose name is the IP address of the
given target host and the file is located in the folder "Goodpoint".

It should be noted that if we have two good points, the intersection region
will be the circle whose diameter is the geographic distance of these two
points.

"""

import os
import numpy as np
from Utils import *


def goodpoints(target, delays, lonlats):
    nlandmark = len(delays)

    goodpoint = "Goodpoint/zonecross-"
    goodpoint_ip = "".join([goodpoint, target])
    if os.path.exists(goodpoint_ip):
        os.remove(goodpoint_ip)
    gp = open(goodpoint_ip, 'a')

    pos = np.argmin(delays)     # 最小半径
    bondist = delays[pos]       # 对应的索引

    crossingpoint = "Crossingpoint/crossingpoint-"
    crossingpoint_ip = "".join([crossingpoint, target])
    crosspts = np.genfromtxt(crossingpoint_ip, delimiter=',')

    test = 0
    for point in crosspts:
        count = 0
        D1 = calcdist(lonlats[pos, 0], lonlats[pos, 1], point[0], point[1])
        if D1 <= bondist:
            # 判断 point 是否在所有交集里
            for i in range(nlandmark):
                D2 = calcdist(lonlats[i, 0], lonlats[i, 1], point[0], point[1])
                if D2 <= delays[i]:
                    count += 1
        # 如果点 j 到所有的 landmark 之间的距离都小于 landmark 到 target 的距离，
        # 那么该点就属于 goodpoint，存入 ficrslt，用于最后的绘图
        if count == nlandmark:
            gp.write("".join([str(point)[1:-1], "\n"]))
            test += 1
    # 如果一个都没有……
    if test == 0:
        gp.write('')

    gp.close()
    return

if __name__ == '__main__':
    TEST_TARGET = '60.247.13.2'
    TEST_DELAYS = np.array([13.0, 10.0, 7.0, 6.0])
    TEST_LONLATNODES = np.array([
        [116.2068, 39.9074],
        [116.3084, 39.9305],
        [116.3826, 39.8583],
        [116.4104, 39.9050]
    ])
    goodpoints(TEST_TARGET, TEST_DELAYS, TEST_LONLATNODES)
