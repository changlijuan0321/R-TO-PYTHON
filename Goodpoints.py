#!usr/bin/env python3
# -*- coding:utf-8 -*-

"""Newcrosscircle.R 找到了所有圆两两之间的交集，Goodpoints.R 负责找到出现在所有
圆中的点 This function seeks the "good points" among the crossing points found
before. These good points represent the intersection region of all
circles. Good points are any points that are located in all circles. These good
points allow us to find the area where we expect that the target host is
located. These good points will be the vertices of the polygon. Note that the
polygon is the approximation of this intersection region. The longitude and the
latitude of good points are saved in a file whose name is the IP address of the
given target host and the file is located in the folder "Bonpoint". It should
be noted that if we have two good points, the intersection region will be the
circle whose diameter is the geographic distance of these two points.

"""

import sys
import os
import numpy as np
from io import StringIO
from Utils import *


def goodpoints(target, delays, lonlats):
    nlandmark = len(delays)
    crossingpoint = "Crossingpoint/crossingpoint - "
    nomfic = "".join([crossingpoint, target])
    bonpoint = "Bonpoint/zonecross - "
    ficrslt = "".join([bonpoint, target])

    if os.path.exists(ficrslt):
        os.remove(ficrslt)

    with open(nomfic, 'r') as f:
        npoints = np.arange(0, nlandmark)
        # seeks the smallest radius
        mindelaysNodes = min(delays)
        # the least distance between any landmark toward a target.
        pos = delays.index(mindelaysNodes)
        bondist = delays[pos]
        c = StringIO(f.read())
        ptcrossing = np.loadtxt(c, delimiter=',')
        long = len(ptcrossing)
        test = 0
        for j in range(0, long):
            compteur = 0
            D1 = calculdist(lonlats[pos, 0], lonlats[pos, 1],
                            ptcrossing[j, 0], ptcrossing[j, 1])
            if D1 <= bondist:
                # seek if the give b point is inside all circles
                for z in npoints:
                    D2 = calculdist(lonlats[z, 0], lonlats[z, 1],
                                    ptcrossing[j, 0], ptcrossing[j, 1])
                    if D2 <= delays[z]:
                        compteur += 1
            # 如果点 j 到所有的 landmark 之间的距离都小于 landmark 到 target 的
            # 距离，那么该点就属于 goodpoint，存入 ficrslt，用于最后的绘图
            if compteur == len(npoints):
                with open(ficrslt, 'a') as f:
                    f.write(str(ptcrossing[j, :])[1:-1])
                    test += 1
        if test == 0:
            with open(ficrslt, 'a') as f:
                    f.write('')

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
