#!usr/bin/env python3
# -*- coding:utf-8 -*-

"""最终对目标实现定位

"""


import os
import numpy as np
from Polygon import *
from Utils import *


def locate(target, delays, lonlats):
    # nblandmark = len(delays)

    # rtt contains the min RTT between target host and landmarks
    # rttNnodes <- rtt[goodlandRTT]
    # posminrtt <- which.min(rttNnodes)
    # 这两个好像并没有被用到，所以暂时被禁用了

    # we output the result in png format, and create the file .png
    # figure = "Figure/location_estimate_of_"
    # figure_ip = "".join([figure, target, ".png"])

    # boolpolygon is a variable to test whether it has a polygon or not
    # boolpolygon = 0

    pos = np.argmin(delays)
    goodpoint = "Goodpoint/zonecross-"
    goodpoint_ip = "".join([goodpoint, target])

    ploc = 'Results/estimated-loc.dat'
    pcen = 'Results/errors.dat'
    ppol = 'Results/area-polygons.dat'
    floc = open(ploc, 'a')
    fcen = open(pcen, 'a')
    fpol = open(ppol, 'a')

    # 若 goodpoint 文件不存在，则无法定位，直接 return
    if not os.path.exists(goodpoint_ip):
        dumpfile(floc, fcen, fpol, target, 0, 0, 0, 0)
        return
    # 若文件存在，但是实际上并没有交点，不存在多边形
    if os.stat(goodpoint_ip).st_size == 0:
        CX = lonlats[pos, 0]
        CY = lonlats[pos, 1]
        RADIUS = delays[pos]
        AREA = np.pi * RADIUS ** 2
        dumpfile(floc, fcen, fpol, target, CX, CY, RADIUS, AREA)
        return

    # 若文件存在，且存在交点
    goodpts = np.genfromtxt(goodpoint_ip)
    if len(goodpts) >= 3:
        # 计算起点，交点中经度最小的点 (x0,y0)，到其他所有交点的角度
        # angles = np.zeros((len(goodpts), 1))
        # for j in range(ngoodpts):
        #     angle[j] = calcangle(xzero, yzero, goodpts[j, 0], goodpts[j, 1])
        # 对角度进行排序
        # angles = sorted(angles)
        # count = 0
        # while count <= ngoodpts:
        #     for k in range(1, ngoodpts):
        #         ang = rad2deg(calcangle(xzero, yzero,
        #                                 goodpts[k, 0],
        #                                 goodpts[k, 1]))
        #         if ang == angles:
        #             ptorder[count, 0] = goodpts[k, 0]
        #             ptorder[count, 1] = goodpts[k, 1]
        #     count = count + 1
        # 原作者的这段代码是为了把点进行排序，python 里面可以很容易实现
        posmin = np.argmin(goodpts[:, 0])
        xzero, yzero = goodpts[posmin]
        ptorder = np.array(
            sorted(goodpts,
                   key=lambda p: calcangle(xzero, yzero, p[0], p[1])))

        # calculate (Cx, Cy) the center of the polygon, R = sqrt(area /
        # pi). (lon, lat)->(dlon * lon, dlat * lat) dlon and dlat represent the
        # conversion coefficients for each polygon
        ptkm = degre2km(ptorder)
        # calculate the centroid of the polygon and the radius of the
        # circle of this polygon, airepolygone is the function that
        # returns the centroid, radius and area of the polygon, and
        # traces the vertices of the polygon
        CX, CY, RADIUS, AREA = polygon(ptorder, ptkm)

        # 如果面积不为 0
        if AREA != 0 and AREA != 1:
            # we trace the centroid of the polygon (estimate of the location of
            # the point)
            # points(CxCyRAIRE[1], CxCyRAIRE[2], pch=4, col="red", cex=0.5)
            dumpfile(floc, fcen, fpol, target, CX, CY, RADIUS, AREA)
            return
            # we have a polygon we put the variable a 1
            # boolpolygon = 1
        # 如果外边界不是一个凸包
        elif AREA == 1:
            AREA = polyarea(goodpts)
            RADIUS = np.sqrt(abs(AREA) / np.pi)
            CX, CY = polycent(goodpts)
            dumpfile(floc, fcen, fpol, target, CX, CY, RADIUS, AREA)
            return
        # 如果多点共线，不存在多边形
        else:
            xzero, yzero = goodpts[0]
            DIAMETER = 0
            for point in goodpts:
                tmp = calcdist(xzero, yzero, point[0], point[1])
                if DIAMETER == 0 or DIAMETER > tmp:
                    DIAMETER = tmp
            RADIUS = DIAMETER / 2
            AREA = np.pi * RADIUS ** 2
            CX, CY = np.mean(goodpts, 0)
            dumpfile(floc, fcen, fpol, target, CX, CY, RADIUS, AREA)
            return
            # boolpolygon = 2
    # 只有两个交点，无法构成一个多边形
    elif len(goodpts) == 2:
        RADIUS = calcdist(goodpts[0, 0], goodpts[0, 1],
                          goodpts[1, 0], goodpts[1, 1]) / 2
        AREA = np.pi * RADIUS ** 2
        CX, CY = np.mean(goodpts, 0)
        dumpfile(floc, fcen, fpol, target, CX, CY, RADIUS, AREA)
        return
        # boolpolygon = 2
    # 只有一个交点
    else:
        # we have a point belonging to all the circles
        CX = lonlats[pos, 0]
        CY = lonlats[pos, 1]
        RADIUS = delays[pos]
        AREA = np.pi * RADIUS ** 2
        dumpfile(floc, fcen, fpol, target, CX, CY, RADIUS, AREA)
    return


def dumpfile(floc, fcen, fpol, *args):
    target = args[0]
    CX = args[1]
    CY = args[2]
    RADIUS = args[3]
    AREA = args[4]
    floc.write(", ".join([target, str(CX), str(CY)]) + "\n")
    fcen.write(", ".join([target, str(RADIUS)] + "\n"))
    fpol.write(", ".join([target, str(AREA)] + "\n"))
    floc.close()
    fcen.close()
    fpol.close()


if __name__ == '__main__':
    TEST_TARGET = '60.247.13.2'
    TEST_DELAYS = np.array([13.0, 10.0, 7.0, 6.0])
    TEST_LONLATNODES = np.array([
        [116.2068, 39.9074],
        [116.3084, 39.9305],
        [116.3826, 39.8583],
        [116.4104, 39.9050]
    ])
    # rtt = 12
    locate(TEST_TARGET, TEST_DELAYS, TEST_LONLATNODES)
