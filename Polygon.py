#!usr/bin/env python3
# -*- coding:utf-8 -*-


import numpy as np
from Utils import *
from scipy.spatial import ConvexHull


def polygon(ptorder, ptkm):
    """
    A = sigma(x_i * y_i+1 - x_i+1 * y_i) / 2
    Cx = sigma((x_i + x_i+1) * (x_i * y_i+1 - x_i+1 * y_i)) / 6A
    Cy = sigma((y_i + y_i+1) * (x_i * y_i+1 - x_i+1 * y_i)) / 6A
    """
    if len(ConvexHull(ptorder).vertices) != len(ptorder):
        return [0, 0, 0, 1]

    AREA = polyarea(ptkm[:-1])
    Cx, Cy = polycent(ptkm[:-1]) / ptkm[-1]
    R = np.sqrt(abs(AREA) / np.pi)

    return [Cx, Cy, R, abs(AREA)]


def polyarea(pol):
    x = pol[:, 0]
    y = pol[:, 1]
    return 0.5 * (np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def polycent(pol):
    x = pol[:, 0]
    y = pol[:, 1]
    a = polyarea(pol)
    Cx = np.sum((x + np.roll(x, -1)) *
                (x * np.roll(y, -1) - y * np.roll(x, -1))) / (6 * a)
    Cy = np.sum((y + np.roll(y, -1)) *
                (x * np.roll(y, -1) - y * np.roll(x, -1))) / (6 * a)
    return [Cx, Cy]
