#!/usr/bin/python3

import numpy as np
import glob
import os
import os.path
import sys
from sklearn.mixture import GaussianMixture as GM
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

from shapely.geometry import Point
from shapely.ops import unary_union
import matplotlib.patches as ptc

def draw_ellipse(position, covariance, ax=None, **kwargs):
    """Draw an ellipse with a given position and covariance"""
    ax = ax or plt.gca()
    
    # Convert covariance to principal axes
    if covariance.shape == (2, 2):
        U, s, Vt = np.linalg.svd(covariance)
        angle = np.degrees(np.arctan2(U[1, 0], U[0, 0]))
        width, height = 2 * np.sqrt(s)
    else:
        angle = 0
        width, height = 2 * np.sqrt(covariance)
    
    # Draw the Ellipse
    for nsig in range(1, 4):
        ax.add_patch(Ellipse(position, nsig * width, nsig * height,
                             angle, **kwargs))


def readQ(input_file):
    pos_tmp=np.genfromtxt(input_file,skip_header=9)
    q_tmp=pos_tmp[:,6:]
    return q_tmp

testStructure = []
dataPoints = dict()
for base, dirs, files in os.walk('QValues/'):
    for directories in dirs:
        if 'surf' in directories:
            continue
        list_file = sorted(glob.glob(base+directories+"/*Q.trj"))
        if len(list_file) == 0:
            continue
        dataPoints[directories] = np.empty((0, 14))
        for input_file in list_file:
            dataPoints[directories] = np.append(dataPoints[directories], readQ(input_file), axis=0)
        if 'LIQ' in list_file[0]:
            continue
        print(list_file[0])
        testStructure.append(readQ(list_file[0]))

for key in dataPoints:
    print(key)


plt.figure()
#plt.title('Database plot (Averaged)')
ax = plt.gca()
# for key in dataPoints:
#     ax.scatter(dataPoints[key][:,2], dataPoints[key][:,6], label=key, alpha=0.25)
for s in range(len(testStructure)):
    ax.scatter(testStructure[s][:,2], testStructure[s][:,4], c='black', marker='x', s=16.0, linewidths=1)


polygons = dict()
for key in dataPoints:
    p = [Point(dataPoints[key][i,2], dataPoints[key][i,4]).buffer(0.002) for i in range(len(dataPoints[key]))]
    # print(len(p))
    # print(type(p[0]))
    # print(p[0].is_valid)
    # for poly in p:
    #     if not p[0].is_valid:
    #         print('non valid')
    #         exit()
    polygons[key] = unary_union(p)
# print(type(polygons['BCT']))

colorCount = 0
first = True
for key in polygons:
    for polygon in polygons[key]:
        if first:
            polygon = ptc.Polygon(np.array(polygon.exterior), facecolor="C"+str(colorCount), lw=0, alpha=0.6, label=key)
            first = False
        else:
            polygon = ptc.Polygon(np.array(polygon.exterior), facecolor="C"+str(colorCount), lw=0, alpha=0.6)
        ax.add_patch(polygon)
    first = True
    colorCount += 1

plt.legend(loc='center left', bbox_to_anchor=(-0.2, 1.4), ncol=3)
plt.ylabel('$q_6$')
plt.xlabel('$q_4$')
plt.ylim([-0.01, 0.60])
plt.xlim([-0.01, 0.78])
plt.tight_layout()
plt.savefig('q46VoronoiLegend.png')

# fig, axs = plt.subplots(14, 14, figsize=(100, 100))
# fig.suptitle('Database plot [q2 ... q8] (Voronoi Cutoff)', fontsize=140)
# for i in range(14):
#     for b in range(14):
#         if b < i:
#             continue
#         for key in dataPoints:
#             axs[i, b].scatter(dataPoints[key][:,i], dataPoints[key][:,b], label=key)
#         for s in range(len(testStructure)):
#             axs[i, b].scatter(testStructure[s][:,i], testStructure[s][:,b], c='black', marker='x')
#         axs[i, b].legend()
#         if i+2 < 9:
#             axs[i, b].set_xlabel('q'+str(i+2))
#         else:
#             axs[i, b].set_xlabel('q'+str(i+2-7)+' mono')
#         if b+2 < 9:
#             axs[i, b].set_ylabel('q'+str(b+2))
#         else:
#             axs[i, b].set_ylabel('q'+str(b+2-7)+' mono')
#         print(i, b)
# fig.savefig('databaseQAverage.png', bbox_inches='tight')
