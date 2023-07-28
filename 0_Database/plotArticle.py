#!/usr/bin/python3

import numpy as np
import glob
import os
import os.path
import sys
from sklearn.mixture import GaussianMixture as GM
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib as mpl

from shapely.geometry import Point
from shapely.ops import unary_union
import matplotlib.patches as ptc
from matplotlib.colors import ListedColormap

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
    #for nsig in range(1, 4):
    nsig = 3
    ax.add_patch(Ellipse(position, nsig * width, nsig * height,
                             angle, **kwargs))
    ax.scatter(position[0], position[1], marker='o', color='black', s=15)


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
# for s in range(len(testStructure)):
#     ax.scatter(testStructure[s][:,2], testStructure[s][:,4], c='black', marker='x', s=16.0, linewidths=1)


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

ax = plt.gca()






####################################################
####### PLOT KMEANS INITIALIZATION ###################
def readQ(input_file):
    pos_tmp=np.genfromtxt(input_file,skip_header=9)
    q_tmp=pos_tmp[:,6:]
    return q_tmp

categories = []
testStructure = []
numDirs = 0
X = np.empty((0,14))
for base, dirs, files in os.walk('QValues/'):
    numDirs += len(dirs)
    for directories in dirs:
        list_file = sorted(glob.glob(base+directories+"/*Q.trj"))
        if len(list_file) == 0:
            numDirs -= 1
            continue
        for input_file in list_file:
            X = np.append(X, readQ(input_file), axis=0)
        testStructure.append(readQ(list_file[0]))
        categories.append(directories)
    
print('Number of clusters used in the model: ', numDirs)

model = GM(numDirs, init_params='k-means++', max_iter=0, tol=1e-9)
model.fit(X)
means = model.means_
covars = model.covariances_
weights = model.weights_

for s in range(len(means)):
    i = 2
    b = 4
    mcov = np.array([[covars[s][i, i], covars[s][i, b]], [covars[s][b, i], covars[s][b, b]]])
    draw_ellipse(np.array([means[s][i], means[s][b]]), mcov, edgecolor='black', facecolor='None', linewidth=1.0)



####################################################
####### PLOT CLASSIFICATION AREAS ###################
# plt.figure()
# ax = plt.gca()

# #Creating and training the GM model
# dirModel = '../0_GMModel/Voronoi/'
# model = GM(8, n_init=100)
# means = np.load(dirModel+'means.npy')
# covar = np.load(dirModel+'covariances.npy')
# model.precisions_cholesky_ = np.linalg.cholesky(np.linalg.inv(covar))
# model.weights_ = np.load(dirModel+'weights.npy')
# model.means_ = means
# model.covariances_ = covar

# means = model.means_
# covars = model.covariances_
# weights = model.weights_


# # for s in range(len(means)):
# #     i = 2
# #     b = 4
# #     mcov = np.array([[covars[s][i, i], covars[s][i, b]], [covars[s][b, i], covars[s][b, b]]])
# #     draw_ellipse(np.array([means[s][i], means[s][b]]), mcov, edgecolor='black', facecolor='None', linewidth=1.0)

# covs2d = []
# means2d = []
# weights2d = weights

# for s in range(len(means)):
#     i = 2
#     b = 4
#     covs2d.append(np.array([[covars[s][i, i], covars[s][i, b]], [covars[s][b, i], covars[s][b, b]]]))
#     means2d.append(np.array([means[s][i], means[s][b]]))
# covs2d = np.array(covs2d)
# means2d = np.array(means2d)

# model2d = GM(8, n_init=100)
# model2d.precisions_cholesky_ = np.linalg.cholesky(np.linalg.inv(covs2d))
# model2d.weights_ = weights2d
# model2d.means_ = means2d
# model2d.covariances_ = covs2d

# q4 = np.linspace(0.0, 0.8, 1000)
# q6 = np.linspace(0.0, 0.6, 1000)
# X, Y = np.meshgrid(q4, q6)
# Z = np.zeros(X.shape)

# for i in range(len(q4)):
#     for b in range(len(q6)):
#         Z[b][i] = model2d.predict(np.array([q4[i], q6[b]]).reshape(1, -1)) - 0.5

# setZ = set()
# for x in Z:
#     for y in x:
#         setZ.add(y)
# print(setZ)

# p = ax.contourf(X, Y, Z, cmap = ListedColormap(('C6', 'C4', 'C1', 'C7', 'C0', 'C2', 'C5', 'C3')))

# #plt.colorbar(p, ax=ax)

# for s in range(len(means)):
#     i = 2
#     b = 4
#     mcov = np.array([[covars[s][i, i], covars[s][i, b]], [covars[s][b, i], covars[s][b, b]]])
#     draw_ellipse(np.array([means[s][i], means[s][b]]), mcov, edgecolor='black', facecolor='None', linewidth=1.0)

#plt.legend(loc='center left', bbox_to_anchor=(-0.2, 1.4), ncol=3)
# plt.ylabel('$q_6$')
# plt.xlabel('$q_4$')

plt.ylim([0.0, 0.60])
plt.xlim([0.0, 0.8])
# Turn off tick labels
ax.set_yticklabels([])
ax.set_xticklabels([])
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
plt.tight_layout()
plt.savefig('databaseKMeans++.png', bbox_inches='tight')