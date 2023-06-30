import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import glob
import os
import os.path
import sys
import re
import matplotlib.transforms as transforms

import ovito
from ovito.io import *
from ovito.modifiers import *
from ovito.data import *
from ovito.pipeline import *

#Natural sorting
def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]


#Variables growth part
ratio = []
sRatio = []
t2 = []
    
dirSeeds = glob.glob('../3_GrowthRandom/WRZ_390N/Output/T_1375K/seed*/')
dirSeeds.sort(key=natural_keys)
# Get temperature values
files = glob.glob(dirSeeds[0]+'dump*.PROB.trj')
files.sort(key=natural_keys)
for f in files:
    t2.append(int(os.path.basename(f).replace('dump', '').replace('.PROB.trj', ''))/1000+16)
# Read cluster sizes for different seeds and compute average
nSurf = []
nSeed = []
gyrationRadius = []
for dS in dirSeeds:
    files = glob.glob(dS+'dump*.PROB.trj')
    files.sort(key=natural_keys)
    nWRZ = []
    nBCTHBN = []
    nGyr = []
    for f in files:
        #Open file and count number of atoms
        pipeline = import_file(f, multiple_frames=True)
        pipeline.modifiers.append(ExpressionSelectionModifier(expression='pliq < 0.5'))
        pipeline.modifiers.append(ClusterAnalysisModifier(cutoff=4, sort_by_size=True, only_selected=True))
        pipeline.modifiers.append(ExpressionSelectionModifier(expression='pwrz > 0.5  && Cluster == 1'))
        data = pipeline.compute()
        nWRZ.append(data.attributes['ExpressionSelection.count.2'])
        #Open file and count number of atoms
        pipeline = import_file(f, multiple_frames=True)
        pipeline.modifiers.append(ExpressionSelectionModifier(expression='pliq < 0.5'))
        pipeline.modifiers.append(ClusterAnalysisModifier(cutoff=4, sort_by_size=True, only_selected=True))
        pipeline.modifiers.append(ExpressionSelectionModifier(expression='pbct > 0.5  ||  phbn > 0.5  && Cluster == 1'))
        data = pipeline.compute()
        nBCTHBN.append(data.attributes['ExpressionSelection.count.2'])
        #Open file and count number of atoms
        pipeline = import_file(f, multiple_frames=True)
        pipeline.modifiers.append(ExpressionSelectionModifier(expression='pliq < 0.5'))
        pipeline.modifiers.append(ClusterAnalysisModifier(cutoff=4, sort_by_size=True, only_selected=True, compute_gyration=True))
        data = pipeline.compute()
        #cluster_table = data.tables['clusters']
        #print(cluster_table)
        #nGyr.append(cluster_table['Radius of Gyration'][0])
    nSeed.append(nWRZ)
    nSurf.append(nBCTHBN)
    #gyrationRadius.append(nGyr)

nSeed = np.array(nSeed)
nSurf = np.array(nSurf)
division = nSeed / nSurf
root = (nSeed + nSurf) ** (1./3.)
ratio = np.mean(division, axis=0)
sRatio = np.std(division, axis=0)
radius = np.mean(root, axis=0)
sRadius = np.std(root, axis=0)
#gyration = np.mean(gyrationRadius, axis=0)
#sGyration = np.std(gyrationRadius, axis=0)

plt.figure()
plt.title('Ratio of $N_{WRZ}$ by $N_{BCT-HBN}$')
plt.ylabel('$N_{seed} / N_{surf}$')
plt.xlabel(r'$N^{\frac{1}{3}}$')
plt.errorbar(radius, ratio, xerr=sRadius, yerr=sRatio, fmt='o', ecolor='gray', lw=0.75, capsize=3, capthick=0.75)

# fig, ax1 = plt.subplots()
# plt.title('Ratio of $N_{WRZ}$ by $N_{BCT-HBN}$')
# ax1.set_xlabel('$t$ / ps')
# ax1.set_ylabel('$N_{seed} / N_{surf}$', color='C0')
# line = ax1.plot(t2, ratio)
# ax1.fill_between(t2, ratio+sRatio, ratio-sRatio, alpha=0.25, color=line[0]._color)
# ax1.tick_params(axis='y', labelcolor='C0')

# ax2 = ax1.twinx()

# ax2.set_ylabel(r'$N^{\frac{1}{3}}$', color='C1')
# line = ax2.plot(t2, radius, c='C1')
# ax2.fill_between(t2, radius+sRadius, radius-sRadius, alpha=0.25, color=line[0]._color)
# ax2.tick_params(axis='y', labelcolor='C1')

# line = ax2.plot(t2, gyration, c='C2')
# ax2.fill_between(t2, gyration+sGyration, gyration-sGyration, alpha=0.25, color=line[0]._color)

#plt.legend()
plt.tight_layout()
plt.savefig('bigWRZgrowing.png')


# #Plotting NVT simulations for different random seeds
# nTs = dict()
# t2s = dict()
# tempDirs = glob.glob('outputVoronoi/step900Seeds/*')
# tempDirs.sort(key=natural_keys)
# for tt in tempDirs:
#     nTs[tt] = []
#     t2s[tt] = []
# for dir in tempDirs:
#     files = glob.glob(dir+'/dump*.PROB.trj')
#     files.sort(key=natural_keys)
#     for f in files:
#         #Getting the number of column for the probability
#         header = ''
#         with open(d+'/step'+str(i)+'/dump'+str(b*1000)+'.PROB.trj') as file:
#             for item in file:
#                 if 'ITEM: ATOMS' in item:
#                     header = item.split(' ')
#                     break
#         index = 0
#         for s in range(len(header)):
#             if 'pLIQ' in header[s]:
#                 index = s - 2
#                 break

#         t2s[dir].append(int(os.path.basename(f).replace('dump', '').replace('.PROB.trj', ''))+16000)
#         #Open file and count number of atoms
#         a = np.genfromtxt(f, skip_header=9)
#         nTs[dir].append(sum(x < 0.5 for x in a[:, index]))
# seeds = []
# for l in nTs:
#     seeds.append(nTs[l])
# averageSeeds = np.average(seeds, axis=0)
# deviationSeeds = np.std(seeds, axis=0)