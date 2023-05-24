import numpy as np
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


directories = sorted(glob.glob('4_Growth/*/Output/*K/', recursive=True))
directories.sort(key=natural_keys)

#Variables relaxation part
nC = dict()
t = []
for d in directories:
    nC[d] = []

#Variables growth part
nT = dict()
t2 = dict()
for tt in directories:
    nT[tt] = []
    t2[tt] = []

for d in directories:
    #Get directory of corresponding relaxation
    dirRelax = d.replace('4_Growth', '2_Relaxation')
    crystal, temperature = re.findall(r"\w+K", dirRelax)[0].split('_')
    _, inserted = re.findall(r"N_\w+", dirRelax)[0].split('_')
    dirRelax = dirRelax.replace(crystal+'_'+inserted+'N', crystal+'_'+temperature)
    dirRelax = dirRelax.replace('T_'+temperature, 'N_'+inserted)
    #Plotting relaxation part
    t = []
    for i in range(1, 3):
        for b in range(9):
            if i == 1:
                t.append(b*1000)
            elif i == 2:
                t.append(b*1000 + 8000)

            pipeline = import_file(dirRelax+'step'+str(i)+'/dump'+str(b*1000)+'.PROB.trj', multiple_frames=True)
            pipeline.modifiers.append(ExpressionSelectionModifier(expression='pliq < 0.5'))
            pipeline.modifiers.append(ClusterAnalysisModifier(cutoff=4, sort_by_size=True, only_selected=True))
            data = pipeline.compute()
            nC[dirRelax].append(data.attributes['ClusterAnalysis.largest_size'])
    t = np.array(t)

    crystal, temperature = re.findall(r"\w+K", dirRelax)[0].split('_')
    _, inserted = re.findall(r"N_\w+", dirRelax)[0].split('_')
    plt.figure()
    plt.title('Cluster size '+crystal+' (Non-liquid at '+temperature+')\nInserted atoms: '+inserted)
    plt.xlabel('t / ps')
    plt.ylabel('N')
    plt.axvline(8000, ls='-.', color='black', alpha=0.5)
    plt.axvline(16000, ls='-.', color='black', alpha=0.5)
    ax = plt.gca()
    trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)
    plt.text(8000-1800, 0.99, 'Step 1', va='top', transform=trans)
    plt.text(16000-1800, 0.99, 'Step 2', va='top', transform=trans)
    plt.plot(t, nC[d])
    plt.axhline(nC[d][-1], 0, 50000, ls='--', color='black')
    plt.text(16000-150, nC[d][-1]+1, '$N_c = '+str(nC[d][-1])+'$', ha='right')
    #plt.legend()
    # plt.savefig(d+'relaxation'+temperature+'.png')

    #Plotting NVT simulations at different temperatures
    files = glob.glob(d+'dump*.PROB.trj')
    files.sort(key=natural_keys)
    for f in files:
        # #Getting the number of column for the probability
        # header = ''
        # with open(d+'/step'+str(i)+'/dump'+str(b*1000)+'.PROB.trj') as file:
        #     for item in file:
        #         if 'ITEM: ATOMS' in item:
        #             header = item.split(' ')
        #             break
        # index = 0
        # for s in range(len(header)):
        #     if 'pLIQ' in header[s]:
        #         index = s - 2
        #         break

        t2[dir].append(int(os.path.basename(f).replace('dump', '').replace('.PROB.trj', ''))+16000)
        #Open file and count number of atoms
        pipeline = import_file(f, multiple_frames=True)
        pipeline.modifiers.append(ExpressionSelectionModifier(expression='pliq < 0.5'))
        pipeline.modifiers.append(ClusterAnalysisModifier(cutoff=4, sort_by_size=True, only_selected=True))
        data = pipeline.compute()
        nT[d].append(data.attributes['ClusterAnalysis.largest_size'])
        # a = np.genfromtxt(f, skip_header=9)
        # nT[dir].append(sum(x < 0.5 for x in a[:, index]))
    
    #Drawing the plot
    plt.figure()
    plt.title('Cluster size WRZ')
    plt.xlabel('t / ps')
    plt.ylabel('N')
    plt.axvline(8000, ls='--', color='black')
    plt.axvline(16000, ls='--', color='black')
    for l in nC:
        plt.plot(t, nC[l])
    for l in nT:
        plt.plot(t2[l], nT[l], label=os.path.basename(l).replace('step', ''))
    # plt.plot(t2['outputVoronoi/step900K'], nT['outputVoronoi/step900K'], label='900K')
    # plt.plot(t2s['outputVoronoi/step900Seeds/900K_1'], averageSeeds, label='Seeds 900')
    # plt.fill_between(t2s['outputVoronoi/step900Seeds/900K_1'], averageSeeds-deviationSeeds, averageSeeds+deviationSeeds, alpha=0.5, color='green')
    plt.legend()
    plt.savefig('average900K.png')


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