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

plotDir = '4_GrowthResetVel/'
print(len(sys.argv))
if len(sys.argv) > 1:
    plotDir = os.path.join(sys.argv[1], '')

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


directories = sorted(glob.glob(plotDir+'*/Output/', recursive=True))
directories.sort(key=natural_keys)

#Variables relaxation part
nC = dict()
t = []

#Variables growth part
nT = dict()
t2 = dict()

for d in directories:
    #Get directory of corresponding relaxation
    dirRelax = d.replace(plotDir, '2_RelaxDnvt/')
    crystal, inserted = re.findall(r"\w+N", dirRelax)[0].split('_')
    dirRelax = dirRelax.replace(crystal+'_'+inserted, crystal+'_1350K')
    dirRelax = dirRelax + 'N_'+inserted[:-1]+'/'
    nC[dirRelax] = []
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
    t = np.array(t) / 1000

    #Get values for title and figure name
    crystal, temperature = re.findall(r"\w+K", dirRelax)[0].split('_')
    _, inserted = re.findall(r"N_\w+", dirRelax)[0].split('_')
    #Draw plot
    plt.figure()
    plt.title('Cluster size '+crystal+' (Relaxed at '+temperature+')\nInserted atoms: '+inserted)
    plt.xlabel('t / ps')
    plt.ylabel('N')
    plt.axvline(8, ls='-.', color='black', alpha=0.5)
    plt.axvline(16, ls='-.', color='black', alpha=0.5)
    ax = plt.gca()
    trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)
    plt.text(8-6, 0.99, 'Step 1', va='top', transform=trans)
    plt.text(16-6, 0.99, 'Step 2', va='top', transform=trans)
    plt.plot(t, nC[dirRelax])
    plt.axhline(nC[dirRelax][-1], 0, 50, ls='--', color='black')
    plt.text(16-0.2, nC[dirRelax][-1]+3, '$N_c = '+str(nC[dirRelax][-1])+'$', ha='right')

    
    #Plotting NVT simulations at different temperatures
    dirTemp = glob.glob(d+'T_*K/')
    dirTemp.sort(key=natural_keys)
    for dT in dirTemp:
        files = glob.glob(dT+'dump*.PROB.trj')
        files.sort(key=natural_keys)
        nT[dT] = []
        t2[dT] = []
        for f in files:
            t2[dT].append(int(os.path.basename(f).replace('dump', '').replace('.PROB.trj', ''))/1000+16)
            #Open file and count number of atoms
            pipeline = import_file(f, multiple_frames=True)
            pipeline.modifiers.append(ExpressionSelectionModifier(expression='pliq < 0.5'))
            pipeline.modifiers.append(ClusterAnalysisModifier(cutoff=4, sort_by_size=True, only_selected=True))
            data = pipeline.compute()
            nT[dT].append(data.attributes['ClusterAnalysis.largest_size'])

        plt.plot(t2[dT], nT[dT], label=os.path.basename(dT[:-1]).replace('T_', ''))
        # print('plotting: '+str(dT)+'         last value: '+str(nT[dT]))
        # plt.plot(t2['outputVoronoi/step900K'], nT['outputVoronoi/step900K'], label='900K')
        # plt.plot(t2s['outputVoronoi/step900Seeds/900K_1'], averageSeeds, label='Seeds 900')
        # plt.fill_between(t2s['outputVoronoi/step900Seeds/900K_1'], averageSeeds-deviationSeeds, averageSeeds+deviationSeeds, alpha=0.5, color='green')
    plt.legend()
    plt.savefig(d+'growthR800K.png')


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