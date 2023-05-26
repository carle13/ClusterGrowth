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

directories = sorted(glob.glob('2_Relaxation/*/Output/N_*/', recursive=True))
nC = dict()
t = []
for d in directories:
    nC[d] = []

for d in directories:
    t = []
    for i in range(1, 3):
        for b in range(9):
            if i == 1:
                t.append(b*1000)
            elif i == 2:
                t.append(b*1000 + 8000)
            # header = ''
            # with open(d+'/Relaxation750K/step'+str(i)+'_750K/dump'+str(b*1000)+'.PROB.trj') as file:
            #     for item in file:
            #         if 'ITEM: ATOMS' in item:
            #             header = item.split(' ')
            #             break
            # index = 0
            # for s in range(len(header)):
            #     if 'pLIQ' in header[s]:
            #         index = s - 2
            #         break

            pipeline = import_file(d+'step'+str(i)+'/dump'+str(b*1000)+'.PROB.trj', multiple_frames=True)
            pipeline.modifiers.append(ExpressionSelectionModifier(expression='pliq < 0.5'))
            pipeline.modifiers.append(ClusterAnalysisModifier(cutoff=4, sort_by_size=True, only_selected=True))
            data = pipeline.compute()
            nC[d].append(data.attributes['ClusterAnalysis.largest_size'])
            # a = np.genfromtxt(d+'/Relaxation750K/step'+str(i)+'_750K/dump'+str(b*1000)+'.PROB.trj', skip_header=9)
            # nC[d].append(sum(x < 0.5 for x in a[:, index]))
    t = np.array(t)

    crystal, temperature = re.findall(r"\w+K", d)[0].split('_')
    _, inserted = re.findall(r"N_\w+", d)[0].split('_')
    plt.figure()
    plt.title('Cluster size '+crystal+' (Non-liquid at '+temperature+')\nInserted atoms: '+inserted)
    plt.xlabel('t / ps')
    plt.ylabel('N')
    plt.axvline(8000, ls='-.', color='black', alpha=0.5)
    plt.axvline(16000, ls='-.', color='black', alpha=0.5)
    ax = plt.gca()
    trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)
    plt.text(8000-1800, 0.01, 'Step 1', va='bottom', transform=trans)
    plt.text(16000-1800, 0.01, 'Step 2', va='bottom', transform=trans)
    plt.plot(t, nC[d])
    plt.axhline(nC[d][-1], 0, 16000, ls='--', color='black')
    plt.text(16000-150, nC[d][-1]+1, '$N_c = '+str(nC[d][-1])+'$', ha='right')
    #plt.legend()
    plt.savefig(d+'relaxation'+temperature+'.png')
    plt.close()