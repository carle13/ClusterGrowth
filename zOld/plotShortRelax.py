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

directories = sorted(glob.glob('3_ShortRelaxation/*/Output/N_*/', recursive=True))
nC = dict()
t = []
for d in directories:
    nC[d] = []

for d in directories:
    t = []
    for i in range(1, 3):
        for b in range(5):
            if i == 1:
                t.append(b*1000)
            elif i == 2:
                t.append(b*1000 + 4000)
            #Count particles belonging to the cluster
            pipeline = import_file(d+'step'+str(i)+'/dump'+str(b*1000)+'.PROB.trj', multiple_frames=True)
            pipeline.modifiers.append(ExpressionSelectionModifier(expression='pliq < 0.5'))
            pipeline.modifiers.append(ClusterAnalysisModifier(cutoff=4, sort_by_size=True, only_selected=True))
            data = pipeline.compute()
            nC[d].append(data.attributes['ClusterAnalysis.largest_size'])
    t = np.array(t) / 1000

    #Get values for title and figure name
    crystal, temperature = re.findall(r"\w+K", d)[0].split('_')
    _, inserted = re.findall(r"N_\w+", d)[0].split('_')
    #Draw plot
    plt.figure()
    plt.title('Cluster size '+crystal+' (Relaxed at '+temperature+')\nInserted atoms: '+inserted)
    plt.xlabel('t / ps')
    plt.ylabel('N')
    plt.axvline(4, ls='-.', color='black', alpha=0.5)
    plt.axvline(8, ls='-.', color='black', alpha=0.5)
    ax = plt.gca()
    trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)
    plt.text(4-1.8, 0.99, 'Step 1', va='top', transform=trans)
    plt.text(8-1.8, 0.99, 'Step 2', va='top', transform=trans)
    plt.plot(t, nC[d])
    plt.axhline(nC[d][-1], 0, 8, ls='--', color='black')
    plt.text(8-0.15, nC[d][-1]+1, '$N_c = '+str(nC[d][-1])+'$', ha='right')
    #plt.legend()
    plt.savefig(d+'relaxation'+temperature+'.png')
    plt.close()