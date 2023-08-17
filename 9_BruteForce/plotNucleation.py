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

probs = ['pbct', 'pwrz', 'phbn']
nC = dict()
t = [None]*5
for i in range(len(t)):
    t[i] = []
for d in probs:
    nC[d] = []

cm = 1/2.54  #centimeters in inches
matplotlib.rcParams.update({'font.size': 7.370078})
fig = plt.figure()
fig.set_figheight(4.2*cm)
fig.set_figwidth(17*cm)

axs = fig.subplots(1, 5)

ranges = [[0,300], [0,300], [0,300], [0,300], [0,300]]
posV = [[100,200,250], [100,200,250], [100,200,250], [100,200,250], [100,200,250]]

for s in range(1, 6):
    ax = axs[s - 1]
    seed = 11111 * s
    for d in probs:
        vals = []
        for b in range(ranges[s-1][0], ranges[s-1][1], 5):
            if d == 'pbct':
                t[s-1].append(b*1000)
            #Count particles belonging to the cluster
            pipeline = import_file('2_500000steps/Output/T_1000K/seed'+str(seed)+'/dump'+str(b*1000)+'.PROB.trj', multiple_frames=True)
            # pipeline.modifiers.append(ExpressionSelectionModifier(expression='pbct > 0.5 || pwrz > 0.5 || phbn > 0.5'))
            # pipeline.modifiers.append(ClusterAnalysisModifier(cutoff=4, sort_by_size=True, only_selected=True))
            pipeline.modifiers.append(ExpressionSelectionModifier(expression=d+' > 0.5'))
            data = pipeline.compute()
            vals.append(data.attributes['ExpressionSelection.count'])
        nC[d].append(vals)
    t[s-1] = np.array(t[s-1]) / 1000


    ax.set_xlabel('$t$ / ps')
    if s-1 == 0:
        ax.set_ylabel('$N$')

    nTotal = np.array(nC['pbct'][s-1]) + np.array(nC['pwrz'][s-1]) + np.array(nC['phbn'][s-1])

    for d in probs:
        ax.plot(t[s-1], nC[d][s-1], label=d.replace('p', '').upper())

    ax.plot(t[s-1], nTotal, label='Total Cluster')
    ax.axvline(posV[s-1][0], c='r', ls='--')
    ax.axvline(posV[s-1][1], c='g', ls='--')
    ax.axvline(posV[s-1][2], c='b', ls='--')
#plt.axhline(nC[d][-1], 0, 16, ls='--', color='black')
#plt.text(16-0.15, nC[d][-1]+1, '$N_c = '+str(nC[d][-1])+'$', ha='right')
plt.legend(ncol=4, loc='upper center', bbox_to_anchor=(-2.5, 1.3))
plt.subplots_adjust(left=0.07, right=0.99, top=0.8, bottom=0.2, wspace = 0.35)
#plt.tight_layout()
#plt.ylim([0, 300])
fig.savefig('compositionBrute.pdf',
    pad_inches = 0.05)