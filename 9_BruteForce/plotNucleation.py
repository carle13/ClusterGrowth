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
matplotlib.rcParams.update({'font.size': 9.98655939})
fig = plt.figure()
fig.set_figheight(4.2*cm)
fig.set_figwidth(17*cm)

axs = fig.subplots(1, 5)

ranges = [[0,700], [0,950], [0,700], [0,900], [0,950]]
posV = [[137,237,324,600], [280,348,500,845], [131,177,250,587], [201,316,400,800], [161,244,302,850]]

for s in range(1, 6):
    ax = axs[s - 1]
    seed = 11111 * s
    for d in probs:
        vals = []
        for b in range(ranges[s-1][0], ranges[s-1][1], 20):
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


    ax.set_xlabel('t [ps]', labelpad=1)
    ax.set_ylim([-30, 650])
    if s-1 == 0:
        ax.set_ylabel('$N$', labelpad=1)
    else:
        ax.set_yticklabels([])

    nTotal = np.array(nC['pbct'][s-1]) + np.array(nC['pwrz'][s-1]) + np.array(nC['phbn'][s-1])

    colorsRGB = [(228/255,26/255,28/255), (77/255,175/255,74/255), (55/255,126/255,184/255)]
    indexColor = 0
    for d in probs:
        ax.plot(t[s-1], nC[d][s-1], label=d.replace('p', '').upper(), c=colorsRGB[indexColor], lw=1)
        indexColor += 1

    ax.plot(t[s-1], nTotal, label='Total Cluster', c='black', lw=1)
    ax.axvline(posV[s-1][0], c='C7', ls='--', dashes=(2, 2))
    ax.axvline(posV[s-1][1], c='C7', ls='--', dashes=(2, 2))
    ax.axvline(posV[s-1][2], c='C7', ls='--', dashes=(2, 2))
    ax.axvline(posV[s-1][3], c='C7', ls='--', dashes=(2, 2))
#plt.axhline(nC[d][-1], 0, 16, ls='--', color='black')
#plt.text(16-0.15, nC[d][-1]+1, '$N_c = '+str(nC[d][-1])+'$', ha='right')
plt.legend(ncol=4, loc='upper center', bbox_to_anchor=(-2.5, 1.4), handlelength=1)
plt.subplots_adjust(left=0.065, right=0.985, top=0.8, bottom=0.205, wspace = 0.3)
#plt.tight_layout()
#plt.ylim([0, 300])
fig.savefig('compositionBrute.pdf', pad_inches = 0.05)