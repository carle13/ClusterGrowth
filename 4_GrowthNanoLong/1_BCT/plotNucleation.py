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

# ranges = [[100,700], [200,845], [100,700], [150,800], [150,850]]
# posV = [[137,237,600], [280,500,845], [131,250,587], [201,400,800], [161,302,850]]

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


directories = sorted(glob.glob('BCT_*N/Output/', recursive=True))
directories.sort(key=natural_keys)

plot = 0
for dir in directories:
    ax = axs[plot]
    for d in probs:
        nC[d].append([])
    for s in range(1, 6):
        seed = 11111 * s
        for d in probs:
            vals = []
            for b in range(85):
                if d == 'pbct' and s == 1:
                    t[plot].append(b*1000)
                #Count particles belonging to the cluster
                pipeline = import_file(dir+'T_1000K/seed'+str(seed)+'/dump'+str(b*1000)+'.PROB.trj', multiple_frames=True)
                pipeline.modifiers.append(ExpressionSelectionModifier(expression='pbct > 0.5 || pwrz > 0.5 || phbn > 0.5'))
                pipeline.modifiers.append(ClusterAnalysisModifier(cutoff=4, sort_by_size=True, only_selected=True))
                pipeline.modifiers.append(ExpressionSelectionModifier(expression=d+' > 0.5'))
                data = pipeline.compute()
                vals.append(data.attributes['ExpressionSelection.count.2'])
            nC[d][plot].append(vals)
    t[plot] = (np.array(t[plot]) / 1000) + 16


    ax.set_ylim([-30, 650])
    ax.set_xlabel('t [ps]', labelpad=1)
    if plot == 0:
        ax.set_ylabel('$N$', labelpad=1)
    else:
        ax.set_yticklabels([])

    nTotal = []
    for j in range(len(nC['pbct'][plot])):
        nTotal.append(np.array(nC['pbct'][plot][j]) + np.array(nC['pwrz'][plot][j]) + np.array(nC['phbn'][plot][j]))

    colorsRGB = [(228/255,26/255,28/255), (77/255,175/255,74/255), (55/255,126/255,184/255)]
    indexColor = 0
    for d in probs:
        average = np.mean(nC[d][plot], axis=0)
        standDev = np.std(nC[d][plot], axis=0)
        line = ax.plot(t[plot], average, label=d.replace('p', '').upper(), c=colorsRGB[indexColor], lw=1)
        ax.fill_between(t[plot], average+standDev, average-standDev, alpha=0.25, color=line[0]._color)
        indexColor += 1
    
    aveTot = np.mean(nTotal, axis=0)
    stdTot = np.std(nTotal, axis=0)
    line = ax.plot(t[plot], aveTot, label='Total Cluster', c='black', lw=1)
    ax.fill_between(t[plot], aveTot+stdTot, aveTot-stdTot, alpha=0.25, color=line[0]._color)
    # ax.axvline(posV[s-1][0], c='C5', ls='--')
    # ax.axvline(posV[s-1][1], c='C6', ls='--')
    # ax.axvline(posV[s-1][2], c='C7', ls='--')
    plot += 1
#plt.axhline(nC[d][-1], 0, 16, ls='--', color='black')
#plt.text(16-0.15, nC[d][-1]+1, '$N_c = '+str(nC[d][-1])+'$', ha='right')
plt.legend(ncol=4, loc='upper center', bbox_to_anchor=(-2.4, 1.3), handlelength=1)
plt.subplots_adjust(left=0.065, right=0.985, top=0.8, bottom=0.205, wspace = 0.3)
#plt.tight_layout()
#plt.ylim([0, 300])
fig.savefig('compositionLong.pdf', pad_inches = 0.05)