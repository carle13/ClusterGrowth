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


directories = sorted(glob.glob('*/Output/', recursive=True))
directories.sort(key=natural_keys)

#Variables relaxation part
nC = dict()
t = []

cm = 1/2.54  #centimeters in inches
matplotlib.rcParams.update({'font.size': 9.98655939})
fig = plt.figure(figsize=(17*cm, 12*cm))
subfigs = fig.subfigures(3, 2, hspace=0, wspace=0).flatten()

axs = []
for s in subfigs:
    axs.append(s.subplots(1, 2))

#axs1 = subfigs[0].subplots(1, 2)
#subfigs[0].suptitle('BCT')

#axs2 = subfigs[1].subplots(1, 2)
#subfigs[1].suptitle('WRZ')

cB = -1
cW = -1
lims = [50, 75, 55, 55, 55]
lims = list(reversed(lims))
tfig = [None]*5
sfig = [None]*5
tempssub = [None]*5
#Variables growth part
t2 = [None]*5

for i in range(5):
    tfig[i] = [None]*2
    sfig[i] = [None]*2
    t2[i] = [None]*2
    tempssub[i] = [None]*2
    for b in range(2):
        tfig[i][b] = dict()
        sfig[i][b] = dict()
        t2[i][b] = dict()
        tempssub[i][b] = []
indexColor = 0
for d in reversed(directories):
    print(d)
    structure = 10
    if 'BCT' in d:
        cB += 1
        structure = 0
        ax = axs[cB][structure]
        indLim = cB
    elif 'WRZ' in d:
        cW += 1
        structure = 1
        ax = axs[cW][structure]
        indLim = cW
        
    #Get directory of corresponding relaxation
    dirRelax = '../../2_RelaxDnvt/' + d
    crystal, inserted = re.findall(r"\w+N", dirRelax)[0].split('_')
    temp = '_1000K'
    # if int(inserted.replace('N', '')) < 100:
    #     temp = '_1250K'
    dirRelax = dirRelax.replace(crystal+'_'+inserted, crystal+temp)
    dirRelax = dirRelax + 'N_'+inserted[:-1]+'/'
    nC[dirRelax] = []
    if not os.path.exists(dirRelax+'step'+str(1)+'/dump'+str(0)+'.PROB.trj'):
        print('does not exist')
        continue
    #Plotting relaxation part
    t = []
    for i in range(1, 3):
        for b in range(6):
            if i == 1:
                t.append(b*1000)
            elif i == 2:
                t.append(b*1000 + 5000)
            
            pipeline = import_file(dirRelax+'step'+str(i)+'/dump'+str(b*1000)+'.PROB.trj', multiple_frames=True)
            pipeline.modifiers.append(ExpressionSelectionModifier(expression='pbct > 0.5 || pwrz > 0.5 || phbn > 0.5'))
            pipeline.modifiers.append(ClusterAnalysisModifier(cutoff=4, sort_by_size=True, only_selected=True))
            data = pipeline.compute()
            nC[dirRelax].append(data.attributes['ClusterAnalysis.largest_size'])
    t = np.array(t) / 1000

    #Get values for title and figure name
    crystal, temperature = re.findall(r"\w+K", dirRelax)[0].split('_')
    _, inserted = re.findall(r"N_\w+", dirRelax)[0].split('_')
    #Draw plot
    #plt.title('Cluster size '+crystal+' (Relaxed at '+temperature+')\nInserted atoms: '+inserted)
    ax.set_xlabel('t [ps]', labelpad=0)
    ax.axvline(5, ls='-.', color='black', alpha=0.3)
    ax.axvline(10, ls='-.', color='black', alpha=0.3)
    trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)
    #plt.text(8-7, 0.99, 'Step 1', va='top', transform=trans)
    #plt.text(16-7, 0.99, 'Step 2', va='top', transform=trans)
    ax.plot(t, nC[dirRelax], c='black', lw=1)
    ax.axhline(nC[dirRelax][-1], 0, 50, ls='--', color='black')
    #plt.text(0, nC[dirRelax][-1]+25, '$N_c = '+str(nC[dirRelax][-1])+'$', ha='left')
    if 'BCT' in d:
        ax.set_ylabel('$N$', labelpad=1)
        ax.set_title('BCT '+'$N_c = '+str(nC[dirRelax][-1])+'$', fontsize=9.98655939)
    else:
        ax.set_yticklabels([])
        ax.set_title('WRZ '+'$N_c = '+str(nC[dirRelax][-1])+'$', fontsize=9.98655939)

    #Plotting NVT simulations at different temperatures
    dirTemp = glob.glob(d+'T_*K/')
    dirTemp.sort(key=natural_keys)
    for dT2 in dirTemp:
        dT = os.path.basename(dT2[:-1]).replace('T_', '')
        #print(dT)
        dirSeeds = glob.glob(dT2+'seed*/')
        dirSeeds.sort(key=natural_keys)
        tfig[indLim][structure][dT] = []
        sfig[indLim][structure][dT] = []
        t2[indLim][structure][dT] = []
        # Get temperature values
        files = glob.glob(dirSeeds[0]+'dump*.PROB.trj')
        files.sort(key=natural_keys)
        for f in files:
            t2[indLim][structure][dT].append(int(os.path.basename(f).replace('dump', '').replace('.PROB.trj', ''))/1000+10)
        # Read cluster sizes for different seeds and compute average
        nSeed = []
        for dS in dirSeeds:
            files = glob.glob(dS+'dump*.PROB.trj')
            files.sort(key=natural_keys)
            nCluster = []
            for f in files:
                #Open file and count number of atoms
                pipeline = import_file(f, multiple_frames=True)
                pipeline.modifiers.append(ExpressionSelectionModifier(expression='pbct > 0.5 || pwrz > 0.5 || phbn > 0.5'))
                pipeline.modifiers.append(ClusterAnalysisModifier(cutoff=4, sort_by_size=True, only_selected=True))
                data = pipeline.compute()
                nCluster.append(data.attributes['ClusterAnalysis.largest_size'])
            nSeed.append(nCluster)
        tfig[indLim][structure][dT] = np.mean(nSeed, axis=0)
        sfig[indLim][structure][dT] = np.std(nSeed, axis=0)
        tempssub[indLim][structure].append(os.path.basename(dT2[:-1]).replace('T_', ''))
    ax.set_ylim([0, lims[indLim]])


cB = 0
cW = 0
for i in range(5):
    tSubfigure = [element for b in range(len(tempssub[i])) for element in tempssub[i][b]]
    tSubfigure = sorted(list(set(tSubfigure)))
    colorsplasma = plt.cm.plasma(np.linspace(0, 1, len(tSubfigure)))
    for s in range(2):
        ax = axs[i][s]
        for dT in range(len(tempssub[i][s])):
            indexColor = 0
            for t in range(len(tSubfigure)):
                if tSubfigure[t] == tempssub[i][s][dT]:
                    indexColor = t
            key = tSubfigure[indexColor]
            line = ax.plot(t2[i][s][key], tfig[i][s][key], label=tempssub[i][s][dT], c=colorsplasma[indexColor], lw=1)
            ax.fill_between(t2[i][s][key], tfig[i][s][key]+sfig[i][s][key], tfig[i][s][key]-sfig[i][s][key], alpha=0.25, color=line[0]._color)


# for d in reversed(directories):
#     if 'BCT' in d:
#         ax = axs[cB][0]
#         indLim = cB
#         cB += 1
#     elif 'WRZ' in d:
#         ax = axs[cW][1]
#         indLim = cW
#         cW += 1
#     dirTemp = glob.glob(d+'T_*K/')
#     dirTemp.sort(key=natural_keys)
#     for dT in sorted(dirTemp):
#         temp = os.path.basename(dT[:-1]).replace('T_', '')
#         indexColor = 0
#         if temp not in tempsfig[indLim]:
#             tempsfig[indLim].append(temp)
#         for t in range(len(tempsfig[indLim])):
#             if tempsfig[indLim][t] == temp:
#                 indexColor = t
#         colorsplasma = plt.cm.plasma(np.linspace(0, 1, 7))
#         line = ax.plot(t2[dT], nT[dT], label=temp, c='C'+str(indexColor))
#         ax.fill_between(t2[dT], nT[dT]+sT[dT], nT[dT]-sT[dT], alpha=0.25, color=line[0]._color)
#     ax.set_ylim([0, lims[indLim]])


handlesWRZ = []
labelsWRZ = []
for sub in subfigs:
    handlesBCT = []
    labelsBCT = []
    axsinsub = sub.get_axes()
    for ax in axsinsub:
        handles, labels = ax.get_legend_handles_labels()
        indRemove = []
        for l in range(len(labels)):
            if labels[l] in labelsBCT:
                indRemove.append(l)
        handles = [handles[i] for i in range(len(handles)) if i not in indRemove]
        labels = [labels[i] for i in range(len(labels)) if i not in indRemove]
        handlesBCT += handles
        labelsBCT += labels
    if len(labelsBCT) != 0:
        labelsBCT, handlesBCT = zip(*sorted(zip(labelsBCT, handlesBCT), key=lambda t: t[0]))
    sub.legend(handlesBCT, labelsBCT, loc='center left', bbox_to_anchor=(0.7, 0.5), handlelength=1)
# for ax in axsWRZ:
#     handles, labels = ax.get_legend_handles_labels()
#     indRemove = []
#     for l in range(len(labels)):
#         if labels[l] in labelsWRZ:
#             indRemove.append(l)
#     handles = [handles[i] for i in range(len(handles)) if i not in indRemove]
#     labels = [labels[i] for i in range(len(labels)) if i not in indRemove]
#     handlesWRZ += handles
#     labelsWRZ += labels
# subfigs[1].legend(handlesWRZ, labelsWRZ, loc='center left', bbox_to_anchor=(0.9, 0.5))
#fig.tight_layout()
#plt.gca().set_axis_off()
plt.subplots_adjust(left=0.13, right=0.7, top=0.87, bottom=0.205, hspace = 2, wspace = 0.1)
#plt.margins(0,0)
#plt.gca().xaxis.set_major_locator(plt.NullLocator())
#plt.gca().yaxis.set_major_locator(plt.NullLocator())
fig.savefig('growthPlots.png',
    pad_inches = 0.05)