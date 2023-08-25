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

plotDir = '3_GrowthRandom/'
# print(len(sys.argv))
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
sT = dict()
t2 = dict()

cm = 1/2.54  #centimeters in inches
matplotlib.rcParams.update({'font.size': 7.370078})
fig = plt.figure(figsize=(17*cm, 6.5*cm))
subfigs = fig.subfigures(2, 1, hspace=1)

axsBCT = subfigs[0].subplots(1, 5)
subfigs[0].suptitle('BCT')

axsWRZ = subfigs[1].subplots(1, 5)
subfigs[1].suptitle('WRZ')

cB = 0
cW = 0
lims = [100, 200, 300, 400, 500]
tempsBCT = []
tempsWRZ = []
for d in directories:
    if 'BCT' in d:
        ax = axsBCT[cB]
        indLim = cB
        cB += 1
    elif 'WRZ' in d:
        ax = axsWRZ[cW]
        indLim = cW
        cW += 1
    #Get directory of corresponding relaxation
    dirRelax = d.replace(plotDir, '2_RelaxDnvt/')
    crystal, inserted = re.findall(r"\w+N", dirRelax)[0].split('_')
    temp = '_1350K'
    if int(inserted.replace('N', '')) < 100:
        temp = '_1250K'
    dirRelax = dirRelax.replace(crystal+'_'+inserted, crystal+temp)
    dirRelax = dirRelax + 'N_'+inserted[:-1]+'/'
    nC[dirRelax] = []
    if not os.path.exists(dirRelax+'step'+str(1)+'/dump'+str(0)+'.PROB.trj'):
        continue
    #Plotting relaxation part
    t = []
    for i in range(1, 3):
        for b in range(9):
            if i == 1:
                t.append(b*1000)
            elif i == 2:
                t.append(b*1000 + 8000)
            
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
    ax.set_xlabel('$t$ / ps')
    if indLim == 0:
        ax.set_ylabel('$N$')
    ax.axvline(8, ls='-.', color='black', alpha=0.3)
    ax.axvline(16, ls='-.', color='black', alpha=0.3)
    trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)
    #plt.text(8-7, 0.99, 'Step 1', va='top', transform=trans)
    #plt.text(16-7, 0.99, 'Step 2', va='top', transform=trans)
    ax.plot(t, nC[dirRelax])
    ax.axhline(nC[dirRelax][-1], 0, 50, ls='--', color='black')
    #plt.text(0, nC[dirRelax][-1]+25, '$N_c = '+str(nC[dirRelax][-1])+'$', ha='left')

    #Plotting NVT simulations at different temperatures
    dirTemp = glob.glob(d+'T_*K/')
    dirTemp.sort(key=natural_keys)
    for dT in dirTemp:
        print(dT)
        dirSeeds = glob.glob(dT+'seed*/')
        dirSeeds.sort(key=natural_keys)
        nT[dT] = []
        sT[dT] = []
        t2[dT] = []
        # Get temperature values
        files = glob.glob(dirSeeds[0]+'dump*.PROB.trj')
        files.sort(key=natural_keys)
        for f in files:
            t2[dT].append(int(os.path.basename(f).replace('dump', '').replace('.PROB.trj', ''))/1000+16)
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
        nT[dT] = np.mean(nSeed, axis=0)
        sT[dT] = np.std(nSeed, axis=0)
        temp = os.path.basename(dT[:-1]).replace('T_', '')
        indexColor = 0
        if 'BCT' in d:
            if temp not in tempsBCT:
                tempsBCT.append(temp)
            for t in range(len(tempsBCT)):
                if tempsBCT[t] == temp:
                    indexColor = t
        elif 'WRZ' in d:
            if temp not in tempsWRZ:
                tempsWRZ.append(temp)
            for t in range(len(tempsWRZ)):
                if tempsWRZ[t] == temp:
                    indexColor = t
        line = ax.plot(t2[dT], nT[dT], label=temp, c='C'+str(indexColor+1))
        ax.fill_between(t2[dT], nT[dT]+sT[dT], nT[dT]-sT[dT], alpha=0.25, color=line[0]._color)
    ax.set_ylim([0, lims[indLim]])

handlesBCT = []
labelsBCT = []
handlesWRZ = []
labelsWRZ = []
for ax in axsBCT:
    handles, labels = ax.get_legend_handles_labels()
    indRemove = []
    for l in range(len(labels)):
        if labels[l] in labelsBCT:
            indRemove.append(l)
    handles = [handles[i] for i in range(len(handles)) if i not in indRemove]
    labels = [labels[i] for i in range(len(labels)) if i not in indRemove]
    handlesBCT += handles
    labelsBCT += labels
subfigs[0].legend(handlesBCT, labelsBCT, loc='center left', bbox_to_anchor=(0.9, 0.5))
for ax in axsWRZ:
    handles, labels = ax.get_legend_handles_labels()
    indRemove = []
    for l in range(len(labels)):
        if labels[l] in labelsWRZ:
            indRemove.append(l)
    handles = [handles[i] for i in range(len(handles)) if i not in indRemove]
    labels = [labels[i] for i in range(len(labels)) if i not in indRemove]
    handlesWRZ += handles
    labelsWRZ += labels
subfigs[1].legend(handlesWRZ, labelsWRZ, loc='center left', bbox_to_anchor=(0.9, 0.5))
#fig.tight_layout()
#plt.gca().set_axis_off()
plt.subplots_adjust(hspace = 2, wspace = 0.35)
#plt.margins(0,0)
#plt.gca().xaxis.set_major_locator(plt.NullLocator())
#plt.gca().yaxis.set_major_locator(plt.NullLocator())
fig.savefig('growthPlots.pdf',
    pad_inches = 0.05)