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

from scipy.optimize import curve_fit

probs = ['pbct', 'pwrz', 'phbn']
nC = dict()
t = []
for d in probs:
    nC[d] = []

for s in range(1, 6):
    seed = 11111 * s
    for d in probs:
        vals = []
        for b in range(0,450,2):
            #Count particles belonging to the cluster
            pipeline = import_file('Output/500ps/seed'+str(seed)+'/dump'+str(b*1000)+'.PROB.trj', multiple_frames=True)
            pipeline.modifiers.append(ExpressionSelectionModifier(expression='pbct > 0.5 || pwrz > 0.5 || phbn > 0.5'))
            pipeline.modifiers.append(ClusterAnalysisModifier(cutoff=4, sort_by_size=True, only_selected=True))
            pipeline.modifiers.append(ExpressionSelectionModifier(expression=d+' > 0.5  && Cluster == 1'))
            data = pipeline.compute()
            vals.append(data.attributes['ExpressionSelection.count.2'])
        nC[d].append(vals)

for b in range(0,450,2):
    t.append(b*1000)
t = (np.array(t) / 1000)

#Draw plot
plt.figure()
#plt.title('BCT Nanoparticle Crystallization (1375K)')
plt.xlabel('$t$ / ps')
plt.ylabel('$N$')

nTotal = []
for s in range(len(nC['pbct'])):
    nTotal.append(np.array(nC['pbct'][s]) + np.array(nC['pwrz'][s]) + np.array(nC['phbn'][s]))

for d in probs:
    average = np.mean(nC[d], axis=0)
    standDev = np.std(nC[d], axis=0)
    line = plt.plot(t, average, label=d.replace('p', '').upper())
    plt.fill_between(t, average+standDev, average-standDev, alpha=0.25, color=line[0]._color)

aveTot = np.mean(nTotal, axis=0)
stdTot = np.std(nTotal, axis=0)
line = plt.plot(t, aveTot, label='Total Cluster')
plt.fill_between(t, aveTot+stdTot, aveTot-stdTot, alpha=0.25, color=line[0]._color)
#plt.axhline(nC[d][-1], 0, 16, ls='--', color='black')
#plt.text(16-0.15, nC[d][-1]+1, '$N_c = '+str(nC[d][-1])+'$', ha='right')
plt.legend()
plt.tight_layout()
#plt.ylim([0, 300])
plt.savefig('meltingcomposition.png')
plt.close()


#Fitting of the curve
def func(x, A,  x0, sigma):                                                                                                  
    return A*(1-np.tanh((x-x0)/sigma))

plt.figure()
plt.plot(t, aveTot, label='Total Cluster')
popt, pcov = curve_fit(func, t, aveTot)
print(popt)
print('x0 = '+str(popt[1]))
plt.plot(t, func(t, popt[0], popt[1], popt[2]), label='Fit')
plt.legend()
plt.savefig('meltingfit.png')