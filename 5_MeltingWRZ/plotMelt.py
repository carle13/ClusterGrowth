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

probs = ['pliq', 'pbct', 'pwrz', 'phbn']
nC = dict()
t = []
for d in probs:
    nC[d] = []

for d in probs:
    t = []
    for b in range(21):
        t.append(b*5000)
        #Count particles belonging to the cluster
        pipeline = import_file('Output/100ps/dump'+str(b*5000)+'.PROB.trj', multiple_frames=True)
        pipeline.modifiers.append(ExpressionSelectionModifier(expression=d+' > 0.5'))
        #pipeline.modifiers.append(ClusterAnalysisModifier(cutoff=4, sort_by_size=True, only_selected=True))
        data = pipeline.compute()
        nC[d].append(data.attributes['ExpressionSelection.count'])
t = np.array(t) / 1000

#Draw plot
plt.figure()
#plt.title('WRZ Melting Simulation (1000K - 3000K)')
plt.xlabel('$t$ / ps')
plt.ylabel('$N$')
# plt.axvline(8, ls='-.', color='black', alpha=0.3)
# plt.axvline(16, ls='-.', color='black', alpha=0.3)
# ax = plt.gca()
# trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)
# plt.text(8-2.5, 0.03, 'Step 1', va='bottom', transform=trans)
# plt.text(16-2.5, 0.03, 'Step 2', va='bottom', transform=trans)
for d in probs:
    plt.plot(t, nC[d], label=d.replace('p', '').upper())
#plt.axhline(nC[d][-1], 0, 16, ls='--', color='black')
#plt.text(16-0.15, nC[d][-1]+1, '$N_c = '+str(nC[d][-1])+'$', ha='right')
#plt.legend()
plt.tight_layout()
#plt.xlim([0, 21])
plt.savefig('wrzmelt100ps.png')
plt.close()