import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import os.path
import sys

directories = ['outputAdaptive', 'outputOriginal', 'outputVoronoi']
nC = dict()
t = []
for d in directories:
    nC[d] = []

for i in range(1, 3):
    for b in range(9):
        if i == 1:
            t.append(b*1000)
        elif i == 2:
            t.append(b*1000 + 8000)
        for d in directories:
            header = ''
            with open(d+'/Relaxation750K/step'+str(i)+'_750K/dump'+str(b*1000)+'.PROB.trj') as file:
                for item in file:
                    if 'ITEM: ATOMS' in item:
                        header = item.split(' ')
                        break
            index = 0
            for s in range(len(header)):
                if 'pLIQ' in header[s]:
                    index = s - 2
                    break

            a = np.genfromtxt(d+'/Relaxation750K/step'+str(i)+'_750K/dump'+str(b*1000)+'.PROB.trj', skip_header=9)
            nC[d].append(sum(x < 0.5 for x in a[:, index]))
t = np.array(t)

plt.figure()
plt.title('Cluster size Non-liquid 750K')
plt.xlabel('t / ps')
plt.ylabel('N')
plt.axvline(8000, ls='--', color='black')
plt.axvline(16000, ls='--', color='black')
for l in nC:
    plt.plot(t, nC[l], label=l.replace('output', ''))
plt.legend()
plt.savefig('relaxation750K.png')