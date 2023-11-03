import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

#sizeBCT = np.array([32, 66, 92, 175, 348])
sizeWRZ = np.array([500, 1000, 1500, 2000])

#critBCT = np.array([1304.6875, 1367.1875, 1398.4375, 1414.0625, 1429.6875])
#errBCT = np.array([7.8125, 7.8125, 7.8125, 7.8125, 7.8125])

# brute force ratio: 0.653146541
meltWRZ = np.array([1212.528782538, 1401.342664, 1503.49301192, 1531.0499834])
#errWRZ = np.array([7.8125, 7.8125, 7.8125, 7.8125, 7.8125])

cm = 1/2.54  #centimeters in inches
plt.figure(figsize=(6.2*cm, 4.75*cm))
matplotlib.rcParams.update({'font.size': 9.98655939})
plt.title('Melting temperatures', fontsize=9.98655939)
plt.xlabel('$N$', labelpad=1)
plt.ylabel(r'$T_{melt}$ [K]', labelpad=1)

#line = plt.errorbar(sizeBCT, critBCT, yerr=errBCT, label='BCT', lw=1, capsize=5.0)
#plt.fill_between(sizeBCT, critBCT+errBCT, critBCT-errBCT, alpha=0.25, color=line[0]._color)

line = plt.plot(sizeWRZ, meltWRZ, label='WRZ', lw=1, marker='o')
#plt.fill_between(sizeWRZ, critWRZ+errWRZ, critWRZ-errWRZ, alpha=0.25, color=line[0]._color)

plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), handlelength=1)
plt.tight_layout()
plt.savefig('melttemps.png', pad_inches = 0.0)