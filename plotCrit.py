import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sizeBCT = np.array([72, 177, 364])
sizeWRZ = np.array([71, 186, 366])

critBCT = np.array([1367.1875, 1421.875, 1445.3125])
errBCT = np.array([7.8125, 15.625, 7.8125])

critWRZ = np.array([1328.125, 1437.5, 1445.3125])
errWRZ = np.array([15.625, 15.625, 7.8125])

cm = 1/2.54  #centimeters in inches
plt.figure(figsize=(8.6*cm, 8.6*cm))
matplotlib.rcParams.update({'font.size': 7.370078})
#plt.title('Critical temperatures')
plt.xlabel('$N_c$')
plt.ylabel(r'$T_{crit}$ / K')

line = plt.errorbar(sizeBCT, critBCT, yerr=errBCT, label='BCT', capsize=5.0)
#plt.fill_between(sizeBCT, critBCT+errBCT, critBCT-errBCT, alpha=0.25, color=line[0]._color)

line = plt.errorbar(sizeWRZ, critWRZ, yerr=errWRZ, label='WRZ', capsize=5.0)
#plt.fill_between(sizeWRZ, critWRZ+errWRZ, critWRZ-errWRZ, alpha=0.25, color=line[0]._color)

plt.legend()
plt.tight_layout()
plt.savefig('critTemps.pdf', bbox_inches='tight')