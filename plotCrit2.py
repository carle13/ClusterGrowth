import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sizeBCT = np.array([32, 66, 92, 175, 348])
sizeWRZ = np.array([33, 64, 93, 164, 338])

critBCT = np.array([1304.6875, 1367.1875, 1398.4375, 1414.0625, 1429.6875])
errBCT = np.array([7.8125, 7.8125, 7.8125, 7.8125, 7.8125])

critWRZ = np.array([1367.1875, 1382.8125, 1429.6875, 1445.3125, 1460.9375])
errWRZ = np.array([7.8125, 7.8125, 7.8125, 7.8125, 7.8125])

cm = 1/2.54  #centimeters in inches
plt.figure(figsize=(6*cm, 4.63*cm))
matplotlib.rcParams.update({'font.size': 7.370078})
plt.title('Critical temperatures')
plt.xlabel('$N_c$')
plt.ylabel(r'$T_{crit}$ / K')

line = plt.errorbar(sizeBCT, critBCT, yerr=errBCT, label='BCT', capsize=5.0)
#plt.fill_between(sizeBCT, critBCT+errBCT, critBCT-errBCT, alpha=0.25, color=line[0]._color)

line = plt.errorbar(sizeWRZ, critWRZ, yerr=errWRZ, label='WRZ', capsize=5.0)
#plt.fill_between(sizeWRZ, critWRZ+errWRZ, critWRZ-errWRZ, alpha=0.25, color=line[0]._color)

plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.tight_layout()
plt.savefig('critTemps2.pdf', pad_inches = 0.0)