import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sizeBCT2000 = np.array([32, 66, 92, 175, 348])
sizeWRZ2000 = np.array([33, 64, 93, 164, 338])

critBCT2000 = np.array([1304.6875, 1367.1875, 1398.4375, 1414.0625, 1429.6875])
errBCT2000 = np.array([7.8125, 7.8125, 7.8125, 7.8125, 7.8125])

critWRZ2000 = np.array([1367.1875, 1382.8125, 1429.6875, 1445.3125, 1460.9375])
errWRZ2000 = np.array([7.8125, 7.8125, 7.8125, 7.8125, 7.8125])
################
sizeBCT1500 = np.array([29, 88, 136])
sizeWRZ1500 = np.array([32, 83, 133])

critBCT1500 = np.array([1304.6875, 1351.5625, 1398.4375])
errBCT1500 = np.array([7.8125, 7.8125, 7.8125])

critWRZ1500 = np.array([1335.9375, 1382.8125, 1414.0625])
errWRZ1500 = np.array([7.8125, 7.8125, 7.8125])
################
sizeBCT1000 = np.array([27, 52, 87])
sizeWRZ1000 = np.array([24, 51, 86])

critBCT1000 = np.array([1242.1875, 1289.0625, 1257.8125])
errBCT1000 = np.array([7.8125, 7.8125, 7.8125])

critWRZ1000 = np.array([1289.0625, 1273.4375, 1335.9375])
errWRZ1000 = np.array([7.8125, 7.8125, 7.8125])
################
sizeBCT500 = np.array([22, 35, 40])
sizeWRZ500 = np.array([22, 32, 38])

critBCT500 = np.array([1039.0625, 1101.5625, 1117.1875])
errBCT500 = np.array([15.625, 15.625, 15.625])

critWRZ500 = np.array([1023.4375, 1117.1875, 1101.5625])
errWRZ500 = np.array([15.625, 15.625, 15.625])

cm = 1/2.54  #centimeters in inches
plt.figure(figsize=(8.6*cm, 6*cm))
matplotlib.rcParams.update({'font.size': 9.98655939})
plt.title('Critical temperatures', fontsize=9.98655939)
plt.xlabel('$N_c$', labelpad=1)
plt.ylabel(r'$T_{crit}$ [K]', labelpad=1)

line = plt.errorbar(sizeBCT2000, critBCT2000, yerr=errBCT2000, label='BCT - 2000', lw=1, capsize=2.0, color=(153/255,0/255,13/255))
#plt.fill_between(sizeBCT, critBCT+errBCT, critBCT-errBCT, alpha=0.25, color=line[0]._color)
line = plt.errorbar(sizeWRZ2000, critWRZ2000, yerr=errWRZ2000, label='WRZ - 2000', lw=1, capsize=2.0, color=(0/255,90/255,50/255))
#plt.fill_between(sizeWRZ, critWRZ+errWRZ, critWRZ-errWRZ, alpha=0.25, color=line[0]._color)

line = plt.errorbar(sizeBCT1500, critBCT1500, yerr=errBCT1500, label='BCT - 1500', lw=1, capsize=2.0, color=(203/255,24/255,29/255))
#plt.fill_between(sizeBCT, critBCT+errBCT, critBCT-errBCT, alpha=0.25, color=line[0]._color)
line = plt.errorbar(sizeWRZ1500, critWRZ1500, yerr=errWRZ1500, label='WRZ - 1500', lw=1, capsize=2.0, color=(35/255,139/255,69/255))
#plt.fill_between(sizeWRZ, critWRZ+errWRZ, critWRZ-errWRZ, alpha=0.25, color=line[0]._color)

line = plt.errorbar(sizeBCT1000, critBCT1000, yerr=errBCT1000, label='BCT - 1000', lw=1, capsize=2.0, color=(239/255,59/255,44/255))
#plt.fill_between(sizeBCT, critBCT+errBCT, critBCT-errBCT, alpha=0.25, color=line[0]._color)
line = plt.errorbar(sizeWRZ1000, critWRZ1000, yerr=errWRZ1000, label='WRZ - 1000', lw=1, capsize=2.0, color=(65/255,171/255,93/255))
#plt.fill_between(sizeWRZ, critWRZ+errWRZ, critWRZ-errWRZ, alpha=0.25, color=line[0]._color)

line = plt.errorbar(sizeBCT500, critBCT500, yerr=errBCT500, label='BCT - 500', lw=1, capsize=2.0, color=(251/255,106/255,74/255))
#plt.fill_between(sizeBCT, critBCT+errBCT, critBCT-errBCT, alpha=0.25, color=line[0]._color)
line = plt.errorbar(sizeWRZ500, critWRZ500, yerr=errWRZ500, label='WRZ - 500', lw=1, capsize=2.0, color=(116/255,196/255,118/255))
#plt.fill_between(sizeWRZ, critWRZ+errWRZ, critWRZ-errWRZ, alpha=0.25, color=line[0]._color)

plt.semilogx()
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), handlelength=1)
plt.tight_layout()
plt.savefig('critNano.pdf', pad_inches = 0.0)