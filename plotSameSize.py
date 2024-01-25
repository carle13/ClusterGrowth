import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# sizeBCT2000 = np.array([32, 66, 92, 175, 348])
# sizeWRZ2000 = np.array([33, 64, 93, 164, 338])

# critBCT2000 = np.array([1304.6875, 1367.1875, 1398.4375, 1414.0625, 1429.6875])
# errBCT2000 = np.array([7.8125, 7.8125, 7.8125, 7.8125, 7.8125])

# critWRZ2000 = np.array([1367.1875, 1382.8125, 1429.6875, 1445.3125, 1460.9375])
# errWRZ2000 = np.array([7.8125, 7.8125, 7.8125, 7.8125, 7.8125])
# ################
# sizeBCT1500 = np.array([29, 88, 136])
# sizeWRZ1500 = np.array([32, 83, 133])

# critBCT1500 = np.array([1304.6875, 1351.5625, 1398.4375])
# errBCT1500 = np.array([7.8125, 7.8125, 7.8125])

# critWRZ1500 = np.array([1335.9375, 1382.8125, 1414.0625])
# errWRZ1500 = np.array([7.8125, 7.8125, 7.8125])
# ################
# sizeBCT1000 = np.array([27, 52, 87])
# sizeWRZ1000 = np.array([24, 51, 86])

# critBCT1000 = np.array([1242.1875, 1289.0625, 1257.8125])
# errBCT1000 = np.array([7.8125, 7.8125, 7.8125])

# critWRZ1000 = np.array([1289.0625, 1273.4375, 1335.9375])
# errWRZ1000 = np.array([7.8125, 7.8125, 7.8125])
# ################
# sizeBCT500 = np.array([22, 35, 40])
# sizeWRZ500 = np.array([22, 32, 38])

# critBCT500 = np.array([1039.0625, 1101.5625, 1117.1875])
# errBCT500 = np.array([15.625, 15.625, 15.625])

# critWRZ500 = np.array([1023.4375, 1117.1875, 1101.5625])
# errWRZ500 = np.array([15.625, 15.625, 15.625])


#Critical clusters used: 500atoms - 35,32,   1000atoms - 27,24,    1500atoms - 29,32,   2000atoms - 32,33
dropletSizes = [500, 1000, 1500, 2000]
critTempsWRZ = [1117.1875, 1289.0625, 1335.9375, 1367.1875]
critTempsBCT = [1101.5625, 1242.1875, 1304.6875, 1304.6875]
errorTempsWRZ = [15.625, 7.8125, 7.8125, 7.8125]
errorTempsBCT = [15.625, 7.8125, 7.8125, 7.8125]


cm = 1/2.54  #centimeters in inches
plt.figure(figsize=(6.6*cm, 6*cm))
matplotlib.rcParams.update({'font.size': 9.98655939})
#plt.title('Critical temperatures', fontsize=9.98655939)
plt.xlabel('$N$', labelpad=1)
plt.ylabel(r'$T_{crit}$ [K]', labelpad=1)

line = plt.errorbar(dropletSizes, critTempsBCT, yerr=errorTempsBCT, label='BCT', lw=1, capsize=2.0, color=(228/255,26/255,28/255))
#plt.fill_between(sizeBCT, critBCT+errBCT, critBCT-errBCT, alpha=0.25, color=line[0]._color)
line = plt.errorbar(dropletSizes, critTempsWRZ, yerr=errorTempsWRZ, label='WRZ', lw=1, capsize=2.0, color=(77/255,175/255,74/255))
#plt.fill_between(sizeWRZ, critWRZ+errWRZ, critWRZ-errWRZ, alpha=0.25, color=line[0]._color)

#plt.semilogx()
#plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), handlelength=1)
plt.legend()
plt.tight_layout()
plt.savefig('critSameSize.pdf', pad_inches = 0.0)