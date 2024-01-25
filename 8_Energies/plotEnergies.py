import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

matplotlib.rcParams['axes.unicode_minus'] = False

sizeBCT = np.array([500, 1000, 1500, 2000])
sizeWRZ = np.array([500, 1000, 1500, 2000])
sizeBCT2 = np.array([500, 1000, 2000])

growthBCT = np.array([-1885.15664288/500, -3813.57488895/1000, -7685.71398441/2000])
#seeds    1,     4,     ,     3
#errBCT = np.array([7.8125, 7.8125, 7.8125, 7.8125, 7.8125])

growthWRZ = np.array([-1887.62511173/500, -3788.86114318/1000, -5726.52260482/1500, -7704.63536491/2000])
#seeds   3,   3,    4,     4
#errWRZ = np.array([7.8125, 7.8125, 7.8125, 7.8125, 7.8125])

perfectBCT = np.array([-1861.86549334/500, -3780.31740223/1000, -5711.82271648/1500, -7666.16883946/2000])
perfectWRZ = np.array([-1866.83599912/500, -3812.44812559/1000, -5720.19247942/1500, -7688.46147222/2000])

cm = 1/2.54  #centimeters in inches
plt.figure(figsize=(8.6*cm, 7.25*cm))
matplotlib.rcParams.update({'font.size': 9.98655939})
#plt.title('Nanoparticle Energy', fontsize=9.98655939)
plt.xlabel('$N$', labelpad=1)
plt.ylabel(r'$E$ [eV]', labelpad=1)

line = plt.plot(sizeBCT2, growthBCT, '--o', label='BCT Seed', lw=1, color=(203/255,24/255,29/255))
#plt.fill_between(sizeBCT, critBCT+errBCT, critBCT-errBCT, alpha=0.25, color=line[0]._color)

line = plt.plot(sizeWRZ, growthWRZ, '--o', label='WRZ Seed', lw=1, color=(35/255,139/255,69/255))
#plt.fill_between(sizeWRZ, critWRZ+errWRZ, critWRZ-errWRZ, alpha=0.25, color=line[0]._color)

line = plt.plot(sizeBCT, perfectBCT, '--o', label='BCT Perf', lw=1, color=(251/255,106/255,74/255))
line = plt.plot(sizeWRZ, perfectWRZ, '--o', label='WRZ Perf', lw=1, color=(116/255,196/255,118/255))

#plt.semilogy()
#plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), handlelength=1)
plt.legend()
plt.tight_layout()
plt.savefig('energies.pdf', pad_inches = 0.0)