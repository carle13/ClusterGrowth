import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sizeBCT = np.array([500, 1000, 1500, 2000])
sizeWRZ = np.array([500, 1000, 1500, 2000])
sizeBCT2 = np.array([500, 1000, 2000])

growthBCT = np.array([-1885.15664288, -3813.57488895, -7685.71398441]) *-1
#seeds    1,     4,     ,     3
#errBCT = np.array([7.8125, 7.8125, 7.8125, 7.8125, 7.8125])

growthWRZ = np.array([-1887.62511173, -3788.86114318, -5726.52260482, -7704.63536491]) *-1
#seeds   3,   3,    4,     4
#errWRZ = np.array([7.8125, 7.8125, 7.8125, 7.8125, 7.8125])

perfectBCT = np.array([-1861.86549334, -3780.31740223, -5711.82271648, -7666.16883946]) *-1
perfectWRZ = np.array([-1866.83599912, -3812.44812559, -5720.19247942, -7688.46147222]) *-1

cm = 1/2.54  #centimeters in inches
plt.figure(figsize=(12*cm, 8*cm))
matplotlib.rcParams.update({'font.size': 9.98655939})
plt.title('Nanoparticle Energy', fontsize=9.98655939)
plt.xlabel('$N$', labelpad=1)
plt.ylabel(r'$E$ [eV]', labelpad=1)

line = plt.plot(sizeBCT2, growthBCT, '--o', label='BCT Growth', lw=1, c='r')
#plt.fill_between(sizeBCT, critBCT+errBCT, critBCT-errBCT, alpha=0.25, color=line[0]._color)

line = plt.plot(sizeWRZ, growthWRZ, '--o', label='WRZ Growth', lw=1, c='g')
#plt.fill_between(sizeWRZ, critWRZ+errWRZ, critWRZ-errWRZ, alpha=0.25, color=line[0]._color)

line = plt.plot(sizeBCT, perfectBCT, '--o', label='BCT Perfect', lw=1, c='b')
line = plt.plot(sizeWRZ, perfectWRZ, '--o', label='WRZ Perfect', lw=1, c='y')

plt.semilogy()
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), handlelength=1)
plt.tight_layout()
plt.savefig('energies.png', pad_inches = 0.0)