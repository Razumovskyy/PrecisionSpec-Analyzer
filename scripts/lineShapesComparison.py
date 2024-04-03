import numpy as np
import matplotlib.pyplot as plt

wavenumber, lineShape1 = np.loadtxt(fname='./lineShape2Data.dat', unpack=True, delimiter=',')
refWavenumber, lineShape2 = np.loadtxt(fname='./lineShape1Data.dat', unpack=True, delimiter=' ')

# wavenumber, lineShape1, lineShape2 = np.loadtxt(fname='./lineShapesData.dat', unpack=True, delimiter=',')

fig, ax = plt.subplots()

small_value = 1e-10
lineShape2[lineShape2 < small_value] = small_value

ax.plot(wavenumber, np.log10(lineShape1), color='r', label='LOR')
ax.plot(refWavenumber, np.log10(lineShape2), color='b', label='DOP')

ax.set_xlabel(r'wavenumber, $\mathregular{cm^{-1}}$')
ax.set_ylabel(r'$l,\,\mathregular{cm}$')
ax.grid(which='major', axis='both', color='gray', alpha=0.5)

# Uncomment for setting a title
# plt.title()

plt.tight_layout()
plt.show()