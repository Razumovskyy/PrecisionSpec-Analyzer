import numpy as np
import matplotlib.pyplot as plt

wavenumber, absorptionCoeff = np.loadtxt(fname='./ObtainedData.dat', unpack=True, delimiter=',')
refWavenumber, refAbsorptionCoeff = np.loadtxt(fname='./ReferenceData.dat', unpack=True, delimiter=' ')

fig, ax = plt.subplots()

ax.plot(wavenumber, absorptionCoeff, color='b', label='Obtained data')
ax.plot(refWavenumber, refAbsorptionCoeff, color='r', label='Reference data')

ax.set_xlabel(r'wavenumber, $\mathregular{cm^{-1}}$')
ax.set_ylabel(r'$k,\,\mathregular{km^{-1}}$')
ax.grid(which='major', axis='both', color='gray', alpha=0.5)

ax.legend()

# Uncomment for setting X-Y limits on the plot
# ax.set_xlim()
# ax.set_ylim(bottom=1.1, top=1.4)

# Uncomment for setting a title
# plt.title(r'(species) Absorption Coefficient in $\mathregular{km^{-1}}$ at (pressure value) in the (range)' 
#             r'$\mathregular{cm^{-1}}$ range for (line shape)')

plt.tight_layout()
plt.show()