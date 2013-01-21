"""
Purpose: to recreate the univariate pdfs of Conlisk et al. (2007)
in figure 2. This provides a test that the univariate pdfs are 
working correctly
"""
import numpy as np
import mete
import matplotlib.pyplot as plt

n0 = 617
c = 256

sing_pdfs = np.zeros((4, 8))

psi = [0.01, 0.25, 0.5, 0.75]
for i in range(0, len(psi)):
    sing_pdfs[i, : ] = [mete.single_prob(n, n0, psi[i], c) for n in range(0, 8)]

n = range(0, 8)

plt.plot(n, sing_pdfs[0, :], color='black', linewidth=1)
plt.plot(n, sing_pdfs[1, :], color='black', linewidth=1, ls='--')
plt.plot(n, sing_pdfs[2, :], color='black', linewidth=2)
plt.plot(n, sing_pdfs[3, :], color='lightgrey', linewidth=2)
plt.axis([0, 7, 0, 0.9])
plt.xlabel('n')
plt.ylabel('P(n)')
plt.legend(['psi = 0.01', 'psi = 0.25', 'psi = 0.50','psi = 0.75'])
plt.savefig('../figs/conlisk_univaritate_fig2a.pdf')

