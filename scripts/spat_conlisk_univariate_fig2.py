"""
Purpose: to recreate the univariate pdfs of Conlisk et al. (2007)
in figure 2. This provides a test that the univariate pdfs are 
working correctly
"""
import numpy as np
import mete
import csv

n0 = 617
c = 256

sing_pdfs = np.zeros((4, 7))

psi = [0.01, 0.25, 0.5, 0.75]
for i in range(0, len(psi)):
    sing_pdfs[i, : ] = [mete.single_prob(n, n0, psi[i], c) for n in range(0, 7)]

writer = open('./data/conlisk_data_fig2a.csv', 'wb') 
datawriter = csv.writer(writer)
     
for i in range(0, np.shape(sing_pdfs)[0]):
    datawriter.writerow(sing_pdfs[i, ])

writer.close()
