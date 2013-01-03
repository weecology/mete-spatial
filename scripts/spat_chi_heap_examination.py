"""
Purpose: to recreate the values from Harte 2007 in Fig 6.7 to enable
checking of the chi_heap recursion equations
"""

import mete
import numpy as np
import csv


def chi_heap_approx(i, j, n0):
    """ CHI approximation formula, Harte 2007, Eq. 6.12"""
    return mete.get_lambda_heap(i, n0) ** 2 / mete.get_lambda_heap(j, n0)

out = np.empty((50, 4))
chi = []
jvals = []
ivals = []

irow = 0
for j in range(1, 6):
    for i in range(j + 3, 16):
        out[irow, 0] = i
        out[irow, 1] = j
        out[irow, 2] = mete.chi_heap(i, j, 100)
        out[irow, 3] = chi_heap_approx(i, j, 100)
        irow += 1

filename = '../sorensen/harte_2007_chi_heap_results.txt'
writer = open(filename, 'wb') 
datawriter = csv.writer(writer)
 
datawriter.writerow(['i', 'j', 'chi', 'chi_appr'])
for i in range(0, np.shape(out)[0]):
    datawriter.writerow(out[i, ])
     
writer.close()
