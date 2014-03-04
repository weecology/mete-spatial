"""
Purpose: to compute the chi_heap recursion equation of Harte 2007 
and to export the results as a .csv file
"""

print 'Compute the analytical HEAP/METE DDR, ...'

import mete
import numpy as np
import csv
import sys

if len(sys.argv) > 1:
    bisect_fine = int(sys.argv[1]) 
    bisect_coarse = int(sys.argv[2])
    sadType = sys.argv[3]
    abu_filepath = sys.argv[4]
    out_filepath = sys.argv[5]
    if len(sys.argv) == 7:
        unit_distance = float(sys.argv[6])
    else:
        unit_distance = 1
else:
    print 'Error: Specify bisect_fine, bisect_coarse, sadType,  abu file path, and out filepath at the command line'

if abu_filepath  != 'None':
    datafile = open(abu_filepath, 'r')
    datareader = csv.reader(datafile)
    data = []
    for row in datareader:
        data.append(row)
    abu = [int(x) for x in data[0]]
else:
    print 'Error: Abundance data not specified'

bisections = range(bisect_fine, bisect_coarse - 2, -2)
A0 = 2 ** bisect_fine
A_vals = [A0 / 2 ** i for i in bisections]
shape = 'golden' # so that results for every j value are calculated

if sadType == 'meteSAD':
    S0 = len(abu)
    N0 = sum(abu)
    sor = [mete.sor_heap(A, A0, S0, N0, shape, unit_distance) for A in A_vals]
    sor = np.concatenate((sor), axis=0)
else:
    sor = [mete.sor_heap_fixed_abu(A, abu, A0, shape, unit_distance) for A in A_vals]
    sor = np.concatenate((sor), axis=0)

writer = open(out_filepath, 'wb') 
datawriter = csv.writer(writer)
 
datawriter.writerow(['i', 'j', 'dist', 'sor'])
for i in range(0, np.shape(sor)[0]):
    datawriter.writerow(sor[i, ])
     
writer.close()

print 'Compute the analytical HEAP/METE DDR, complete!'

