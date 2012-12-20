"""
Purpose: to compute the chi_heap recursion equation of Harte 2007 
and to export the results as a .csv file
"""

import mete
import numpy as np
import csv
import sys

if len(sys.argv) > 1:
    A = int(sys.argv[1]) 
    A0 = int(sys.argv[2])
    shape = sys.argv[3]
    abu_filename = sys.argv[4]
    out_filename = sys.argv[5]
else:
    print 'Error: Specify A, A0, shape, and abu file at the command line'

if abu_filename  != 'None':
    datafile = open(abu_filename, 'r')
    datareader = csv.reader(datafile)
    data = []
    for row in datareader:
        data.append(row)
    abu = [int(x) for x in data[0]]
else:
    print 'Error: Abundance data not specified'

sor = mete.sor_heap(A, abu, A0, shape)

writer = open(out_filename, 'wb') 
datawriter = csv.writer(writer)
 
datawriter.writerow(['j', 'dist', 'sor'])
for i in range(0, np.shape(sor)[0]):
    datawriter.writerow(sor[i, ])
     
writer.close()
