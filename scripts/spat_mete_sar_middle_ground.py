"""
Author: Dan McGlinn
Date: 2/11/13
Purpose: to calculate the recursive downscaling of Pi but not psi to predict the
SAR, this is a middle ground model between the recursive and non-recursive models. Uses only the METE predicted SADS.
"""

import numpy as np
import csv
import sys
import os

import mete

print 'Computing METE SAR, ...'

datafile = open('./data/shrtnames.txt', 'r')
datareader = csv.reader(datafile)
shrtnames = []
for row in datareader:
    shrtnames.append(row)

shrtnames = shrtnames[0][0].split()

if len(sys.argv) > 1:
    site_index = int(sys.argv[1])
    shrtnames = [shrtnames[site_index]]

# set number of permutations to compute results for
nperm = 200

for shrt_name in shrtnames:
    
    datafile = open('./sar/' + shrt_name + '_empir_sar.csv', 'r')
    datareader = csv.reader(datafile)
    data = []
    for row in datareader:
        data.append(row)
    
    # drop the header row
    site_data = []
    for i in range(1, len(data)):
        site_data.append(data[i][0:3])
    
    # convert strings to floats            
    site_data = [map(float, x) for x in site_data]
    
    # enforce a minimum individual density of 2
    indices = mete.which([site_data[i][2] > 2 for i in range(0, len(site_data))])
    site_data = [site_data[i] for i in indices]
    
    site_data = np.array(site_data)
    
    Amin = min(site_data[ : , 0])
    Amax = max(site_data[ : , 0])
    S0 = int(max(site_data[ : , 1]))
    N0 = int(max(site_data[ : , 2]))
    
    sar_down_iterative = []
    for i in range(0, nperm):
        p = mete.exp(-mete.get_beta(S0, N0))
        n0_rvs = mete.trunc_logser_rvs(p, N0, S0)
        sar_down_iterative.append(mete.downscale_sar_fixed_abu(Amax, n0_rvs, Amin))
    
    Avals = sar_down_iterative[0][0][ : ]
    len_A = len(Avals)
    
    out = np.empty((nperm * len_A, 2))
    
    for j in range(0, nperm):
        for i in range(0, 2):
            start = (j * len_A)
            stop = start + len_A
            out[start : stop, i] = sar_down_iterative[j][i] 
        
    filename = './sar/' + shrt_name + '_mete_sar_middle_ground.txt'
    
    writer = open(filename, 'wb') 
    datawriter = csv.writer(writer)
    datawriter.writerow(['area', 'sr'])
    for i in range(0, np.shape(out)[0]):
        datawriter.writerow(out[i, ])
    
    writer.close()

print 'Computing METE SAR, complete!'