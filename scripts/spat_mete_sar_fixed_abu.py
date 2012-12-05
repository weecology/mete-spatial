"""
Author: Dan McGlinn
Date: 9/21/11
Purpose: to calculate the METE prediction for the SAR when abundance is fixed
"""

import numpy as np
import csv
import sys
import os

import mete

if len(sys.argv) > 1:
    site_names = [sys.argv[1]]
else:
    site_names = ['bci','cocoli1','cocoli2','cross','sherman1','sherman2',
                  'sherman3', 'serp', 'oosting', 'ferp', 'luquillo', 'graveyard',
                  'landsend', 'rocky', 'bormann', 'woodbridge', 'baldmnt', 'bryan',
                  'bigoak']

for shrt_name in site_names:
    datafile = open('../sar/' + shrt_name + '_empir_sar.csv', 'r')
    datareader = csv.reader(datafile)
    data = []
    for row in datareader:
        data.append(row)
   
    # drop the header row
    site_data = []
    for i in range(1, len(data)):
        site_data.append(data[i][0:3])

    site_data = [map(float, x) for x in site_data]

    site_data = np.array(site_data)

    # read in SAD information
    sad_file = open('../data/' + shrt_name + '_sad.csv', 'r')
    sadreader = csv.reader(sad_file)
    sad = []
    for row in sadreader:
        sad.append(row)
        
    # get parameters needed for computing the mete sar
    Amin = min(site_data[ : , 0])
    Amax = max(site_data[ : , 0])    
    n0vals = [int(n0) for n0 in sad[0]]
    
    sar_down_iterative = mete.downscale_sar_fixed_abu(Amax, n0vals, Amin)
    Avals = sar_down_iterative[0][ : ]
    
    sar_down_noniterative = mete.sar_noniterative_fixed_abu(Avals, Amax, n0vals)

    # add values at Amax
    sar_down_iterative[0].append(Amax)
    sar_down_iterative[1].append(len(n0vals))
    
    # Make an array so that the data is easier to output
    out = np.empty((len(sar_down_iterative[0]), 3)) 
    for i in range(0, 2):
        out[ : , i] = sar_down_iterative[i] 
    
    out[ : , 2] = sar_down_noniterative[1]
    
    filename = '../sar/' + shrt_name + '_empirSAD_mete_sar.txt'
    writer = open(filename, 'wb') 
    datawriter = csv.writer(writer)
    datawriter.writerow(['area', 'sr_iter', 'sr_noniter'])
    for i in range(0, np.shape(out)[0]):
        datawriter.writerow(out[i, ])
    
    writer.close()

