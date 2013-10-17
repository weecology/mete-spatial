"""
Author: Dan McGlinn
Date: 9/21/11
Purpose: to calculate the METE prediction for the SAR
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


for shrt_name in shrtnames:
    for sadType in ['meteSAD', 'empirSAD']:
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
        
        if sadType == 'meteSAD':
            # enforce a minimum individual density of 2
            indices = mete.which([site_data[i][2] > 2 for i in range(0, len(site_data))])
            site_data = [site_data[i] for i in indices]
        
        site_data = np.array(site_data)

        # get parameters needed for computing the mete sar
        if sadType == 'empirSAD':
            # read in SAD information  
            sad_file = open('./data/' + shrt_name + '_sad.csv', 'r')
            sadreader = csv.reader(sad_file)
            sad = []
            for row in sadreader:
                sad.append(row)

        Amin = min(site_data[ : , 0])
        Amax = max(site_data[ : , 0])
        if sadType == 'meteSAD':
            S0 = int(max(site_data[ : , 1]))
            N0 = int(max(site_data[ : , 2]))
        else:
            n0vals = [int(n0) for n0 in sad[0]]

        if sadType == 'meteSAD':
            sar_down_iterative = mete.downscale_sar(Amax, S0, N0, Amin)
        else:
            sar_down_iterative = mete.downscale_sar_fixed_abu(Amax, n0vals, Amin)
            
        Avals = sar_down_iterative[0][ : ]
        
        try:
            if sadType == 'meteSAD':
                sar_down_noniterative = mete.sar_noniterative(Avals, Amax, S0, N0, 'precise')
            else:
                sar_down_noniterative = mete.sar_noniterative_fixed_abu(Avals, Amax, n0vals)
            sar_noniterative_worked = True
        except ValueError:
            print "Downscaling non-iterative SAR failed"
            sar_noniterative_worked = False
    
        if sadType == 'meteSAD':
            # add values at Amax
            sar_down_iterative[0].append(Amax)
            sar_down_iterative[1].append(S0)

        # Make an array so that the data is easier to output
        if sar_noniterative_worked:    
            out = np.empty((len(sar_down_iterative[0]), 3)) 
        else:
            out = np.empty((len(sar_down_iterative[0]), 2))       
        
        for i in range(0, 2):
            out[ : , i] = sar_down_iterative[i] 

        if sar_noniterative_worked:                 
            out[ : , 2] = sar_down_noniterative[1]
              
        if sadType == 'meteSAD':
            filename = './sar/' + shrt_name + '_mete_sar.txt'
        else:
            filename = './sar/' + shrt_name + '_empirSAD_mete_sar.txt'
        
        writer = open(filename, 'wb') 
        datawriter = csv.writer(writer)
    
        if sar_noniterative_worked:
            datawriter.writerow(['area', 'sr_iter', 'sr_noniter'])
        else:
            datawriter.writerow(['area', 'sr_iter'])

        for i in range(0, np.shape(out)[0]):
            datawriter.writerow(out[i, ])
    
        writer.close()

print 'Computing METE SAR, complete!'