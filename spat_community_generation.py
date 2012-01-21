"""
Author: Dan McGlinn
Date: 9/21/11
Purpose: to generate spatial community data using the random bisection
scheme and the MaxEnt (i.e., log-series) SAD. The output is a comma
delimited .txt file site x species matrix the first three rows specify
community #, quadrat x coordinate, and quadrat y coordinate.  This is
denoted in the file's header row. Note: this command can be run from the
command line specifying the following argument: number of species,
number of individuals, number of communities and number of bisections.
The file is output into the directory where this script is called from
by default.

"""

import numpy as np
import csv
import sys

import mete
import spat

if(len(sys.argv) > 1):
    S = int(sys.argv[1]) 
    N = int(sys.argv[2]) 
    ncomm = int(sys.argv[3]) 
    bisec = int(sys.argv[4])
    transect = str2bool(sys.argv[5])
    abu = sys.argv[6]
    shrt_name = sys.argv[7]
else:
    S = 10
    N = 100
    ncomm = 1
    bisec = 9
    transect = False
    abu = 'None'
    shrt_name = None
    
if(abu != 'None'):
    datafile = open(abu,'r')
    datareader = csv.reader(datafile)
    data = []
    for row in datareader:
        data.append(row)
    abu = [int(x) for x in data[0]]
else:
    abu = None
 
nquad = 2 ** (bisec - 1) # number of quadrats per community 
 
comms = [mete.sim_spatial_whole(S, N, bisec, transect, abu) for i in range(0, ncomm)]

# Make an array so that the data is easier to output
out = np.empty((ncomm, nquad, S + 3)) 
for i in range(0,ncomm):
    for j in range(0, nquad):
        out[i,j,] = [i + 1] + comms[i][j][0:2] + comms[i][j][2]

# create a data header
header = []
for i in range(1, S+1):
    header.append('sp' + str(i))

header = ['comm','x','y'] + header

filename = spat.comm_filename(S,N,ncomm,bisec,transect,abu,shrt_name)

writer = open(filename,'wb') 
datawriter = csv.writer(writer)
datawriter.writerow(header)
for i in range(0, ncomm):
    datawriter.writerows(out[i,:,:])

writer.close()  

