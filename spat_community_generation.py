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
import os

import mete

from itertools import repeat
from math import exp

def str2bool(v):
    return v.lower() in ("yes", "true", "t", "1")

def comm_filename(S, N, ncomm, bisec, transect=False, abu=None, comm_name=None):
    """Get filename for simulated community produced by mete.sim_spatial_whole()
    
    Keyword arguments:
    S -- the number of species
    N -- the number of individuals
    ncomm -- the number of communities that were generated
    bisec -- the number of bisections
    transect -- a boolean if False then it is assumed a 2D grid was generated
    abu -- the path to an abundance file else this should be None
    comm_name -- name for the output community
    """
                 
    if not comm_name:
        comm_name = 'S%s_N%s' % (S, N)
    if abu:
        empir = '_empirSAD'
    else:
        empir = ''
    if transect:
        runtype = 'transect'
    else:
        runtype = 'grid'
    return './comms/simulated_comms_%s%s_C%s_B%s_%s.txt' % (comm_name,
                                                                empir, ncomm,
                                                                bisec, runtype)

def output_results(comms, S, N, ncomm, bisec, transect, abu, shrt_name):
    
    nquad = 2 ** (bisec - 1) # number of quadrats per community 
    if 'comms' not in os.listdir(os.path.curdir):
            os.mkdir('comms')    

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
    filename = comm_filename(S,N,ncomm,bisec,transect,abu,shrt_name)
    writer = open(filename,'wb') 
    datawriter = csv.writer(writer)
    datawriter.writerow(header)
    for i in range(0, ncomm):
        datawriter.writerows(out[i,:,:])
    writer.close()
    
def explore_parameter_space(Svals, Nvals, ncomm, bisec, transect=False):
    for S in Svals:
        for N in Nvals:
            beta = mete.get_beta(S, N)
            if exp(-beta) < 1:
                comms = [mete.sim_spatial_whole(S, N, bisec, transect=transect,
                                            beta=beta) for i in range(0, ncomm)]
                output_results(comms, S, N, ncomm, bisec, transect, None, None)
            print S, N

if __name__ == "__main__":
        
    if len(sys.argv) > 1:
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
        
    if abu != 'None':
        datafile = open(abu,'r')
        datareader = csv.reader(datafile)
        data = []
        for row in datareader:
            data.append(row)
        abu = [int(x) for x in data[0]]
    else:
        abu = None
     
    comms = [mete.sim_spatial_whole(S, N, bisec, transect, abu) for i in range(0, ncomm)]
    output_results(comms, S, N, ncomm, bisec, transect, abu, shrt_name)