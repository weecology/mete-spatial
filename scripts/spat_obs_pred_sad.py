import numpy as np
import csv
import sys
import os

from mete import *

def get_sad_data(site_name, sad_path='../data/'):
    datafile = open(sad_path + site_name + '_sad.csv', 'r')
    datareader = csv.reader(datafile)
    n0 = []
    for row in datareader:
        n0.append(row)
    
    n0 = n0[0]
    # convert strings to floats            
    n0 = map(int, n0)
    return n0

def get_obs_pred_sad(site_name, sad_path='../data/', S_cutoff=1,
                     to_file=False, output_dir='../sad/'):
    if to_file:    
        if os.path.exists(output_dir) is False:
            os.mkdir('../sad')
    n0 = get_sad_data(site_name, sad_path)
    S0 = len(n0)
    if S0 <= S_cutoff:
        print 'S0 less than or equal to cutoff value'
        return None
    N0 = sum(n0)
    rad = get_mete_rad(S0, N0)[0]
    out = np.zeros((len(n0), ), dtype = ('S10, i8, i8'))
    out['f0'] = site_name
    out['f1'] = n0
    out['f2'] = rad
    if to_file:
        filename = output_dir + site_name + '_obs_pred_rad.csv'
        writer = open(filename, 'wb') 
        datawriter = csv.writer(writer)
        for i in range(0, np.shape(out)[0]):
            datawriter.writerow(out[i, ])
        writer.close()
    else:
        return out

print 'Comparing METE and empirical SADs, ...'

site_names = ['bci','cocoli1','cocoli2','cross','sherman1','sherman2',
              'sherman3', 'serp', 'oosting', 'ferp', 'luquillo', 'graveyard',
              'landsend', 'rocky', 'bormann', 'woodbridge', 'baldmnt', 'bryan',
              'bigoak']

if len(sys.argv) > 1:
    site_index = int(sys.argv[1])
    site_names = [site_names[site_index]]

map(lambda x: get_obs_pred_sad(x, to_file=True), site_names)

print 'Comparing METE and empirical SADs, complete!'
