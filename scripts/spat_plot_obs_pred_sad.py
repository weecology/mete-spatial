from __future__ import division
import numpy as np
import csv
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.axes_grid.inset_locator import inset_axes
from mete import *
import macroecotools

def import_obs_pred_data(input_filename):
    data = np.genfromtxt(input_filename, dtype = "S15,f8,f8",
                                  names = ['site','obs','pred'],
                                  delimiter = ",")
    return data


def get_obs_pred_from_file(datasets, data_dir, filename):
    """Read obs and pred value from a file"""
    sites= []
    obs = []
    pred = []
    for i, dataset in enumerate(datasets):
        obs_pred_data = import_obs_pred_data(data_dir + dataset + filename)
        site_list = [dataset + '_' + x for x in list(obs_pred_data['site'])]
        sites.extend(site_list)
        obs.extend(list(obs_pred_data['obs']))
        pred.extend(list(obs_pred_data['pred']))
    return np.array(sites), np.array(obs), np.array(pred)

def plot_obs_pred(obs, pred, radius, loglog, ax = None, inset = False, sites = None):
    """Generic function to generate an observed vs predicted figure with 1:1 line"""
    if not ax:
        fig = plt.figure(figsize = (3.5, 3.5))
        ax = plt.subplot(111)

    axis_min = 0.9 * min(list(obs[obs > 0]) + list(pred[pred > 0]))
    if loglog:
        axis_max = 3 * max(list(obs)+list(pred))
    else:
        axis_max = 1.1 * max(list(obs)+list(pred))
    macroecotools.plot_color_by_pt_dens(np.array(pred), np.array(obs), radius, loglog=loglog, plot_obj = ax)      
    plt.plot([axis_min, axis_max],[axis_min, axis_max], 'k-')
    plt.xlim(axis_min, axis_max)
    plt.ylim(axis_min, axis_max)
    ax.tick_params(axis = 'both', which = 'major', labelsize = 6)
    if loglog:
        plt.annotate(r'$R^2$ = %0.2f' %macroecotools.obs_pred_rsquare(np.log10(obs[(obs != 0) * (pred != 0)]), np.log10(pred[(obs != 0) * (pred != 0)])),
                     xy = (0.05, 0.85), xycoords = 'axes fraction', fontsize = 7)
    else:
        plt.annotate(r'$R^2$ = %0.2f' %macroecotools.obs_pred_rsquare(obs, pred),
                     xy = (0.05, 0.85), xycoords = 'axes fraction', fontsize = 7)
    if inset:
        axins = inset_axes(ax, width="30%", height="30%", loc=4)
        if loglog:
            hist_mete_r2(sites[(obs != 0) * (pred != 0)], np.log10(obs[(obs != 0) * (pred != 0)]), 
                         np.log10(pred[(obs != 0) * (pred != 0)]))
        else:
            hist_mete_r2(sites, obs, pred)
        plt.setp(axins, xticks=[], yticks=[])
    return ax

def plot_obs_pred_sad(datasets, data_dir = "./data/", dest_dir = "./", radius = 2, inset = False):
    """Plot the observed vs predicted abundance for each species for multiple datasets."""
    rad_sites, rad_obs, rad_pred = get_obs_pred_from_file(datasets, data_dir, '_obs_pred_rad.csv')
    if inset:
        fig = plot_obs_pred(rad_obs, rad_pred, radius, 1, inset = True, sites = rad_sites)
    else: 
        fig = plot_obs_pred(rad_obs, rad_pred, radius, 1)
    fig.set_xlabel('Predicted abundance', labelpad = 4, size = 8)
    fig.set_ylabel('Observed abundance', labelpad = 4, size = 8)
    plt.savefig(dest_dir + 'obs_pred_sad.png', dpi = 400)

print 'Plotting obs-pred plot for SAD, ...'

datafile = open('./data/shrtnames.txt', 'r')
datareader = csv.reader(datafile)
shrtnames = []
for row in datareader:
    shrtnames.append(row)

shrtnames = shrtnames[0][0].split()

plot_obs_pred_sad(shrtnames, data_dir='./sad/',
                  dest_dir='./figs/sad_figs/')

for shrt_name in shrtnames:
    plot_obs_pred_sad([shrt_name], data_dir='./sad/', 
                      dest_dir='./figs/sad_figs/' + shrt_name + '_')
    
print 'Plotting obs-pred plot for SAD, complete!'
