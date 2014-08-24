# -*- coding: utf-8 -*-
"""
Created on Tue Aug 19 15:32:59 2014

@author: luis
"""

## This script will run through a few different SNR values between two
## circular gaussians. With each separation, we run many trials
## to assess bias, uncertainty, and the correlation matrix.

from __future__ import division
from pandas import Series, DataFrame
import pandas as pd
import numpy as np
import drawLibrary
import noiseLibrary
import galsim
import lmfit
import matplotlib.pyplot as plt
import scipy.linalg as scla
import matplotlib.gridspec as gridspec
import scipy.interpolate

# Parameters for object a
flux_a = 0            # We will loop through this value later
hlr_a = 1             # arcsec
e1_a = 0.0
e2_a = 0.0
x0_a = -1
y0_a = 0
n_a = 0.5

# Parameters for object b
flux_b = 0            # We will loop through this value later
hlr_b = hlr_a         # arcsec
e1_b = 0.0
e2_b = 0.0
x0_b = 1
y0_b = 0
n_b = 0.5

true_val = np.array([flux_a,hlr_a,e1_a,e2_a,x0_a,y0_a,flux_b,hlr_b,e1_b,e2_b,x0_b,y0_b])

# Track the separation
sep = x0_b - x0_a

# Galsim function definitions
sersic_func = galsim.Sersic

# Set the RNG
seed_1 = galsim.BaseDeviate(1)
seed_2 = galsim.BaseDeviate(2)
seed_3 = galsim.BaseDeviate(3)

# Image properties
pixel_scale = 1/5     # arcsec / pixel
x_len = y_len = 100            # pixel

# Use LSST defined sky noise for r-band
add_noise_flag = True
texp = 6900 # seconds
sbar = 26.8 # sky photons per second per pixel
sky_level = texp*sbar
sky_noise = np.sqrt(sky_level)

# psf properties
psf_flag = False
beta = 3
fwhm_psf = 0.6

# Calc function for SNR --> flux
SNR_to_flux, snr_points, flux_pts = noiseLibrary.calc_SNR_to_flux(hlr_a,e1_a,e2_a,x0_a,y0_a,n_a,
                                                        hlr_b,e1_b,e2_b,x0_b,y0_b,n_b,
                                                        psf_flag,beta,fwhm_psf,
                                                        x_len,y_len,pixel_scale,sersic_func,sersic_func,seed_1,seed_2,seed_3,
                                                        False,sky_level,sbar,texp,
                                                        1000,1000,10)
                                            
#plt.scatter(snr_points,flux_pts,c='b',alpha=0.5)
#plt.title('Flux vs SNR'); plt.xlabel('SNR'); plt.ylabel('Flux')
snr_points = np.array(snr_points); flux_pts = np.array(flux_pts) 
cond = np.logical_and(snr_points > 0, snr_points < 150)
flux_pts = flux_pts[cond]
snr_points = snr_points[cond]
#plt.xlim([0,np.max(snr_points)]); plt.ylim([0,np.max(flux_pts)])
SNR_to_flux = scipy.interpolate.interp1d(snr_points,flux_pts,kind='cubic')


# SNR range to loop through                                            
SNR_range = [100,40,30,20,15,10,5]
# Flux range to loop through
Flux_range = [1e6,5e5,1e5,1e4,1e3,1e2]
# number of trials
num_trials = 50

# Data to keep track of
resid_matrix = []
error_matrix = []
pull_matrix = []

for SNR in SNR_range:
    tot_flux = SNR_to_flux(SNR)
    flux_a = tot_flux/2
    flux_b = tot_flux/2
    true_val[0] = flux_a; true_val[6] = flux_b
    
    resid_vector = []
    error_vector = []
    pull_vector = []
    
    for i in xrange(0,num_trials):
        # Obtain instantiation
        im_no_noise, im_noise, best_fit, result = noiseLibrary.run_2_galaxy_full_params_simple(flux_a,hlr_a,e1_a,e2_a,x0_a,y0_a,n_a,
                                                                                               flux_b,hlr_b,e1_b,e2_b,x0_b,y0_b,n_b,
                                                                                               psf_flag,beta,fwhm_psf,
                                                                                               x_len,y_len,pixel_scale,sersic_func,sersic_func,seed_1,seed_2,seed_3,
                                                                                               add_noise_flag,sky_level)
                                                                                               
        # Obtain the covariance and correlation matrix.
        if result.covar is None:
            print "Covariance Matrix is Null"
            im_no_noise, im_noise, best_fit, result = noiseLibrary.run_2_galaxy_full_params_simple(flux_a,hlr_a,e1_a,e2_a,x0_a,y0_a,n_a,
                                                                                                   flux_b,hlr_b,e1_b,e2_b,x0_b,y0_b,n_b,
                                                                                                   psf_flag,beta,fwhm_psf,
                                                                                                   x_len,y_len,pixel_scale,sersic_func,sersic_func,seed_1,seed_2,seed_3,
                                                                                                   add_noise_flag,sky_level)
                                                                                       
        error_diag = (np.sqrt(np.diag(result.covar)))
        error_mat = np.outer(error_diag,error_diag)
        correlation_mat = result.covar/error_mat
        error_diag = error_diag.tolist()
        
        # Obtain the residuals
        resid = noiseLibrary.calc_resid(result.params,true_val)
        
        # Put residuals and error bars in a matrix and obtain the pull
        pull = np.array(resid)/error_diag
        # Put pull in this form in order to plot with imshow
        pull_mat = np.array([[pull[0]],[pull[1]],[pull[2]],[pull[3]],[pull[4]],[pull[5]],
                             [pull[6]],[pull[7]],[pull[8]],[pull[9]],[pull[10]],[pull[11]]])                                                                                       
        pull = pull.tolist()
        # Store the data in vectors
        resid_vector.append(resid)
        error_vector.append(error_diag)
        pull_vector.append(pull)

    resid_matrix.append(resid_vector)
    error_matrix.append(error_vector)
    pull_matrix.append(pull_vector)                                    

resid_matrix = np.array(resid_matrix)
error_matrix = np.array(error_matrix)
pull_matrix = np.array(pull_matrix)

# Title of plots
title = 'Residuals'

# If you want to graph pulls without renaming everything
graph_pulls = False
if graph_pulls:
    resid_copy = resid_matrix
    resid_matrix = pull_matrix
    title = 'Pulls'

# Obtain the matrices for each SNR from the 3d matrix
Resid_SNR_100 = resid_matrix[0,:,:]
Resid_SNR_40 = resid_matrix[1,:,:]
Resid_SNR_30 = resid_matrix[2,:,:]
Resid_SNR_20 = resid_matrix[3,:,:]
Resid_SNR_15 = resid_matrix[4,:,:]
Resid_SNR_10 = resid_matrix[5,:,:]
Resid_SNR_5 = resid_matrix[6,:,:]

gs = gridspec.GridSpec(2,2)
data_pts = False
# Plotting for e1 and e2 for objects a and b
fontsize = 13
fig = plt.figure(figsize=(20,11))
suptitle = 'Bias Anaylsis For Two Profiles With Sersic Index: $%.2f$\n Sep: $%.2f\/arcs$; $S=%.2f $; $t_{exp}=%.2fsec$; Trials = $%.2f$'%(n_a,sep,sbar,texp,num_trials)
plt.suptitle(suptitle,fontsize=fontsize)
ax1 = fig.add_subplot(gs[0,0])

col = 2
mark = 'x'
alpha = 0.05
mean_linewidth = 1
bar_linewidth = 2

max_pts = np.max([np.mean(Resid_SNR_100[:,col]),np.mean(Resid_SNR_40[:,col]),np.mean(Resid_SNR_30[:,col]),
                  np.mean(Resid_SNR_20[:,col]),np.mean(Resid_SNR_15[:,col]),np.mean(Resid_SNR_10[:,col]),np.mean(Resid_SNR_5[:,col])])
min_pts = np.min([np.mean(Resid_SNR_100[:,col]),np.mean(Resid_SNR_40[:,col]),np.mean(Resid_SNR_30[:,col]),
                  np.mean(Resid_SNR_20[:,col]),np.mean(Resid_SNR_15[:,col]),np.mean(Resid_SNR_10[:,col]),np.mean(Resid_SNR_5[:,col])])
max_std = np.max([np.std(Resid_SNR_100[:,col]),np.std(Resid_SNR_40[:,col]),np.std(Resid_SNR_30[:,col]),
                  np.std(Resid_SNR_20[:,col]),np.std(Resid_SNR_15[:,col]),np.std(Resid_SNR_10[:,col]),np.std(Resid_SNR_5[:,col])])

plt.ylim([min_pts-max_std,max_pts+max_std])

plt.title(title + ' of $e1$ for Object a vs SNR',fontsize=fontsize); plt.ylabel('Residuals of $e1_a$',fontsize=fontsize)
plt.xlabel('SNR',fontsize=fontsize)

if data_pts:
    plt.scatter(100*np.ones(num_trials),Resid_SNR_100[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([100],np.mean(Resid_SNR_100[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([100],np.mean(Resid_SNR_100[:,col]),yerr=np.std(Resid_SNR_100[:,col]),ecolor='b',elinewidth=bar_linewidth)
plt.errorbar([100],np.mean(Resid_SNR_100[:,col]),yerr=np.std(Resid_SNR_100[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')

if data_pts:
    plt.scatter(40*np.ones(num_trials),Resid_SNR_40[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([40],np.mean(Resid_SNR_40[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([40],np.mean(Resid_SNR_40[:,col]),yerr=np.std(Resid_SNR_40[:,col]),ecolor='b',elinewidth=bar_linewidth)
plt.errorbar([40],np.mean(Resid_SNR_40[:,col]),yerr=np.std(Resid_SNR_40[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')

if data_pts:    
    plt.scatter(30*np.ones(num_trials),Resid_SNR_30[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([30],np.mean(Resid_SNR_30[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([30],np.mean(Resid_SNR_30[:,col]),yerr=np.std(Resid_SNR_30[:,col]),ecolor='b',elinewidth=bar_linewidth)
plt.errorbar([30],np.mean(Resid_SNR_30[:,col]),yerr=np.std(Resid_SNR_30[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')

if data_pts:
    plt.scatter(20*np.ones(num_trials),Resid_SNR_20[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([20],np.mean(Resid_SNR_20[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([20],np.mean(Resid_SNR_20[:,col]),yerr=np.std(Resid_SNR_20[:,col]),ecolor='b',elinewidth=bar_linewidth)
plt.errorbar([20],np.mean(Resid_SNR_20[:,col]),yerr=np.std(Resid_SNR_20[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


if data_pts:    
    plt.scatter(15*np.ones(num_trials),Resid_SNR_15[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([15],np.mean(Resid_SNR_15[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([15],np.mean(Resid_SNR_15[:,col]),yerr=np.std(Resid_SNR_15[:,col]),ecolor='b',elinewidth=bar_linewidth)
plt.errorbar([15],np.mean(Resid_SNR_15[:,col]),yerr=np.std(Resid_SNR_15[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


if data_pts:
    plt.scatter(10*np.ones(num_trials),Resid_SNR_10[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([10],np.mean(Resid_SNR_10[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([10],np.mean(Resid_SNR_10[:,col]),yerr=np.std(Resid_SNR_10[:,col]),ecolor='b',elinewidth=bar_linewidth)
plt.errorbar([10],np.mean(Resid_SNR_10[:,col]),yerr=np.std(Resid_SNR_10[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')

if data_pts:
    plt.scatter(5*np.ones(num_trials),Resid_SNR_5[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([5],np.mean(Resid_SNR_5[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([5],np.mean(Resid_SNR_5[:,col]),yerr=np.std(Resid_SNR_5[:,col]),ecolor='b',elinewidth=bar_linewidth)
plt.errorbar([5],np.mean(Resid_SNR_5[:,col]),yerr=np.std(Resid_SNR_5[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


ax2 = fig.add_subplot(gs[0,1])
col = 3
mark = 'x'

max_pts = np.max([np.mean(Resid_SNR_100[:,col]),np.mean(Resid_SNR_40[:,col]),np.mean(Resid_SNR_30[:,col]),
                  np.mean(Resid_SNR_20[:,col]),np.mean(Resid_SNR_15[:,col]),np.mean(Resid_SNR_10[:,col]),np.mean(Resid_SNR_5[:,col])])
min_pts = np.min([np.mean(Resid_SNR_100[:,col]),np.mean(Resid_SNR_40[:,col]),np.mean(Resid_SNR_30[:,col]),
                  np.mean(Resid_SNR_20[:,col]),np.mean(Resid_SNR_15[:,col]),np.mean(Resid_SNR_10[:,col]),np.mean(Resid_SNR_5[:,col])])
max_std = np.max([np.std(Resid_SNR_100[:,col]),np.std(Resid_SNR_40[:,col]),np.std(Resid_SNR_30[:,col]),
                  np.std(Resid_SNR_20[:,col]),np.std(Resid_SNR_15[:,col]),np.std(Resid_SNR_10[:,col]),np.std(Resid_SNR_5[:,col])])

plt.ylim([min_pts-max_std,max_pts+max_std])


plt.title(title + ' of $e2$ for Object a vs SNR',fontsize=fontsize); plt.ylabel('Residuals of $e2_a$',fontsize=fontsize)
plt.xlabel('SNR',fontsize=fontsize)

if data_pts:
    plt.scatter(100*np.ones(num_trials),Resid_SNR_100[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([100],np.mean(Resid_SNR_100[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([100],np.mean(Resid_SNR_100[:,col]),yerr=np.std(Resid_SNR_100[:,col]),ecolor='b',elinewidth=bar_linewidth)
plt.errorbar([100],np.mean(Resid_SNR_100[:,col]),yerr=np.std(Resid_SNR_100[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


if data_pts:
    plt.scatter(40*np.ones(num_trials),Resid_SNR_40[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([40],np.mean(Resid_SNR_40[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([40],np.mean(Resid_SNR_40[:,col]),yerr=np.std(Resid_SNR_40[:,col]),ecolor='b',elinewidth=bar_linewidth)
plt.errorbar([40],np.mean(Resid_SNR_40[:,col]),yerr=np.std(Resid_SNR_40[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


if data_pts:    
    plt.scatter(30*np.ones(num_trials),Resid_SNR_30[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([30],np.mean(Resid_SNR_30[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([30],np.mean(Resid_SNR_30[:,col]),yerr=np.std(Resid_SNR_30[:,col]),ecolor='b',elinewidth=bar_linewidth)
plt.errorbar([30],np.mean(Resid_SNR_30[:,col]),yerr=np.std(Resid_SNR_30[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


if data_pts:
    plt.scatter(20*np.ones(num_trials),Resid_SNR_20[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([20],np.mean(Resid_SNR_20[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([20],np.mean(Resid_SNR_20[:,col]),yerr=np.std(Resid_SNR_20[:,col]),ecolor='b',elinewidth=bar_linewidth)
plt.errorbar([20],np.mean(Resid_SNR_20[:,col]),yerr=np.std(Resid_SNR_20[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


if data_pts:    
    plt.scatter(15*np.ones(num_trials),Resid_SNR_15[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([15],np.mean(Resid_SNR_15[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([15],np.mean(Resid_SNR_15[:,col]),yerr=np.std(Resid_SNR_15[:,col]),ecolor='b',elinewidth=bar_linewidth)
plt.errorbar([15],np.mean(Resid_SNR_15[:,col]),yerr=np.std(Resid_SNR_15[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


if data_pts:
    plt.scatter(10*np.ones(num_trials),Resid_SNR_10[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([10],np.mean(Resid_SNR_10[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([10],np.mean(Resid_SNR_10[:,col]),yerr=np.std(Resid_SNR_10[:,col]),ecolor='b',elinewidth=bar_linewidth)
plt.errorbar([10],np.mean(Resid_SNR_10[:,col]),yerr=np.std(Resid_SNR_10[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')

if data_pts:
    plt.scatter(5*np.ones(num_trials),Resid_SNR_10[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([5],np.mean(Resid_SNR_5[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([5],np.mean(Resid_SNR_5[:,col]),yerr=np.std(Resid_SNR_5[:,col]),ecolor='b',elinewidth=bar_linewidth)
plt.errorbar([5],np.mean(Resid_SNR_5[:,col]),yerr=np.std(Resid_SNR_5[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


ax3 = fig.add_subplot(gs[1,0])
col = 2 + 6
mark = 'o'

max_pts = np.max([np.mean(Resid_SNR_100[:,col]),np.mean(Resid_SNR_40[:,col]),np.mean(Resid_SNR_30[:,col]),
                  np.mean(Resid_SNR_20[:,col]),np.mean(Resid_SNR_15[:,col]),np.mean(Resid_SNR_10[:,col]),np.mean(Resid_SNR_5[:,col])])
min_pts = np.min([np.mean(Resid_SNR_100[:,col]),np.mean(Resid_SNR_40[:,col]),np.mean(Resid_SNR_30[:,col]),
                  np.mean(Resid_SNR_20[:,col]),np.mean(Resid_SNR_15[:,col]),np.mean(Resid_SNR_10[:,col]),np.mean(Resid_SNR_5[:,col])])
max_std = np.max([np.std(Resid_SNR_100[:,col]),np.std(Resid_SNR_40[:,col]),np.std(Resid_SNR_30[:,col]),
                  np.std(Resid_SNR_20[:,col]),np.std(Resid_SNR_15[:,col]),np.std(Resid_SNR_10[:,col]),np.std(Resid_SNR_5[:,col])])

plt.ylim([min_pts-max_std,max_pts+max_std])


plt.title(title + ' of $e1$ for Object b vs SNR',fontsize=fontsize); plt.ylabel('Residuals of $e1_b$',fontsize=fontsize)
plt.xlabel('SNR',fontsize=fontsize)

if data_pts:
    plt.scatter(100*np.ones(num_trials),Resid_SNR_100[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([100],np.mean(Resid_SNR_100[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([100],np.mean(Resid_SNR_100[:,col]),yerr=np.std(Resid_SNR_100[:,col]),ecolor='b',elinewidth=bar_linewidth)
plt.errorbar([100],np.mean(Resid_SNR_100[:,col]),yerr=np.std(Resid_SNR_100[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


if data_pts:
    plt.scatter(40*np.ones(num_trials),Resid_SNR_40[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([40],np.mean(Resid_SNR_40[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([40],np.mean(Resid_SNR_40[:,col]),yerr=np.std(Resid_SNR_40[:,col]),ecolor='b',elinewidth=bar_linewidth)
plt.errorbar([40],np.mean(Resid_SNR_40[:,col]),yerr=np.std(Resid_SNR_40[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


if data_pts:    
    plt.scatter(30*np.ones(num_trials),Resid_SNR_30[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([30],np.mean(Resid_SNR_30[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([30],np.mean(Resid_SNR_30[:,col]),yerr=np.std(Resid_SNR_30[:,col]),ecolor='b',elinewidth=bar_linewidth)
plt.errorbar([30],np.mean(Resid_SNR_30[:,col]),yerr=np.std(Resid_SNR_30[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


if data_pts:
    plt.scatter(20*np.ones(num_trials),Resid_SNR_20[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([20],np.mean(Resid_SNR_20[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([20],np.mean(Resid_SNR_20[:,col]),yerr=np.std(Resid_SNR_20[:,col]),ecolor='b',elinewidth=bar_linewidth)
plt.errorbar([20],np.mean(Resid_SNR_20[:,col]),yerr=np.std(Resid_SNR_20[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


if data_pts:    
    plt.scatter(15*np.ones(num_trials),Resid_SNR_15[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([15],np.mean(Resid_SNR_15[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([15],np.mean(Resid_SNR_15[:,col]),yerr=np.std(Resid_SNR_15[:,col]),ecolor='b',elinewidth=bar_linewidth)
plt.errorbar([15],np.mean(Resid_SNR_15[:,col]),yerr=np.std(Resid_SNR_15[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


if data_pts:
    plt.scatter(10*np.ones(num_trials),Resid_SNR_10[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([10],np.mean(Resid_SNR_10[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([10],np.mean(Resid_SNR_10[:,col]),yerr=np.std(Resid_SNR_10[:,col]),ecolor='b',elinewidth=bar_linewidth)
plt.errorbar([10],np.mean(Resid_SNR_10[:,col]),yerr=np.std(Resid_SNR_10[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')

if data_pts:
    plt.scatter(5*np.ones(num_trials),Resid_SNR_10[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([5],np.mean(Resid_SNR_5[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([5],np.mean(Resid_SNR_5[:,col]),yerr=np.std(Resid_SNR_5[:,col]),ecolor='b',elinewidth=bar_linewidth)
plt.errorbar([5],np.mean(Resid_SNR_5[:,col]),yerr=np.std(Resid_SNR_5[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')

ax4 = fig.add_subplot(gs[1,1])
col = 3 + 6
mark = 'x'

max_pts = np.max([np.mean(Resid_SNR_100[:,col]),np.mean(Resid_SNR_40[:,col]),np.mean(Resid_SNR_30[:,col]),
                  np.mean(Resid_SNR_20[:,col]),np.mean(Resid_SNR_15[:,col]),np.mean(Resid_SNR_10[:,col]),np.mean(Resid_SNR_5[:,col])])
min_pts = np.min([np.mean(Resid_SNR_100[:,col]),np.mean(Resid_SNR_40[:,col]),np.mean(Resid_SNR_30[:,col]),
                  np.mean(Resid_SNR_20[:,col]),np.mean(Resid_SNR_15[:,col]),np.mean(Resid_SNR_10[:,col]),np.mean(Resid_SNR_5[:,col])])
max_std = np.max([np.std(Resid_SNR_100[:,col]),np.std(Resid_SNR_40[:,col]),np.std(Resid_SNR_30[:,col]),
                  np.std(Resid_SNR_20[:,col]),np.std(Resid_SNR_15[:,col]),np.std(Resid_SNR_10[:,col]),np.std(Resid_SNR_5[:,col])])

plt.ylim([min_pts-max_std,max_pts+max_std])

plt.title(title + ' of $e2$ for Object b vs SNR',fontsize=fontsize); plt.ylabel('Residuals of $e2_b$',fontsize=fontsize)
plt.xlabel('SNR',fontsize=fontsize)

if data_pts:
    plt.scatter(100*np.ones(num_trials),Resid_SNR_100[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([100],np.mean(Resid_SNR_100[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([100],np.mean(Resid_SNR_100[:,col]),yerr=np.std(Resid_SNR_100[:,col]),ecolor='b',elinewidth=bar_linewidth)
plt.errorbar([100],np.mean(Resid_SNR_100[:,col]),yerr=np.std(Resid_SNR_100[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


if data_pts:
    plt.scatter(40*np.ones(num_trials),Resid_SNR_40[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([40],np.mean(Resid_SNR_40[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([40],np.mean(Resid_SNR_40[:,col]),yerr=np.std(Resid_SNR_40[:,col]),ecolor='b',elinewidth=bar_linewidth)
plt.errorbar([40],np.mean(Resid_SNR_40[:,col]),yerr=np.std(Resid_SNR_40[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


if data_pts:    
    plt.scatter(30*np.ones(num_trials),Resid_SNR_30[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([30],np.mean(Resid_SNR_30[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([30],np.mean(Resid_SNR_30[:,col]),yerr=np.std(Resid_SNR_30[:,col]),ecolor='b',elinewidth=bar_linewidth)
plt.errorbar([30],np.mean(Resid_SNR_30[:,col]),yerr=np.std(Resid_SNR_30[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


if data_pts:
    plt.scatter(20*np.ones(num_trials),Resid_SNR_20[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([20],np.mean(Resid_SNR_20[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([20],np.mean(Resid_SNR_20[:,col]),yerr=np.std(Resid_SNR_20[:,col]),ecolor='b',elinewidth=bar_linewidth)
plt.errorbar([20],np.mean(Resid_SNR_20[:,col]),yerr=np.std(Resid_SNR_20[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


if data_pts:    
    plt.scatter(15*np.ones(num_trials),Resid_SNR_15[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([15],np.mean(Resid_SNR_15[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([15],np.mean(Resid_SNR_15[:,col]),yerr=np.std(Resid_SNR_15[:,col]),ecolor='b',elinewidth=bar_linewidth)
plt.errorbar([15],np.mean(Resid_SNR_15[:,col]),yerr=np.std(Resid_SNR_15[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


if data_pts:
    plt.scatter(10*np.ones(num_trials),Resid_SNR_10[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([10],np.mean(Resid_SNR_10[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([10],np.mean(Resid_SNR_10[:,col]),yerr=np.std(Resid_SNR_10[:,col]),ecolor='b',elinewidth=bar_linewidth)
plt.errorbar([10],np.mean(Resid_SNR_10[:,col]),yerr=np.std(Resid_SNR_10[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')

if data_pts:
    plt.scatter(5*np.ones(num_trials),Resid_SNR_10[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([5],np.mean(Resid_SNR_5[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([5],np.mean(Resid_SNR_5[:,col]),yerr=np.std(Resid_SNR_5[:,col]),ecolor='b',elinewidth=bar_linewidth)
plt.errorbar([5],np.mean(Resid_SNR_5[:,col]),yerr=np.std(Resid_SNR_5[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')



# Plotting for HLR and Flux for a and b----------------------------------------
fig = plt.figure(figsize=(20,11))
plt.suptitle(suptitle,fontsize=fontsize)
ax1 = fig.add_subplot(gs[0,0])
col = 0
mark = 'o'

max_pts = np.max([np.mean(Resid_SNR_100[:,col]),np.mean(Resid_SNR_40[:,col]),np.mean(Resid_SNR_30[:,col]),
                  np.mean(Resid_SNR_20[:,col]),np.mean(Resid_SNR_15[:,col]),np.mean(Resid_SNR_10[:,col]),np.mean(Resid_SNR_5[:,col])])
min_pts = np.min([np.mean(Resid_SNR_100[:,col]),np.mean(Resid_SNR_40[:,col]),np.mean(Resid_SNR_30[:,col]),
                  np.mean(Resid_SNR_20[:,col]),np.mean(Resid_SNR_15[:,col]),np.mean(Resid_SNR_10[:,col]),np.mean(Resid_SNR_5[:,col])])
max_std = np.max([np.std(Resid_SNR_100[:,col]),np.std(Resid_SNR_40[:,col]),np.std(Resid_SNR_30[:,col]),
                  np.std(Resid_SNR_20[:,col]),np.std(Resid_SNR_15[:,col]),np.std(Resid_SNR_10[:,col]),np.std(Resid_SNR_5[:,col])])

plt.ylim([min_pts-max_std,max_pts+max_std])


plt.title(title + ' of $Flux$ for Object a vs SNR',fontsize=fontsize); plt.ylabel('Residuals of $Flux_a$',fontsize=fontsize)
plt.xlabel('SNR',fontsize=fontsize)


if data_pts:
    plt.scatter(100*np.ones(num_trials),Resid_SNR_100[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([100],np.mean(Resid_SNR_100[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([100],np.mean(Resid_SNR_100[:,col]),yerr=np.std(Resid_SNR_100[:,col]),ecolor='b',elinewidth=bar_linewidth)
plt.errorbar([100],np.mean(Resid_SNR_100[:,col]),yerr=np.std(Resid_SNR_100[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


if data_pts:
    plt.scatter(40*np.ones(num_trials),Resid_SNR_40[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([40],np.mean(Resid_SNR_40[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([40],np.mean(Resid_SNR_40[:,col]),yerr=np.std(Resid_SNR_40[:,col]),ecolor='b',elinewidth=bar_linewidth)
plt.errorbar([40],np.mean(Resid_SNR_40[:,col]),yerr=np.std(Resid_SNR_40[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


if data_pts:    
    plt.scatter(30*np.ones(num_trials),Resid_SNR_30[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([30],np.mean(Resid_SNR_30[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([30],np.mean(Resid_SNR_30[:,col]),yerr=np.std(Resid_SNR_30[:,col]),ecolor='b',elinewidth=bar_linewidth)
plt.errorbar([30],np.mean(Resid_SNR_30[:,col]),yerr=np.std(Resid_SNR_30[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


if data_pts:
    plt.scatter(20*np.ones(num_trials),Resid_SNR_20[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([20],np.mean(Resid_SNR_20[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([20],np.mean(Resid_SNR_20[:,col]),yerr=np.std(Resid_SNR_20[:,col]),ecolor='b',elinewidth=bar_linewidth)
plt.errorbar([20],np.mean(Resid_SNR_20[:,col]),yerr=np.std(Resid_SNR_20[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


if data_pts:    
    plt.scatter(15*np.ones(num_trials),Resid_SNR_15[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([15],np.mean(Resid_SNR_15[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([15],np.mean(Resid_SNR_15[:,col]),yerr=np.std(Resid_SNR_15[:,col]),ecolor='b',elinewidth=bar_linewidth)
plt.errorbar([15],np.mean(Resid_SNR_15[:,col]),yerr=np.std(Resid_SNR_15[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


if data_pts:
    plt.scatter(10*np.ones(num_trials),Resid_SNR_10[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([10],np.mean(Resid_SNR_10[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([10],np.mean(Resid_SNR_10[:,col]),yerr=np.std(Resid_SNR_10[:,col]),ecolor='b',elinewidth=bar_linewidth)
plt.errorbar([10],np.mean(Resid_SNR_10[:,col]),yerr=np.std(Resid_SNR_10[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')

if data_pts:
    plt.scatter(5*np.ones(num_trials),Resid_SNR_10[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([5],np.mean(Resid_SNR_5[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([5],np.mean(Resid_SNR_5[:,col]),yerr=np.std(Resid_SNR_5[:,col]),ecolor='b',elinewidth=bar_linewidth)
plt.errorbar([5],np.mean(Resid_SNR_5[:,col]),yerr=np.std(Resid_SNR_5[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


ax2 = fig.add_subplot(gs[0,1])
col = 1
mark = 'x'

max_pts = np.max([np.mean(Resid_SNR_100[:,col]),np.mean(Resid_SNR_40[:,col]),np.mean(Resid_SNR_30[:,col]),
                  np.mean(Resid_SNR_20[:,col]),np.mean(Resid_SNR_15[:,col]),np.mean(Resid_SNR_10[:,col]),np.mean(Resid_SNR_5[:,col])])
min_pts = np.min([np.mean(Resid_SNR_100[:,col]),np.mean(Resid_SNR_40[:,col]),np.mean(Resid_SNR_30[:,col]),
                  np.mean(Resid_SNR_20[:,col]),np.mean(Resid_SNR_15[:,col]),np.mean(Resid_SNR_10[:,col]),np.mean(Resid_SNR_5[:,col])])
max_std = np.max([np.std(Resid_SNR_100[:,col]),np.std(Resid_SNR_40[:,col]),np.std(Resid_SNR_30[:,col]),
                  np.std(Resid_SNR_20[:,col]),np.std(Resid_SNR_15[:,col]),np.std(Resid_SNR_10[:,col]),np.std(Resid_SNR_5[:,col])])

plt.ylim([min_pts-max_std,max_pts+max_std])


plt.title(title + ' of $Hlr$ for Object a vs SNR',fontsize=fontsize); plt.ylabel('Residuals of $Hlr_a$',fontsize=fontsize)
plt.xlabel('SNR',fontsize=fontsize)

if data_pts:
    plt.scatter(100*np.ones(num_trials),Resid_SNR_100[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([100],np.mean(Resid_SNR_100[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([100],np.mean(Resid_SNR_100[:,col]),yerr=np.std(Resid_SNR_100[:,col]),ecolor='b',elinewidth=bar_linewidth)
plt.errorbar([100],np.mean(Resid_SNR_100[:,col]),yerr=np.std(Resid_SNR_100[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


if data_pts:
    plt.scatter(40*np.ones(num_trials),Resid_SNR_40[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([40],np.mean(Resid_SNR_40[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([40],np.mean(Resid_SNR_40[:,col]),yerr=np.std(Resid_SNR_40[:,col]),ecolor='b',elinewidth=bar_linewidth)
plt.errorbar([40],np.mean(Resid_SNR_40[:,col]),yerr=np.std(Resid_SNR_40[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


if data_pts:    
    plt.scatter(30*np.ones(num_trials),Resid_SNR_30[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([30],np.mean(Resid_SNR_30[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([30],np.mean(Resid_SNR_30[:,col]),yerr=np.std(Resid_SNR_30[:,col]),ecolor='b',elinewidth=bar_linewidth)
plt.errorbar([30],np.mean(Resid_SNR_30[:,col]),yerr=np.std(Resid_SNR_30[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


if data_pts:
    plt.scatter(20*np.ones(num_trials),Resid_SNR_20[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([20],np.mean(Resid_SNR_20[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([20],np.mean(Resid_SNR_20[:,col]),yerr=np.std(Resid_SNR_20[:,col]),ecolor='b',elinewidth=bar_linewidth)
plt.errorbar([20],np.mean(Resid_SNR_20[:,col]),yerr=np.std(Resid_SNR_20[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


if data_pts:    
    plt.scatter(15*np.ones(num_trials),Resid_SNR_15[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([15],np.mean(Resid_SNR_15[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([15],np.mean(Resid_SNR_15[:,col]),yerr=np.std(Resid_SNR_15[:,col]),ecolor='b',elinewidth=bar_linewidth)
plt.errorbar([15],np.mean(Resid_SNR_15[:,col]),yerr=np.std(Resid_SNR_15[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


if data_pts:
    plt.scatter(10*np.ones(num_trials),Resid_SNR_10[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([10],np.mean(Resid_SNR_10[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([10],np.mean(Resid_SNR_10[:,col]),yerr=np.std(Resid_SNR_10[:,col]),ecolor='b',elinewidth=bar_linewidth)
plt.errorbar([10],np.mean(Resid_SNR_10[:,col]),yerr=np.std(Resid_SNR_10[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')

if data_pts:
    plt.scatter(5*np.ones(num_trials),Resid_SNR_10[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([5],np.mean(Resid_SNR_5[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([5],np.mean(Resid_SNR_5[:,col]),yerr=np.std(Resid_SNR_5[:,col]),ecolor='b',elinewidth=bar_linewidth)
plt.errorbar([5],np.mean(Resid_SNR_5[:,col]),yerr=np.std(Resid_SNR_5[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


ax3 = fig.add_subplot(gs[1,0])
col = 0 + 6
mark = 'o'

max_pts = np.max([np.mean(Resid_SNR_100[:,col]),np.mean(Resid_SNR_40[:,col]),np.mean(Resid_SNR_30[:,col]),
                  np.mean(Resid_SNR_20[:,col]),np.mean(Resid_SNR_15[:,col]),np.mean(Resid_SNR_10[:,col]),np.mean(Resid_SNR_5[:,col])])
min_pts = np.min([np.mean(Resid_SNR_100[:,col]),np.mean(Resid_SNR_40[:,col]),np.mean(Resid_SNR_30[:,col]),
                  np.mean(Resid_SNR_20[:,col]),np.mean(Resid_SNR_15[:,col]),np.mean(Resid_SNR_10[:,col]),np.mean(Resid_SNR_5[:,col])])
max_std = np.max([np.std(Resid_SNR_100[:,col]),np.std(Resid_SNR_40[:,col]),np.std(Resid_SNR_30[:,col]),
                  np.std(Resid_SNR_20[:,col]),np.std(Resid_SNR_15[:,col]),np.std(Resid_SNR_10[:,col]),np.std(Resid_SNR_5[:,col])])

plt.ylim([min_pts-max_std,max_pts+max_std])


plt.title(title + ' of $Flux$ for Object b vs SNR',fontsize=fontsize); plt.ylabel('Residuals of $Flux_b$',fontsize=fontsize)
plt.xlabel('SNR',fontsize=fontsize)

if data_pts:
    plt.scatter(100*np.ones(num_trials),Resid_SNR_100[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([100],np.mean(Resid_SNR_100[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([100],np.mean(Resid_SNR_100[:,col]),yerr=np.std(Resid_SNR_100[:,col]),ecolor='b',elinewidth=bar_linewidth)
plt.errorbar([100],np.mean(Resid_SNR_100[:,col]),yerr=np.std(Resid_SNR_100[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


if data_pts:
    plt.scatter(40*np.ones(num_trials),Resid_SNR_40[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([40],np.mean(Resid_SNR_40[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([40],np.mean(Resid_SNR_40[:,col]),yerr=np.std(Resid_SNR_40[:,col]),ecolor='b',elinewidth=bar_linewidth)
plt.errorbar([40],np.mean(Resid_SNR_40[:,col]),yerr=np.std(Resid_SNR_40[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


if data_pts:    
    plt.scatter(30*np.ones(num_trials),Resid_SNR_30[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([30],np.mean(Resid_SNR_30[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([30],np.mean(Resid_SNR_30[:,col]),yerr=np.std(Resid_SNR_30[:,col]),ecolor='b',elinewidth=bar_linewidth)
plt.errorbar([30],np.mean(Resid_SNR_30[:,col]),yerr=np.std(Resid_SNR_30[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


if data_pts:
    plt.scatter(20*np.ones(num_trials),Resid_SNR_20[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([20],np.mean(Resid_SNR_20[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([20],np.mean(Resid_SNR_20[:,col]),yerr=np.std(Resid_SNR_20[:,col]),ecolor='b',elinewidth=bar_linewidth)
plt.errorbar([20],np.mean(Resid_SNR_20[:,col]),yerr=np.std(Resid_SNR_20[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


if data_pts:    
    plt.scatter(15*np.ones(num_trials),Resid_SNR_15[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([15],np.mean(Resid_SNR_15[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([15],np.mean(Resid_SNR_15[:,col]),yerr=np.std(Resid_SNR_15[:,col]),ecolor='b',elinewidth=bar_linewidth)
plt.errorbar([15],np.mean(Resid_SNR_15[:,col]),yerr=np.std(Resid_SNR_15[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


if data_pts:
    plt.scatter(10*np.ones(num_trials),Resid_SNR_10[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([10],np.mean(Resid_SNR_10[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([10],np.mean(Resid_SNR_10[:,col]),yerr=np.std(Resid_SNR_10[:,col]),ecolor='b',elinewidth=bar_linewidth)
plt.errorbar([10],np.mean(Resid_SNR_10[:,col]),yerr=np.std(Resid_SNR_10[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')

if data_pts:
    plt.scatter(5*np.ones(num_trials),Resid_SNR_10[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([5],np.mean(Resid_SNR_5[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([5],np.mean(Resid_SNR_5[:,col]),yerr=np.std(Resid_SNR_5[:,col]),ecolor='b',elinewidth=bar_linewidth)
plt.errorbar([5],np.mean(Resid_SNR_5[:,col]),yerr=np.std(Resid_SNR_5[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


ax4 = fig.add_subplot(gs[1,1])
col = 1 + 6
mark = 'x'

plt.title(title + ' of $Hlr$ for Object b vs SNR',fontsize=fontsize); plt.ylabel('Residuals of $Hlr_b$',fontsize=fontsize)
plt.xlabel('SNR',fontsize=fontsize)

max_pts = np.max([np.mean(Resid_SNR_100[:,col]),np.mean(Resid_SNR_40[:,col]),np.mean(Resid_SNR_30[:,col]),
                  np.mean(Resid_SNR_20[:,col]),np.mean(Resid_SNR_15[:,col]),np.mean(Resid_SNR_10[:,col]),np.mean(Resid_SNR_5[:,col])])
min_pts = np.min([np.mean(Resid_SNR_100[:,col]),np.mean(Resid_SNR_40[:,col]),np.mean(Resid_SNR_30[:,col]),
                  np.mean(Resid_SNR_20[:,col]),np.mean(Resid_SNR_15[:,col]),np.mean(Resid_SNR_10[:,col]),np.mean(Resid_SNR_5[:,col])])
max_std = np.max([np.std(Resid_SNR_100[:,col]),np.std(Resid_SNR_40[:,col]),np.std(Resid_SNR_30[:,col]),
                  np.std(Resid_SNR_20[:,col]),np.std(Resid_SNR_15[:,col]),np.std(Resid_SNR_10[:,col]),np.std(Resid_SNR_5[:,col])])

plt.ylim([min_pts-max_std,max_pts+max_std])


if data_pts:
    plt.scatter(100*np.ones(num_trials),Resid_SNR_100[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([100],np.mean(Resid_SNR_100[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([100],np.mean(Resid_SNR_100[:,col]),yerr=np.std(Resid_SNR_100[:,col]),ecolor='b',elinewidth=bar_linewidth)
plt.errorbar([100],np.mean(Resid_SNR_100[:,col]),yerr=np.std(Resid_SNR_100[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


if data_pts:
    plt.scatter(40*np.ones(num_trials),Resid_SNR_40[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([40],np.mean(Resid_SNR_40[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([40],np.mean(Resid_SNR_40[:,col]),yerr=np.std(Resid_SNR_40[:,col]),ecolor='b',elinewidth=bar_linewidth)
plt.errorbar([40],np.mean(Resid_SNR_40[:,col]),yerr=np.std(Resid_SNR_40[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


if data_pts:    
    plt.scatter(30*np.ones(num_trials),Resid_SNR_30[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([30],np.mean(Resid_SNR_30[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([30],np.mean(Resid_SNR_30[:,col]),yerr=np.std(Resid_SNR_30[:,col]),ecolor='b',elinewidth=bar_linewidth)
plt.errorbar([30],np.mean(Resid_SNR_30[:,col]),yerr=np.std(Resid_SNR_30[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


if data_pts:
    plt.scatter(20*np.ones(num_trials),Resid_SNR_20[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([20],np.mean(Resid_SNR_20[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([20],np.mean(Resid_SNR_20[:,col]),yerr=np.std(Resid_SNR_20[:,col]),ecolor='b',elinewidth=bar_linewidth)
plt.errorbar([20],np.mean(Resid_SNR_20[:,col]),yerr=np.std(Resid_SNR_20[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


if data_pts:    
    plt.scatter(15*np.ones(num_trials),Resid_SNR_15[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([15],np.mean(Resid_SNR_15[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([15],np.mean(Resid_SNR_15[:,col]),yerr=np.std(Resid_SNR_15[:,col]),ecolor='b',elinewidth=bar_linewidth)
plt.errorbar([15],np.mean(Resid_SNR_15[:,col]),yerr=np.std(Resid_SNR_15[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


if data_pts:
    plt.scatter(10*np.ones(num_trials),Resid_SNR_10[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([10],np.mean(Resid_SNR_10[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([10],np.mean(Resid_SNR_10[:,col]),yerr=np.std(Resid_SNR_10[:,col]),ecolor='b',elinewidth=bar_linewidth)
plt.errorbar([10],np.mean(Resid_SNR_10[:,col]),yerr=np.std(Resid_SNR_10[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')

if data_pts:
    plt.scatter(5*np.ones(num_trials),Resid_SNR_10[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([5],np.mean(Resid_SNR_5[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([5],np.mean(Resid_SNR_5[:,col]),yerr=np.std(Resid_SNR_5[:,col]),ecolor='b',elinewidth=bar_linewidth)
plt.errorbar([5],np.mean(Resid_SNR_5[:,col]),yerr=np.std(Resid_SNR_5[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')

plt.show()

######### Plotting for just the error on the mean -----------------------------

fig = plt.figure(figsize=(20,11))
plt.suptitle(suptitle,fontsize=fontsize)
ax1 = fig.add_subplot(gs[0,0])

col = 2
mark = 'x'
alpha = 0.05
mean_linewidth = 1
bar_linewidth = 2

plt.title(title + ' of $e1$ for Object a vs SNR',fontsize=fontsize); plt.ylabel('Residuals of $e1_a$',fontsize=fontsize)
plt.xlabel('SNR',fontsize=fontsize)

max_pts = np.max([np.mean(Resid_SNR_100[:,col]),np.mean(Resid_SNR_40[:,col]),np.mean(Resid_SNR_30[:,col]),
                  np.mean(Resid_SNR_20[:,col]),np.mean(Resid_SNR_15[:,col]),np.mean(Resid_SNR_10[:,col]),np.mean(Resid_SNR_5[:,col])])
min_pts = np.min([np.mean(Resid_SNR_100[:,col]),np.mean(Resid_SNR_40[:,col]),np.mean(Resid_SNR_30[:,col]),
                  np.mean(Resid_SNR_20[:,col]),np.mean(Resid_SNR_15[:,col]),np.mean(Resid_SNR_10[:,col]),np.mean(Resid_SNR_5[:,col])])
max_std = (np.max([np.std(Resid_SNR_100[:,col]),np.std(Resid_SNR_40[:,col]),np.std(Resid_SNR_30[:,col]),
                  np.std(Resid_SNR_20[:,col]),np.std(Resid_SNR_15[:,col]),np.std(Resid_SNR_10[:,col]),np.std(Resid_SNR_5[:,col])]))/np.sqrt(num_trials)

plt.ylim([min_pts-max_std,max_pts+max_std])

if data_pts:
    plt.scatter(100*np.ones(num_trials),Resid_SNR_100[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([100],np.mean(Resid_SNR_100[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([100],np.mean(Resid_SNR_100[:,col]),yerr=np.std(Resid_SNR_100[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')

if data_pts:
    plt.scatter(40*np.ones(num_trials),Resid_SNR_40[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([40],np.mean(Resid_SNR_40[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([40],np.mean(Resid_SNR_40[:,col]),yerr=np.std(Resid_SNR_40[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')

if data_pts:    
    plt.scatter(30*np.ones(num_trials),Resid_SNR_30[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([30],np.mean(Resid_SNR_30[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([30],np.mean(Resid_SNR_30[:,col]),yerr=np.std(Resid_SNR_30[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')

if data_pts:
    plt.scatter(20*np.ones(num_trials),Resid_SNR_20[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([20],np.mean(Resid_SNR_20[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([20],np.mean(Resid_SNR_20[:,col]),yerr=np.std(Resid_SNR_20[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


if data_pts:    
    plt.scatter(15*np.ones(num_trials),Resid_SNR_15[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([15],np.mean(Resid_SNR_15[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([15],np.mean(Resid_SNR_15[:,col]),yerr=np.std(Resid_SNR_15[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


if data_pts:
    plt.scatter(10*np.ones(num_trials),Resid_SNR_10[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([10],np.mean(Resid_SNR_10[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([10],np.mean(Resid_SNR_10[:,col]),yerr=np.std(Resid_SNR_10[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')

if data_pts:
    plt.scatter(5*np.ones(num_trials),Resid_SNR_10[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([5],np.mean(Resid_SNR_5[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([5],np.mean(Resid_SNR_5[:,col]),yerr=np.std(Resid_SNR_5[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


ax2 = fig.add_subplot(gs[0,1])
col = 3
mark = 'x'

plt.title(title + ' of $e2$ for Object a vs SNR',fontsize=fontsize); plt.ylabel('Residuals of $e2_a$',fontsize=fontsize)
plt.xlabel('SNR',fontsize=fontsize)
max_pts = np.max([np.mean(Resid_SNR_100[:,col]),np.mean(Resid_SNR_40[:,col]),np.mean(Resid_SNR_30[:,col]),
                  np.mean(Resid_SNR_20[:,col]),np.mean(Resid_SNR_15[:,col]),np.mean(Resid_SNR_10[:,col]),np.mean(Resid_SNR_5[:,col])])
min_pts = np.min([np.mean(Resid_SNR_100[:,col]),np.mean(Resid_SNR_40[:,col]),np.mean(Resid_SNR_30[:,col]),
                  np.mean(Resid_SNR_20[:,col]),np.mean(Resid_SNR_15[:,col]),np.mean(Resid_SNR_10[:,col]),np.mean(Resid_SNR_5[:,col])])
max_std = (np.max([np.std(Resid_SNR_100[:,col]),np.std(Resid_SNR_40[:,col]),np.std(Resid_SNR_30[:,col]),
                  np.std(Resid_SNR_20[:,col]),np.std(Resid_SNR_15[:,col]),np.std(Resid_SNR_10[:,col]),np.std(Resid_SNR_5[:,col])]))/np.sqrt(num_trials)

plt.ylim([min_pts-max_std,max_pts+max_std])

if data_pts:
    plt.scatter(100*np.ones(num_trials),Resid_SNR_100[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([100],np.mean(Resid_SNR_100[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([100],np.mean(Resid_SNR_100[:,col]),yerr=np.std(Resid_SNR_100[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')

if data_pts:
    plt.scatter(40*np.ones(num_trials),Resid_SNR_40[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([40],np.mean(Resid_SNR_40[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([40],np.mean(Resid_SNR_40[:,col]),yerr=np.std(Resid_SNR_40[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')

if data_pts:    
    plt.scatter(30*np.ones(num_trials),Resid_SNR_30[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([30],np.mean(Resid_SNR_30[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([30],np.mean(Resid_SNR_30[:,col]),yerr=np.std(Resid_SNR_30[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')

if data_pts:
    plt.scatter(20*np.ones(num_trials),Resid_SNR_20[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([20],np.mean(Resid_SNR_20[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([20],np.mean(Resid_SNR_20[:,col]),yerr=np.std(Resid_SNR_20[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


if data_pts:    
    plt.scatter(15*np.ones(num_trials),Resid_SNR_15[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([15],np.mean(Resid_SNR_15[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([15],np.mean(Resid_SNR_15[:,col]),yerr=np.std(Resid_SNR_15[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


if data_pts:
    plt.scatter(10*np.ones(num_trials),Resid_SNR_10[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([10],np.mean(Resid_SNR_10[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([10],np.mean(Resid_SNR_10[:,col]),yerr=np.std(Resid_SNR_10[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')

if data_pts:
    plt.scatter(5*np.ones(num_trials),Resid_SNR_10[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([5],np.mean(Resid_SNR_5[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([5],np.mean(Resid_SNR_5[:,col]),yerr=np.std(Resid_SNR_5[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')



ax3 = fig.add_subplot(gs[1,0])
col = 2 + 6
mark = 'o'

plt.title(title + ' of $e1$ for Object b vs SNR',fontsize=fontsize); plt.ylabel('Residuals of $e1_b$',fontsize=fontsize)
plt.xlabel('SNR',fontsize=fontsize)
max_pts = np.max([np.mean(Resid_SNR_100[:,col]),np.mean(Resid_SNR_40[:,col]),np.mean(Resid_SNR_30[:,col]),
                  np.mean(Resid_SNR_20[:,col]),np.mean(Resid_SNR_15[:,col]),np.mean(Resid_SNR_10[:,col]),np.mean(Resid_SNR_5[:,col])])
min_pts = np.min([np.mean(Resid_SNR_100[:,col]),np.mean(Resid_SNR_40[:,col]),np.mean(Resid_SNR_30[:,col]),
                  np.mean(Resid_SNR_20[:,col]),np.mean(Resid_SNR_15[:,col]),np.mean(Resid_SNR_10[:,col]),np.mean(Resid_SNR_5[:,col])])
max_std = (np.max([np.std(Resid_SNR_100[:,col]),np.std(Resid_SNR_40[:,col]),np.std(Resid_SNR_30[:,col]),
                  np.std(Resid_SNR_20[:,col]),np.std(Resid_SNR_15[:,col]),np.std(Resid_SNR_10[:,col]),np.std(Resid_SNR_5[:,col])]))/np.sqrt(num_trials)

plt.ylim([min_pts-max_std,max_pts+max_std])

if data_pts:
    plt.scatter(100*np.ones(num_trials),Resid_SNR_100[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([100],np.mean(Resid_SNR_100[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([100],np.mean(Resid_SNR_100[:,col]),yerr=np.std(Resid_SNR_100[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')

if data_pts:
    plt.scatter(40*np.ones(num_trials),Resid_SNR_40[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([40],np.mean(Resid_SNR_40[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([40],np.mean(Resid_SNR_40[:,col]),yerr=np.std(Resid_SNR_40[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')

if data_pts:    
    plt.scatter(30*np.ones(num_trials),Resid_SNR_30[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([30],np.mean(Resid_SNR_30[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([30],np.mean(Resid_SNR_30[:,col]),yerr=np.std(Resid_SNR_30[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')

if data_pts:
    plt.scatter(20*np.ones(num_trials),Resid_SNR_20[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([20],np.mean(Resid_SNR_20[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([20],np.mean(Resid_SNR_20[:,col]),yerr=np.std(Resid_SNR_20[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


if data_pts:    
    plt.scatter(15*np.ones(num_trials),Resid_SNR_15[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([15],np.mean(Resid_SNR_15[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([15],np.mean(Resid_SNR_15[:,col]),yerr=np.std(Resid_SNR_15[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


if data_pts:
    plt.scatter(10*np.ones(num_trials),Resid_SNR_10[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([10],np.mean(Resid_SNR_10[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([10],np.mean(Resid_SNR_10[:,col]),yerr=np.std(Resid_SNR_10[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')

if data_pts:
    plt.scatter(5*np.ones(num_trials),Resid_SNR_10[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([5],np.mean(Resid_SNR_5[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([5],np.mean(Resid_SNR_5[:,col]),yerr=np.std(Resid_SNR_5[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')

ax4 = fig.add_subplot(gs[1,1])
col = 3 + 6
mark = 'x'
plt.title(title + ' of $e2$ for Object b vs SNR',fontsize=fontsize); plt.ylabel('Residuals of $e2_b$',fontsize=fontsize)
plt.xlabel('SNR',fontsize=fontsize)
max_pts = np.max([np.mean(Resid_SNR_100[:,col]),np.mean(Resid_SNR_40[:,col]),np.mean(Resid_SNR_30[:,col]),
                  np.mean(Resid_SNR_20[:,col]),np.mean(Resid_SNR_15[:,col]),np.mean(Resid_SNR_10[:,col]),np.mean(Resid_SNR_5[:,col])])
min_pts = np.min([np.mean(Resid_SNR_100[:,col]),np.mean(Resid_SNR_40[:,col]),np.mean(Resid_SNR_30[:,col]),
                  np.mean(Resid_SNR_20[:,col]),np.mean(Resid_SNR_15[:,col]),np.mean(Resid_SNR_10[:,col]),np.mean(Resid_SNR_5[:,col])])
max_std = (np.max([np.std(Resid_SNR_100[:,col]),np.std(Resid_SNR_40[:,col]),np.std(Resid_SNR_30[:,col]),
                  np.std(Resid_SNR_20[:,col]),np.std(Resid_SNR_15[:,col]),np.std(Resid_SNR_10[:,col]),np.std(Resid_SNR_5[:,col])]))/np.sqrt(num_trials)

plt.ylim([min_pts-max_std,max_pts+max_std])

if data_pts:
    plt.scatter(100*np.ones(num_trials),Resid_SNR_100[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([100],np.mean(Resid_SNR_100[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([100],np.mean(Resid_SNR_100[:,col]),yerr=np.std(Resid_SNR_100[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')

if data_pts:
    plt.scatter(40*np.ones(num_trials),Resid_SNR_40[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([40],np.mean(Resid_SNR_40[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([40],np.mean(Resid_SNR_40[:,col]),yerr=np.std(Resid_SNR_40[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')

if data_pts:    
    plt.scatter(30*np.ones(num_trials),Resid_SNR_30[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([30],np.mean(Resid_SNR_30[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([30],np.mean(Resid_SNR_30[:,col]),yerr=np.std(Resid_SNR_30[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')

if data_pts:
    plt.scatter(20*np.ones(num_trials),Resid_SNR_20[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([20],np.mean(Resid_SNR_20[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([20],np.mean(Resid_SNR_20[:,col]),yerr=np.std(Resid_SNR_20[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


if data_pts:    
    plt.scatter(15*np.ones(num_trials),Resid_SNR_15[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([15],np.mean(Resid_SNR_15[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([15],np.mean(Resid_SNR_15[:,col]),yerr=np.std(Resid_SNR_15[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


if data_pts:
    plt.scatter(10*np.ones(num_trials),Resid_SNR_10[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([10],np.mean(Resid_SNR_10[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([10],np.mean(Resid_SNR_10[:,col]),yerr=np.std(Resid_SNR_10[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')

if data_pts:
    plt.scatter(5*np.ones(num_trials),Resid_SNR_10[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([5],np.mean(Resid_SNR_5[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([5],np.mean(Resid_SNR_5[:,col]),yerr=np.std(Resid_SNR_5[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')



# Plotting for HLR and Flux for a and b----------------------------------------
fig = plt.figure(figsize=(20,11))
plt.suptitle(suptitle,fontsize=fontsize)
ax1 = fig.add_subplot(gs[0,0])
col = 0
mark = 'o'
plt.title(title + ' of $Flux$ for Object a vs SNR',fontsize=fontsize); plt.ylabel('Residuals of $Flux_a$',fontsize=fontsize)
plt.xlabel('SNR',fontsize=fontsize)
max_pts = np.max([np.mean(Resid_SNR_100[:,col]),np.mean(Resid_SNR_40[:,col]),np.mean(Resid_SNR_30[:,col]),
                  np.mean(Resid_SNR_20[:,col]),np.mean(Resid_SNR_15[:,col]),np.mean(Resid_SNR_10[:,col]),np.mean(Resid_SNR_5[:,col])])
min_pts = np.min([np.mean(Resid_SNR_100[:,col]),np.mean(Resid_SNR_40[:,col]),np.mean(Resid_SNR_30[:,col]),
                  np.mean(Resid_SNR_20[:,col]),np.mean(Resid_SNR_15[:,col]),np.mean(Resid_SNR_10[:,col]),np.mean(Resid_SNR_5[:,col])])
max_std = (np.max([np.std(Resid_SNR_100[:,col]),np.std(Resid_SNR_40[:,col]),np.std(Resid_SNR_30[:,col]),
                  np.std(Resid_SNR_20[:,col]),np.std(Resid_SNR_15[:,col]),np.std(Resid_SNR_10[:,col]),np.std(Resid_SNR_5[:,col])]))/np.sqrt(num_trials)

plt.ylim([min_pts-max_std,max_pts+max_std])

if data_pts:
    plt.scatter(100*np.ones(num_trials),Resid_SNR_100[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([100],np.mean(Resid_SNR_100[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([100],np.mean(Resid_SNR_100[:,col]),yerr=np.std(Resid_SNR_100[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')

if data_pts:
    plt.scatter(40*np.ones(num_trials),Resid_SNR_40[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([40],np.mean(Resid_SNR_40[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([40],np.mean(Resid_SNR_40[:,col]),yerr=np.std(Resid_SNR_40[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')

if data_pts:    
    plt.scatter(30*np.ones(num_trials),Resid_SNR_30[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([30],np.mean(Resid_SNR_30[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([30],np.mean(Resid_SNR_30[:,col]),yerr=np.std(Resid_SNR_30[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')

if data_pts:
    plt.scatter(20*np.ones(num_trials),Resid_SNR_20[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([20],np.mean(Resid_SNR_20[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([20],np.mean(Resid_SNR_20[:,col]),yerr=np.std(Resid_SNR_20[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


if data_pts:    
    plt.scatter(15*np.ones(num_trials),Resid_SNR_15[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([15],np.mean(Resid_SNR_15[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([15],np.mean(Resid_SNR_15[:,col]),yerr=np.std(Resid_SNR_15[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


if data_pts:
    plt.scatter(10*np.ones(num_trials),Resid_SNR_10[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([10],np.mean(Resid_SNR_10[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([10],np.mean(Resid_SNR_10[:,col]),yerr=np.std(Resid_SNR_10[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')

if data_pts:
    plt.scatter(5*np.ones(num_trials),Resid_SNR_10[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([5],np.mean(Resid_SNR_5[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([5],np.mean(Resid_SNR_5[:,col]),yerr=np.std(Resid_SNR_5[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


ax2 = fig.add_subplot(gs[0,1])
col = 1
mark = 'x'

plt.title(title + ' of $Hlr$ for Object a vs SNR',fontsize=fontsize); plt.ylabel('Residuals of $Hlr_a$',fontsize=fontsize)
plt.xlabel('SNR',fontsize=fontsize)
max_pts = np.max([np.mean(Resid_SNR_100[:,col]),np.mean(Resid_SNR_40[:,col]),np.mean(Resid_SNR_30[:,col]),
                  np.mean(Resid_SNR_20[:,col]),np.mean(Resid_SNR_15[:,col]),np.mean(Resid_SNR_10[:,col]),np.mean(Resid_SNR_5[:,col])])
min_pts = np.min([np.mean(Resid_SNR_100[:,col]),np.mean(Resid_SNR_40[:,col]),np.mean(Resid_SNR_30[:,col]),
                  np.mean(Resid_SNR_20[:,col]),np.mean(Resid_SNR_15[:,col]),np.mean(Resid_SNR_10[:,col]),np.mean(Resid_SNR_5[:,col])])
max_std = (np.max([np.std(Resid_SNR_100[:,col]),np.std(Resid_SNR_40[:,col]),np.std(Resid_SNR_30[:,col]),
                  np.std(Resid_SNR_20[:,col]),np.std(Resid_SNR_15[:,col]),np.std(Resid_SNR_10[:,col]),np.std(Resid_SNR_5[:,col])]))/np.sqrt(num_trials)

plt.ylim([min_pts-max_std,max_pts+max_std])

if data_pts:
    plt.scatter(100*np.ones(num_trials),Resid_SNR_100[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([100],np.mean(Resid_SNR_100[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([100],np.mean(Resid_SNR_100[:,col]),yerr=np.std(Resid_SNR_100[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')

if data_pts:
    plt.scatter(40*np.ones(num_trials),Resid_SNR_40[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([40],np.mean(Resid_SNR_40[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([40],np.mean(Resid_SNR_40[:,col]),yerr=np.std(Resid_SNR_40[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')

if data_pts:    
    plt.scatter(30*np.ones(num_trials),Resid_SNR_30[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([30],np.mean(Resid_SNR_30[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([30],np.mean(Resid_SNR_30[:,col]),yerr=np.std(Resid_SNR_30[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')

if data_pts:
    plt.scatter(20*np.ones(num_trials),Resid_SNR_20[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([20],np.mean(Resid_SNR_20[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([20],np.mean(Resid_SNR_20[:,col]),yerr=np.std(Resid_SNR_20[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


if data_pts:    
    plt.scatter(15*np.ones(num_trials),Resid_SNR_15[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([15],np.mean(Resid_SNR_15[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([15],np.mean(Resid_SNR_15[:,col]),yerr=np.std(Resid_SNR_15[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


if data_pts:
    plt.scatter(10*np.ones(num_trials),Resid_SNR_10[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([10],np.mean(Resid_SNR_10[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([10],np.mean(Resid_SNR_10[:,col]),yerr=np.std(Resid_SNR_10[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')

if data_pts:
    plt.scatter(5*np.ones(num_trials),Resid_SNR_10[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([5],np.mean(Resid_SNR_5[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([5],np.mean(Resid_SNR_5[:,col]),yerr=np.std(Resid_SNR_5[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


ax3 = fig.add_subplot(gs[1,0])
col = 0 + 6
mark = 'o'

plt.title(title + ' of $Flux$ for Object b vs SNR',fontsize=fontsize); plt.ylabel('Residuals of $Flux_b$',fontsize=fontsize)
plt.xlabel('SNR',fontsize=fontsize)
max_pts = np.max([np.mean(Resid_SNR_100[:,col]),np.mean(Resid_SNR_40[:,col]),np.mean(Resid_SNR_30[:,col]),
                  np.mean(Resid_SNR_20[:,col]),np.mean(Resid_SNR_15[:,col]),np.mean(Resid_SNR_10[:,col]),np.mean(Resid_SNR_5[:,col])])
min_pts = np.min([np.mean(Resid_SNR_100[:,col]),np.mean(Resid_SNR_40[:,col]),np.mean(Resid_SNR_30[:,col]),
                  np.mean(Resid_SNR_20[:,col]),np.mean(Resid_SNR_15[:,col]),np.mean(Resid_SNR_10[:,col]),np.mean(Resid_SNR_5[:,col])])
max_std = (np.max([np.std(Resid_SNR_100[:,col]),np.std(Resid_SNR_40[:,col]),np.std(Resid_SNR_30[:,col]),
                  np.std(Resid_SNR_20[:,col]),np.std(Resid_SNR_15[:,col]),np.std(Resid_SNR_10[:,col]),np.std(Resid_SNR_5[:,col])]))/np.sqrt(num_trials)

plt.ylim([min_pts-max_std,max_pts+max_std])

if data_pts:
    plt.scatter(100*np.ones(num_trials),Resid_SNR_100[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([100],np.mean(Resid_SNR_100[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([100],np.mean(Resid_SNR_100[:,col]),yerr=np.std(Resid_SNR_100[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')

if data_pts:
    plt.scatter(40*np.ones(num_trials),Resid_SNR_40[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([40],np.mean(Resid_SNR_40[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([40],np.mean(Resid_SNR_40[:,col]),yerr=np.std(Resid_SNR_40[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')

if data_pts:    
    plt.scatter(30*np.ones(num_trials),Resid_SNR_30[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([30],np.mean(Resid_SNR_30[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([30],np.mean(Resid_SNR_30[:,col]),yerr=np.std(Resid_SNR_30[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')

if data_pts:
    plt.scatter(20*np.ones(num_trials),Resid_SNR_20[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([20],np.mean(Resid_SNR_20[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([20],np.mean(Resid_SNR_20[:,col]),yerr=np.std(Resid_SNR_20[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


if data_pts:    
    plt.scatter(15*np.ones(num_trials),Resid_SNR_15[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([15],np.mean(Resid_SNR_15[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([15],np.mean(Resid_SNR_15[:,col]),yerr=np.std(Resid_SNR_15[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


if data_pts:
    plt.scatter(10*np.ones(num_trials),Resid_SNR_10[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([10],np.mean(Resid_SNR_10[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([10],np.mean(Resid_SNR_10[:,col]),yerr=np.std(Resid_SNR_10[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')

if data_pts:
    plt.scatter(5*np.ones(num_trials),Resid_SNR_10[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([5],np.mean(Resid_SNR_5[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([5],np.mean(Resid_SNR_5[:,col]),yerr=np.std(Resid_SNR_5[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')



ax4 = fig.add_subplot(gs[1,1])
col = 1 + 6
mark = 'x'

plt.title(title + ' of $Hlr$ for Object b vs SNR',fontsize=fontsize); plt.ylabel('Residuals of $Hlr_b$',fontsize=fontsize)
plt.xlabel('SNR',fontsize=fontsize)
max_pts = np.max([np.mean(Resid_SNR_100[:,col]),np.mean(Resid_SNR_40[:,col]),np.mean(Resid_SNR_30[:,col]),
                  np.mean(Resid_SNR_20[:,col]),np.mean(Resid_SNR_15[:,col]),np.mean(Resid_SNR_10[:,col]),np.mean(Resid_SNR_5[:,col])])
min_pts = np.min([np.mean(Resid_SNR_100[:,col]),np.mean(Resid_SNR_40[:,col]),np.mean(Resid_SNR_30[:,col]),
                  np.mean(Resid_SNR_20[:,col]),np.mean(Resid_SNR_15[:,col]),np.mean(Resid_SNR_10[:,col]),np.mean(Resid_SNR_5[:,col])])
max_std = (np.max([np.std(Resid_SNR_100[:,col]),np.std(Resid_SNR_40[:,col]),np.std(Resid_SNR_30[:,col]),
                  np.std(Resid_SNR_20[:,col]),np.std(Resid_SNR_15[:,col]),np.std(Resid_SNR_10[:,col]),np.std(Resid_SNR_5[:,col])]))/np.sqrt(num_trials)

plt.ylim([min_pts-max_std,max_pts+max_std])

if data_pts:
    plt.scatter(100*np.ones(num_trials),Resid_SNR_100[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([100],np.mean(Resid_SNR_100[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([100],np.mean(Resid_SNR_100[:,col]),yerr=np.std(Resid_SNR_100[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')

if data_pts:
    plt.scatter(40*np.ones(num_trials),Resid_SNR_40[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([40],np.mean(Resid_SNR_40[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([40],np.mean(Resid_SNR_40[:,col]),yerr=np.std(Resid_SNR_40[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')

if data_pts:    
    plt.scatter(30*np.ones(num_trials),Resid_SNR_30[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([30],np.mean(Resid_SNR_30[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([30],np.mean(Resid_SNR_30[:,col]),yerr=np.std(Resid_SNR_30[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')

if data_pts:
    plt.scatter(20*np.ones(num_trials),Resid_SNR_20[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([20],np.mean(Resid_SNR_20[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([20],np.mean(Resid_SNR_20[:,col]),yerr=np.std(Resid_SNR_20[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


if data_pts:    
    plt.scatter(15*np.ones(num_trials),Resid_SNR_15[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([15],np.mean(Resid_SNR_15[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([15],np.mean(Resid_SNR_15[:,col]),yerr=np.std(Resid_SNR_15[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')


if data_pts:
    plt.scatter(10*np.ones(num_trials),Resid_SNR_10[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([10],np.mean(Resid_SNR_10[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([10],np.mean(Resid_SNR_10[:,col]),yerr=np.std(Resid_SNR_10[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')

if data_pts:
    plt.scatter(5*np.ones(num_trials),Resid_SNR_10[:,col],marker=mark,c='g',alpha=alpha)
plt.scatter([5],np.mean(Resid_SNR_5[:,col]),marker='o',c='b',linewidth=mean_linewidth)
plt.errorbar([5],np.mean(Resid_SNR_5[:,col]),yerr=np.std(Resid_SNR_5[:,col])/np.sqrt(num_trials),ecolor='k',elinewidth=bar_linewidth)
plt.axhline(0,color='k',linestyle='--')
