# -*- coding: utf-8 -*-
"""
Created on Wed Aug 13 13:54:38 2014

@author: luis
"""

## This script runs through a simple double fit of two sersic profiles

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

# Parameters for object a
flux_a = 5e5          # total counts on the image
hlr_a = 1             # arcsec
e1_a = 0.0
e2_a = 0.0
x0_a = -1
y0_a = 0
n_a = 0.5

# Parameters for object b
flux_b = 5e5          # total counts on the image
hlr_b = hlr_a         # arcsec
e1_b = 0.0
e2_b = 0.0
x0_b = 1
y0_b = 0
n_b = 0.5

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
Sbar = 26.8 # sky photons per second per pixel
sky_level = 0
sky_noise = np.sqrt(sky_level)

# psf properties
psf_flag = False
beta = 3
fwhm_psf = 0.6

# Obtain instantiation
im, best_fit, result = noiseLibrary.run_2_galaxy_full_params_simple(flux_a,hlr_a,e1_a,e2_a,x0_a,y0_a,n_a,
                                                                    flux_b,hlr_b,e1_b,e2_b,x0_b,y0_b,n_b,
                                                                    psf_flag,beta,fwhm_psf,
                                                                    x_len,y_len,pixel_scale,sersic_func,sersic_func,seed_1,seed_2,seed_3,
                                                                    add_noise_flag,sky_level)

# Calculate the signal-to-noise ratio
threshold = 0.5*np.sqrt(sky_level)
mask = im.array > threshold

nu = np.sqrt(texp/Sbar)*np.sqrt((mask*im.array*im.array).sum())

threshold_j = 0.5*np.sqrt(Sbar * texp) / texp # threshold for contributing to SNR
weight = im.array
nu_j = (im.array * weight * mask * texp).sum() / np.sqrt((weight * weight * mask * texp * Sbar).sum())
                                                                    
# Report errors
lmfit.report_errors(result.params)

# Obtain the covariance and correlation matrix.
error_diag = np.sqrt(np.diag(result.covar))
error_mat = np.outer(error_diag,error_diag)
correlation_mat = result.covar/error_mat

# Provide the fontsize
fonts = 10
# Initialize the geometry of the grid
gs = gridspec.GridSpec(3,5)
# Set the figure size
fig = plt.figure(figsize=(15,10))
plt.suptitle(r'$Image\/Scale:$ $\frac{0.2\rq\rq}{Pixel}$',fontsize=fonts+10)
# Add the first subplot
ax1 = fig.add_subplot(gs[0,0])
plt.title('$Image$',fontsize=fonts+8)
plt.ylabel('$Pixels$',fontsize=fonts+8)
plt.tick_params(axis='both',which='major',labelsize=9)
a = ax1.imshow(im.array,interpolation='none',origin='lower',vmax=im.array.max(),vmin=im.array.min())
plt.colorbar(a,shrink=1)
# Add the second subplot
ax2 = fig.add_subplot(gs[1,0])
plt.ylabel('$Pixels$',fontsize=fonts+8)
plt.tick_params(axis='both',which='major',labelsize=9)
plt.title('$Best Fit$',fontsize=fonts+8)
b = ax2.imshow(best_fit.array,origin='lower',vmax=best_fit.array.max(),vmin=best_fit.array.min())
plt.colorbar(b,shrink=1)
# Add the third subplot
ax3 = fig.add_subplot(gs[2,0])
plt.tick_params(axis='both',which='major',labelsize=9)
plt.title('$Residual$',fontsize=fonts+8)
plt.xlabel('$Pixels$',fontsize=fonts+8)
plt.ylabel('$Pixels$',fontsize=fonts+8)
c = ax3.imshow((best_fit - im).array,interpolation='none',origin='lower')
plt.colorbar(c,shrink=1)
# Add the correlation matrix
ax4 = fig.add_subplot(gs[:,1:])
plt.title('Correlation Coefficient Matrix',fontsize=fonts+8)
plt.xticks([0,1,2,3,4,5,6,7,8,9,10,11],['$Flux_a$','$HLR_a$','$e1_a$','$e2_a$','$x0_a$','$y0_a$',
           '$Flux_b$','$HLR_b$','$e1_b$','$e2_b$','$x0_b$','$y0_b$'],fontsize=15)
plt.yticks([0,1,2,3,4,5,6,7,8,9,10,11],['$Flux_a$','$HLR_a$','$e1_a$','$e2_a$','$x0_a$','$y0_a$',
           '$Flux_b$','$HLR_b$','$e1_b$','$e2_b$','$x0_b$','$y0_b$'],fontsize=15)
d = ax4.imshow(np.around(correlation_mat,decimals=2),interpolation='none',origin='lower',vmin=-1,vmax=1)
plt.colorbar(d,shrink=1)
for (i,j), val in np.ndenumerate(correlation_mat):
    ax4.annotate('%0.2f'%(val), (j,i), ha='center', va='center',size=12)
text = ax4.text(0,-1.5,'Separation: %.2f arcseconds'%sep,fontsize=15)
plt.show()
    