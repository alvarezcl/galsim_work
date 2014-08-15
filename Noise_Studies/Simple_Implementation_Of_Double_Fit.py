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
flux_a = 5e4          # total counts on the image
hlr_a = 1         # arcsec
e1_a = 0.0
e2_a = 0.0
x0_a = -2
y0_a = 0
n_a = 0.5

# Parameters for object b
flux_b = 5e4       # total counts on the image
hlr_b = hlr_a         # arcsec
e1_b = 0.0
e2_b = 0.0
x0_b = 2
y0_b = 0
n_b = 0.5

# Galsim function definitions
sersic_func = galsim.Sersic

# Set the RNG
seed_1 = galsim.BaseDeviate(1)
seed_2 = galsim.BaseDeviate(2)
seed_3 = galsim.BaseDeviate(3)

# Image properties
pixel_scale = 1/5     # arcsec / pixel
x_len = y_len = 100            # pixel
sky_level = 1000         # counts / pixel

# Obtain instantiation
im, best_fit, result = noiseLibrary.run_2_galaxy_full_params_simple(flux_a,hlr_a,e1_a,e2_a,x0_a,y0_a,n_a,
                                                                    flux_b,hlr_b,e1_b,e2_b,x0_b,y0_b,n_b,
                                                                    x_len,y_len,pixel_scale,sersic_func,sersic_func,seed_1,seed_2,seed_3,sky_level)


# Report errors
lmfit.report_errors(result.params)

# Provide the fontsize
fonts = 10
# Initialize the geometry of the grid
gs = gridspec.GridSpec(1,3)
# Set the figure size
fig = plt.figure(figsize=(15,10))
plt.suptitle(r'$Image\/Scale:$ $\frac{0.2\rq\rq}{Pixel}$',fontsize=fonts+10)
# Add the first subplot
ax1 = fig.add_subplot(gs[0,0])
plt.title('$Image$',fontsize=fonts+8)
plt.ylabel('$Pixels$',fontsize=fonts+8)
plt.tick_params(axis='both',which='major',labelsize=9)
a = ax1.imshow(im.array,interpolation='none',origin='lower',vmax=im.array.max(),vmin=im.array.min())
plt.colorbar(a,shrink=0.5)
# Add the second subplot
ax2 = fig.add_subplot(gs[0,1])
plt.ylabel('$Pixels$',fontsize=fonts+8)
plt.tick_params(axis='both',which='major',labelsize=9)
plt.title('$Best Fit$',fontsize=fonts+8)
b = ax2.imshow(best_fit.array,origin='lower',vmax=best_fit.array.max(),vmin=best_fit.array.min())
plt.colorbar(b,shrink=0.5)
# Add the third subplot
ax3 = fig.add_subplot(gs[0,2])
plt.tick_params(axis='both',which='major',labelsize=9)
plt.title('$Residual$',fontsize=fonts+8)
plt.xlabel('$Pixels$',fontsize=fonts+8)
plt.ylabel('$Pixels$',fontsize=fonts+8)
c = ax3.imshow((best_fit - im).array,interpolation='none',origin='lower')
plt.colorbar(c,shrink=0.5)