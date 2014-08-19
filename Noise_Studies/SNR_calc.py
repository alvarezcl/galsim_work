# -*- coding: utf-8 -*-
"""
Created on Tue Aug 19 10:08:16 2014

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
flux_a = 0          # total counts on the image
hlr_a = 1             # arcsec
e1_a = 0.0
e2_a = 0.0
x0_a = -1
y0_a = 0
n_a = 0.5

# Parameters for object b
flux_b = flux_a          # total counts on the image
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
add_noise_flag = False
texp = 6900 # seconds
sbar = 26.8 # sky photons per second per pixel
sky_level = texp*sbar
sky_noise = np.sqrt(sky_level)

# psf properties
psf_flag = False
beta = 3
fwhm_psf = 0.6

data = {'Flux_tot':[],'SNR':[],'Frac_pix':[]}

for i in xrange(0,20):
    flux_a += (i+1)*1000
    flux_b = flux_a
    data['Flux_tot'].append(flux_a+flux_b)
    # Obtain instantiation
    im = noiseLibrary.draw_simple(flux_a,hlr_a,e1_a,e2_a,x0_a,y0_a,n_a,
                                  flux_b,hlr_b,e1_b,e2_b,x0_b,y0_b,n_b,
                                  psf_flag,beta,fwhm_psf,
                                  x_len,y_len,pixel_scale,sersic_func,sersic_func,seed_1,seed_2,seed_3,
                                  add_noise_flag,sky_level)
                                  
    nu_s, mask_per_s, nu, mask = noiseLibrary.calc_SNR(im, texp, sbar, 0.5)
    data['SNR'].append(nu)
    pix_count_image = (im.array > 0).sum()
    pix_count_masked_image = (mask > 0).sum()
    fractional_pix_count = pix_count_masked_image/pix_count_image 
    data['Frac_pix'].append(fractional_pix_count)
    print nu
    print i                            

plt.figure(), plt.imshow(im.array,interpolation='none',origin='lower'); plt.title('Image'); plt.colorbar()
plt.figure(), plt.imshow(mask,interpolation='none',origin='lower'); plt.title('Mask'); plt.colorbar()
plt.figure(), plt.imshow(mask*im.array,interpolation='none',origin='lower',vmin=im.array.min(),vmax=im.array.max()); plt.title('Masked Imaged'); plt.colorbar()
plt.show()

domain = (data.pop('Flux_tot'))
figsize = (20,11)
pos = [(0,0),(1,0)]
row = 2
col = 1
colors = ['g','b']
markers = ['x','o']
x_label = '$\/Total\/Flux$'
x_lim = [np.min(domain),np.max(domain)]
y_lim = [[np.min(data['Frac_pix']),np.max(data['Frac_pix'])],[np.min(data['SNR']),np.max(data['SNR'])]]
y_label = ['$Pixel\/Fraction$','$SNR$']
suptitle = '$Varying\/Total\/Flux;\/Sep:\/%.2f\/arcsec$'%(sep)
t_fontsize = 18
fontsize = 15
l_fontsize = 10
text_p = (0,-2,'$Separation:\/%.2f\/ Hlr:\/%.2f$'%(sep,hlr_a),fontsize)
fig = noiseLibrary.plot(domain,data,figsize,pos,row,col,colors,markers,
                        legend=None,fontsize=fontsize,legend_fontsize=l_fontsize,t_fontsize=t_fontsize,
                        text_x=text_p[0],text_y=text_p[1],text_s='',text_fs=text_p[3],
                        suptitle=suptitle,x_label=x_label,x_lim=x_lim,y_lim=None,y_label=y_label)
                                                                                             
                                                                    
                                                                    