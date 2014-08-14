# -*- coding: utf-8 -*-
"""
Created on Wed Jul  9 14:58:05 2014

@author: luis
"""

## This file looks at the effects of simple galsim objects and methods.

from __future__ import division
import sys
import numpy as np
import galsim
import os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# Shear params
e1 = 0.1
e2 = 0.3

texp = 6900       # sec
sky_noise = np.sqrt(26.8*texp)    # e/pixel

# Parameters
gal_flux = 400000     # total counts on the image
gal_HLR = 4.       # arcsec
psf_sigma = 1.     # arcsec
pixel_scale = 0.2  # arcsec / pixel
size = 200

# Booleans to affect final image
noise_on_image = True
sum_of_galaxy = True
shear_on_image = True
set_seed = True

# Create primary galaxy 
galaxy = galsim.Gaussian(flux=gal_flux*5.0,half_light_radius=gal_HLR).shift((0,3))

# Create secondary galaxy and shift
galaxy_sec = galsim.Gaussian(flux=gal_flux*5,half_light_radius=gal_HLR)
galaxy_sec = galaxy_sec.shift((-5,-7))

# Set seed to non-zero integer value for reproducability
if set_seed:
    dev = galsim.BaseDeviate(213524)

# Produce a shear on the galaxy
if shear_on_image:
    galaxy = galaxy.shear(g1=e1,g2=e2)
    galaxy_sec = galaxy_sec.shear(g1=e1,g2=-e2)

# Sum the galaxies (concatenate the images)
if sum_of_galaxy:
    galaxy = galaxy + galaxy_sec

# Create psf
psf = galsim.Gaussian(flux=1.,sigma=psf_sigma)

# Convolve
final = galsim.Convolve([galaxy, psf])

## Draws Images of each for comparison
image = galsim.ImageD(size,size,scale=pixel_scale)
image_gal = galaxy.drawImage(image=image,method='fft')
image_psf = psf.draw(image=image)
image_final = final.drawImage(image=image,method='fft')
image_final_w_noise = image_final

# Add noise    
if noise_on_image:
    image_psf.addNoise(galsim.PoissonNoise(rng=dev,sky_level=sky_noise))
    image_final_w_noise.addNoise(galsim.PoissonNoise(rng=dev,sky_level=sky_noise))

gs = gridspec.GridSpec(1,2)
fig = plt.figure(figsize=(20,11))
a = plt.imshow(image_final.array,origin='lower')
ax1 = fig.add_subplot(111)
plt.title('Blended Galaxy Profiles w/ Sky Noise.',fontsize=20)
plt.xlabel('$Pixels$',fontsize=20); plt.ylabel('$Pixels$',fontsize=20)
plt.colorbar(a,shrink=1)
ax2 = fig.add_subplot(gs[0,1])
b = plt.imshow(image_gal.array,interpolation='none',origin='lower')
plt.colorbar(b,shrink=1)
plt.title('Blended Galsim Galaxy Profiles w/ Noise.')
plt.xlabel('$Pixels$'); plt.ylabel('$Pixels$')
plt.show()

# Write the image to a file
if not os.path.isdir('output'):
    os.mkdir('output')
file_name_gal = os.path.join('output','gal.fits')
file_name_psf = os.path.join('output','psf.fits')
file_name_final = os.path.join('output','final.fits')
# Note: if the file already exists, this will overwrite it.
image_gal.write(file_name_gal)
image_psf.write(file_name_psf)
image_final.write(file_name_final)
