# -*- coding: utf-8 -*-
"""
Created on Wed Jul  9 14:58:05 2014

@author: luis
"""

## This file looks at the effects of simple galsim objects and methods.

import sys
import numpy as np
import galsim
import os
import matplotlib.pyplot as plt

# Parameters
gal_flux = 1e5     # total counts on the image
gal_HLR = 8.       # arcsec
psf_sigma = 1.     # arcsec
pixel_scale = 0.2  # arcsec / pixel
noise = 2.        # standard deviation of the counts in each pixel

# Shear params
g1 = 0.5
g2 = 0.5

# Booleans to affect final image
noise_on_image = False
sum_of_galaxy = False
shear_on_image = True

# Create primary galaxy 
galaxy = galsim.Gaussian(flux=gal_flux,half_light_radius=gal_HLR)

# Create secondary galaxy and shift
galaxy_sec = galsim.Gaussian(flux=gal_flux/2.0,half_light_radius=gal_HLR)
galaxy_sec = galaxy_sec.shift((10,10))

# Sum the galaxies (concatenate the images)
if sum_of_galaxy:
    galaxy = galaxy + galaxy_sec

# Produce a shear on the galaxy
if shear_on_image:
    galaxy = galaxy.shear(g1=g1,g2=g2)

im = galaxy.drawShoot()
plt.imshow(im.array,interpolation='none')
plt.show()

# Create psf
psf = galsim.Gaussian(flux=1.,sigma=psf_sigma)

# Convolve
final = galsim.Convolve([galaxy, psf])

# Draws Images of each for comparison
image_gal = galaxy.drawShoot(scale=pixel_scale)
image_psf = psf.draw(scale=pixel_scale)
image_final = final.drawShoot(scale=pixel_scale)

# Add noise    
if noise_on_image:
    image_gal.addNoise(galsim.GaussianNoise(sigma=noise))
    image_psf.addNoise(galsim.GaussianNoise(sigma=noise))
    image_final.addNoise(galsim.GaussianNoise(sigma=noise))

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
