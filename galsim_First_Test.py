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


# Parameters
gal_flux = 1e5     # total counts on the image
gal_HLR = 8.       # arcsec
psf_sigma = 1.     # arcsec
pixel_scale = 0.2  # arcsec / pixel
noise = 12.        # standard deviation of the counts in each pixel

g1 = 0.5
g2 = 0.1

noise_on_image = True
shear_on_image = True

# Create galaxy 
galaxy = galsim.Gaussian(flux=gal_flux,half_light_radius=gal_HLR)

# Produce a shear on the galaxy
if shear_on_image:
    galaxy = galaxy.shear(g1=g1,g2=g2)

# Create psf
psf = galsim.Gaussian(flux=1.,sigma=psf_sigma)

# Convolve
final = galsim.Convolve([galaxy, psf])

# Draws Images of each for comparison
image_gal = galaxy.draw(scale=pixel_scale)
image_psf = psf.draw(scale=pixel_scale)
image_final = final.draw(scale=pixel_scale)

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
