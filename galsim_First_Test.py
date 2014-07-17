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
gal_HLR = 4.       # arcsec
psf_sigma = 1.     # arcsec
pixel_scale = 0.2  # arcsec / pixel
noise = 1.        # standard deviation of the counts in each pixel
size = 300

# Shear params
e1 = 0.5
e2 = 0.5

# Booleans to affect final image
noise_on_image = True
sum_of_galaxy = True
shear_on_image = True

# Create primary galaxy 
galaxy = galsim.Gaussian(flux=gal_flux*2.0,half_light_radius=gal_HLR)

# Create secondary galaxy and shift
galaxy_sec = galsim.Gaussian(flux=gal_flux/2.5,half_light_radius=gal_HLR/2)
galaxy_sec = galaxy_sec.shift((5,5))

# Produce a shear on the galaxy
if shear_on_image:
    galaxy = galaxy.shear(e1=e1,e2=e2)
    galaxy_sec = galaxy_sec.shear(e1=e1,e2=-e2)

# Sum the galaxies (concatenate the images)
if sum_of_galaxy:
    galaxy = galaxy + galaxy_sec

# Create psf
psf = galsim.Gaussian(flux=1.,sigma=psf_sigma)

# Convolve
final = galsim.Convolve([galaxy, psf])

# Draws Images of each for comparison
image = galsim.ImageD(size,size,scale=pixel_scale)
image_gal = galaxy.drawShoot(image=image)
image_psf = psf.draw(image=image)
image_final = final.drawShoot(image=image)

# Add noise    
if noise_on_image:
    image_gal.addNoise(galsim.GaussianNoise(sigma=noise))
    image_psf.addNoise(galsim.GaussianNoise(sigma=noise))
    image_final.addNoise(galsim.GaussianNoise(sigma=noise))

plt.imshow(image_final.array,interpolation='none',origin='lower')
plt.colorbar()
plt.title('Blended galsim galaxy profiles.')
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
