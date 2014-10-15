# -*- coding: utf-8 -*-
"""
Created on Wed Jul  9 14:58:05 2014

@author: luis
"""

# Based on galsim_First_Test by LCA
# Modifications by MSSG
# Start: 7/31/2014

## This file looks at the effects of simple galsim Sersic objects and methods.

import sys
import numpy as np
import galsim
import os
import matplotlib.pyplot as plt

# Parameters
gal_flux = 1e6     # total counts on the image
gal_HLR = 4.       # arcsec
psf_sigma = 1     # arcsec
pixel_scale = 0.2  # arcsec / pixel
noise = 3e1        # standard deviation of the counts in each pixel
pstampsize = 100

# Shear params
e1 = 0.4 
e2 = -0.5

# Booleans to affect final image
noise_on_image = True
sum_galaxies = True
shear_image = True
set_seed = True # If we want to set random num seed for reproducability


# Create primary galaxy 
# galaxy = galsim.Gaussian(flux=gal_flux*2.0,half_light_radius=gal_HLR)
deVauc_ix = 4 # deVauc bulge
expl_ix = 1   # expl bulge
galaxy = galsim.Sersic(n=deVauc_ix, flux=gal_flux*2.0,half_light_radius=gal_HLR)
galaxy += galsim.Sersic(n=expl_ix, flux=gal_flux*2.0,half_light_radius=gal_HLR)


# Create secondary galaxy and shift
fluxratio = 2 #  Primary to Secondary ratio
offset_B = 10 
galaxy_B = galsim.Gaussian(flux=gal_flux/fluxratio,half_light_radius=gal_HLR/2)
galaxy_B = galaxy_B.shift((offset_B ,offset_B )) # Shift by this many arcsec in a direction

# Set seed to non-zero integer value for reproducability
if set_seed:
    dev = galsim.BaseDeviate(213524)

# Produce a shear on the galaxy
if shear_image:
    galaxy = galaxy.shear(e1=e1,e2=e2)
    galaxy_B = galaxy_B.shear(e1=e1,e2=-e2)

# Sum the galaxies (concatenate the images)
if sum_galaxies:
    galaxy = galaxy + galaxy_B

# Create psf
psf = galsim.Gaussian(flux=1.,sigma=psf_sigma)

# Convolve
convolved_galaxy = galsim.Convolve([galaxy, psf])

## Draws Images of each for comparison
#image = galsim.ImageD(pstampsize,pstampsize,scale=pixel_scale)
#image_gal = galaxy.drawImage(image=image,method='phot',rng=dev)
#image_psf = psf.draw(image=image)
#image_convolved_galaxy = convolved_galaxy.drawImage(image=image,method='phot',rng=dev)

# Draws Images of each for comparison
image = galsim.ImageD(pstampsize,pstampsize,scale=pixel_scale)
image_gal = galaxy.drawImage(image=image,method='phot',rng=dev)
image_psf = psf.draw(image=image)
image_convolved_galaxy = convolved_galaxy.drawImage(image=image,method='phot',rng=dev)

# Add noise    
if noise_on_image:
    image_gal.addNoise(galsim.GaussianNoise(sigma=noise))
    image_psf.addNoise(galsim.GaussianNoise(sigma=noise))
    image_convolved_galaxy.addNoise(galsim.GaussianNoise(sigma=noise))

plt.imshow(image_convolved_galaxy.array,interpolation='none',origin='lower')
plt.colorbar()
plt.title('Blended galsim galaxy profiles.')
plt.xlabel('$Pixels$'); plt.ylabel('$Pixels$')
plt.show()

# Write the image to a file
if not os.path.isdir('output'):
    os.mkdir('output') # Make the dir if it doesn't already exist

file_name_gal = os.path.join('output','gal.fits')
file_name_psf = os.path.join('output','psf.fits')
file_name_convolved_galaxy = os.path.join('output','convolved_galaxy.fits')
# Note: if the file already exists, this will overwrite it.
image_gal.write(file_name_gal)
image_psf.write(file_name_psf)
image_convolved_galaxy.write(file_name_convolved_galaxy)
