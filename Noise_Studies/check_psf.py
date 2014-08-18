# -*- coding: utf-8 -*-
"""
Created on Fri Aug 15 17:45:36 2014

@author: luis
"""

import noiseLibrary
import galsim
import numpy as np

# Parameters for object a
flux_a_tot = 5e6          # total counts on the image
bulge_a_frac = 0.5
hlr_a_bulge = 1
hlr_a_disk = 1        # arcsec
e1_a = 0.0
e2_a = 0.0
x0_a = -2
y0_a = 0
n_a_bulge = 3.5
n_a_disk = 1

# Parameters for object b
flux_b_tot = 5e4       # total counts on the image
bulge_b_frac = 0.5
hlr_b_bulge = hlr_a_bulge
hlr_b_disk = hlr_a_disk        # arcsec
e1_b = 0.0
e2_b = 0.0
x0_b = 2
y0_b = 0
n_b_bulge = 3.5
n_b_disk = 1

# psf properties
psf_flag = True
beta = 3
fwhm_psf = 0.6

# Galsim function definitions
sersic_func = galsim.Sersic

# Set the RNG
seed_1 = galsim.BaseDeviate(1)
seed_2 = galsim.BaseDeviate(2)
seed_3 = galsim.BaseDeviate(3)

# Image properties
pixel_scale = 1/5     # arcsec / pixel
x_len = y_len = 100            # pixel
sky_level = 0         # counts / pixel

gal_a_bulge = galsim.Sersic(n=n_b_bulge,flux=flux_a_tot*bulge_b_frac,half_light_radius=hlr_b_bulge)
gal_a_disk = galsim.Sersic(n=n_b_disk,flux=flux_b_tot*(1-bulge_b_frac),half_light_radius=hlr_b_disk)

gal_b_bulge = galsim.Sersic(n=n_b_bulge,flux=flux_b_tot*bulge_b_frac,half_light_radius=hlr_b_bulge)
gal_b_disk = galsim.Sersic(n=n_b_disk,flux=flux_b_tot*(1-bulge_b_frac),half_light_radius=hlr_b_disk)

gal = gal_a_bulge + gal_a_disk + gal_b_bulge + gal_b_disk

psf = galsim.Moffat(beta=beta,fwhm=fwhm_psf)

gal_psf = galsim.Convolve([gal,psf])

image_canvas = galsim.ImageD(x_len,y_len,scale=pixel_scale)
image_a = gal_psf.drawImage(image=image_canvas,method='fft')

# Then for the other version

gal_psf_sum = noiseLibrary.convolve_with_psf(gal_a_bulge,beta,fwhm_psf) + noiseLibrary.convolve_with_psf(gal_a_disk,beta,fwhm_psf) + noiseLibrary.convolve_with_psf(gal_b_bulge,beta,fwhm_psf) + noiseLibrary.convolve_with_psf(gal_b_disk,beta,fwhm_psf)
image_b = gal_psf_sum.drawImage(image=image_canvas,method='fft')

