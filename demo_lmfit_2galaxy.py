# -*- coding: utf-8 -*-
"""
Created on Fri Jul 11 16:32:06 2014

@author: luis
"""

## This file creates the functions necessary for usage with
## lmfit and obtains the estimates of HLR,Flux,e1,e2,x0, and
## y0 for a 2 object image.

import lmfit
import galsim
import numpy as np
import matplotlib.pyplot as plt

# The "data" aspect that throws photons at the image for two galaxies
def drawShoot_galaxy_2(flux_1,hlr_1,e1_1,e2_1,x_center1,y_center1,
                       flux_2,hlr_2,e1_2,e2_2,x_center2,y_center2,
                       x_len,y_len,scale):
                           
    gal_1 = galsim.Gaussian(half_light_radius=hlr_1, flux=flux_1)
    gal_1 = gal_1.shear(e1=e1_1, e2=e2_1)
    gal_1 = gal_1.shift(x_center1,y_center1)
    image_1 = galsim.ImageD(x_len, y_len, scale=scale)
    image_1 = gal_1.drawShoot(image=image_1)
 
    gal_2 = galsim.Gaussian(half_light_radius=hlr_2, flux=flux_2)
    gal_2 = gal_2.shear(e1=e1_2, e2=e2_2)
    gal_2 = gal_2.shift(x_center2,y_center2)
    image_2 = galsim.ImageD(x_len, y_len, scale=scale)
    image_2 = gal_2.drawShoot(image=image_2)
    image = image_1 + image_2
    
    return image

# Continous version, or what you would call the "model" for two galaxies
def draw_galaxy_2(flux_1,hlr_1,e1_1,e2_1,x_center1,y_center1,
                  flux_2,hlr_2,e1_2,e2_2,x_center2,y_center2,
                  x_len,y_len,scale):
                           
    gal_1 = galsim.Gaussian(half_light_radius=hlr_1, flux=flux_1)
    gal_1 = gal_1.shear(e1=e1_1, e2=e2_1)
    gal_1 = gal_1.shift(x_center1,y_center1)
    image_1 = galsim.ImageD(x_len, y_len, scale=scale)
    image_1 = gal_1.draw(image=image_1)
 
    gal_2 = galsim.Gaussian(half_light_radius=hlr_2, flux=flux_2)
    gal_2 = gal_2.shear(e1=e1_2, e2=e2_2)
    gal_2 = gal_2.shift(x_center2,y_center2)
    image_2 = galsim.ImageD(x_len, y_len, scale=scale)
    image_2 = gal_2.draw(image=image_2)
    image = image_1 + image_2
    
    return image

# Take the difference of the data and the model for two galaxies
def resid_2(param, target_image,x_len, y_len, scale):
    flux_1 = param['flux_1'].value
    hlr_1 = param['hlr_1'].value
    e1_1 = param['e1_1'].value
    e2_1 = param['e2_1'].value
    x_center1 = param['x_center1'].value
    y_center1 = param['y_center1'].value

    flux_2 = param['flux_2'].value
    hlr_2 = param['hlr_2'].value
    e1_2 = param['e1_2'].value
    e2_2 = param['e2_2'].value
    x_center2 = param['x_center2'].value
    y_center2 = param['y_center2'].value
    
    image = draw_galaxy_2(flux_1,hlr_1,e1_1,e2_1,x_center1,y_center1,
                  flux_2,hlr_2,e1_2,e2_2,x_center2,y_center2,
                  x_len,y_len,scale)
                  
    return (image-target_image).array.ravel()

# Define independent parameters for image not trying to minimize
x_len = 50
y_len = 50
scale = 0.2

# Draw from some true distribution of points with flux, hlr, e1, e2
# x0, y0, and picture parameters.
im = drawShoot_galaxy_2(10000, 1.0, 0.5, 0.3, 0, 0,
                        3000, 0.5, 0.5, -0.5, 2.0, 2.0, 
                        x_len, y_len, scale)

# Define some seed that's far from true values and insert into
# lmfit object for galaxy one and two
p0 = (5500,2.0,0.2,0.3,-2.1,2.1,
      2000,0.3,0.1,-0.2,2.1,2.1)
parameters = lmfit.Parameters()
parameters.add('flux_1', value=p0[0])
parameters.add('hlr_1', value=p0[1], min=0.0)
parameters.add('e1_1', value=p0[2], min=-1.0, max=1.0)
parameters.add('e2_1', value=p0[3], min=-1.0, max=1.0)
parameters.add('x_center1',value=p0[4])
parameters.add('y_center1',value=p0[5])

parameters.add('flux_2', value=p0[6])
parameters.add('hlr_2', value=p0[7], min=0.0)
parameters.add('e1_2', value=p0[8], min=-1.0, max=1.0)
parameters.add('e2_2', value=p0[9], min=-1.0, max=1.0)
parameters.add('x_center2',value=p0[10])
parameters.add('y_center2',value=p0[11])


# Extract params that minimize the difference of the data from the model.
result = lmfit.minimize(resid_2, parameters, args=(im,x_len,y_len,scale))

best_fit = draw_galaxy_2(result.params['flux_1'].value,
                       result.params['hlr_1'].value,
                       result.params['e1_1'].value,
                       result.params['e2_1'].value,
                       result.params['x_center1'].value,
                       result.params['y_center1'].value,
                       result.params['flux_2'].value,
                       result.params['hlr_2'].value,
                       result.params['e1_2'].value,
                       result.params['e2_2'].value,
                       result.params['x_center2'].value,
                       result.params['y_center2'].value,
                       x_len,y_len,scale)

lmfit.report_errors(result.params)

#hsm = im.FindAdaptiveMom()
#print hsm.observed_shape.e1
#print hsm.observed_shape.e2

# Plot results on one figure
fig = plt.figure()
ax1 = fig.add_subplot(131)
a = ax1.imshow(im.array,interpolation='none',origin='lower')
plt.colorbar(a,shrink=0.4)
plt.title('Binned Data')
plt.xlabel('$Pixels$'); plt.ylabel('$Pixels$')
ax2 = fig.add_subplot(132)
b = ax2.imshow(best_fit.array,origin='lower')
plt.colorbar(b,shrink=0.4)
plt.title('Gaussian Surfaces')
plt.xlabel('$Pixels$'); plt.ylabel('$Pixels$')
ax3 = fig.add_subplot(133)
c = ax3.imshow((best_fit-im).array,origin='lower')
plt.colorbar(c,shrink=0.4)
plt.title('Residual')
plt.xlabel('$Pixels$'); plt.ylabel('$Pixels$')
plt.show()