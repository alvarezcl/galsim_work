
## This file creates the functions necessary for usage with
## lmfit and obtains the estimates of HLR,Flux,e1,e2,x0, and
## y0 for a 1 object image. 

import lmfit
import galsim
import numpy as np
import matplotlib.pyplot as plt

# The "data" aspect that throws photons at the image for one galaxy
def drawShoot_galaxy(flux, hlr, e1, e2, x0, y0, x_len, y_len, scale):
    gal = galsim.Gaussian(half_light_radius=hlr, flux=flux)
    gal = gal.shear(e1=e1, e2=e2)
    gal = gal.shift(x0,y0)
    image = galsim.ImageD(x_len, y_len, scale=scale)
    image = gal.drawShoot(image=image)
    return image

# Continous version, or what you would call the "model" for one galaxy
def draw_galaxy(flux, hlr, e1, e2, x0, y0, x_len, y_len, scale):
    gal = galsim.Gaussian(half_light_radius=hlr, flux=flux)
    gal = gal.shear(e1=e1, e2=e2)
    gal = gal.shift(x0,y0)
    image = galsim.ImageD(x_len, y_len, scale=scale)
    image = gal.draw(image=image)
    return image

# Take the difference of the data and the model for one galaxy
def resid_1(param, target_image, x_len, y_len, scale):
    flux = param['flux'].value
    hlr = param['hlr'].value
    e1 = param['e1'].value
    e2 = param['e2'].value
    x0 = param['x0'].value
    y0 = param['y0'].value
    image = draw_galaxy(flux,hlr,e1,e2,x0,y0,x_len,y_len,scale)
    return (image-target_image).array.ravel()
    
# Define independent parameters for image not trying to minimize
x_len = 40
y_len = 40
scale = 0.2

# Draw from some true distribution of points with flux, hlr, e1, e2
# x0, y0, and picture parameters.
im = drawShoot_galaxy(4000, 1.0, 0.2, 0.3, -0.2, 0, x_len, y_len, scale)

# Define some seed that's far from true values and insert into
# lmfit object.
p0 = (3000,1.2,0.0,0.0,2,2)
parameters = lmfit.Parameters()
parameters.add('flux', value=p0[0])
parameters.add('hlr', value=p0[1], min=0.0)
parameters.add('e1', value=p0[2], min=-1.0, max=1.0)
parameters.add('e2', value=p0[3], min=-1.0, max=1.0)
parameters.add('x0',value=p0[4])
parameters.add('y0',value=p0[5])

# Extract params that minimize the difference of the data from the model.
result = lmfit.minimize(resid_1, parameters, args=(im,x_len,y_len,scale))
best_fit = draw_galaxy(result.params['flux'].value,
                       result.params['hlr'].value,
                       result.params['e1'].value,
                       result.params['e2'].value,
                       result.params['x0'].value,
                       result.params['y0'].value,x_len,y_len,scale)

lmfit.report_errors(result.params)

# Use galsim methods to find ellipticities
#hsm = im.FindAdaptiveMom()
#print hsm.observed_shape.e1
#print hsm.observed_shape.e2

# Plot on one figure
fig = plt.figure()
ax1 = fig.add_subplot(131)
a = ax1.imshow(im.array,interpolation='none',origin='lower')
plt.colorbar(a,shrink=0.4)
plt.title('Binned Data')
plt.xlabel('$Pixels$'); plt.ylabel('$Pixels$')
ax2 = fig.add_subplot(132)
b = ax2.imshow(best_fit.array,origin='lower')
plt.colorbar(b,shrink=0.4)
plt.title('Gaussian Surface')
plt.xlabel('$Pixels$'); plt.ylabel('$Pixels$')
ax3 = fig.add_subplot(133)
c = ax3.imshow((best_fit-im).array,origin='lower')
plt.colorbar(c,shrink=0.4)
plt.title('Residual')
plt.xlabel('$Pixels$'); plt.ylabel('$Pixels$')
plt.show()