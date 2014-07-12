## Josh's structure for fitting a 2D gauss onto a galsim generated
## object. I updated for shifted coordinates.

import lmfit
import galsim
import numpy as np
import matplotlib.pyplot as plt

# Essentially the "data" aspect that throws photons at the image for one galaxy
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
def resid(param, target_image, x_len, y_len, scale):
    flux = param['flux'].value
    hlr = param['hlr'].value
    e1 = param['e1'].value
    e2 = param['e2'].value
    x0 = param['x0'].value
    y0 = param['y0'].value
    image = draw_galaxy(flux,hlr,e1,e2,x0,y0,x_len,y_len,scale)
    return (image-target_image).array.ravel()
    
# Define independent parameters for image not trying to minimize
x_len = 100
y_len = 100
scale = 0.2

# Draw from some true distribution of points with flux, hlr, e1, e2
# x0, y0, and picture parameters.
im = drawShoot_galaxy(4000, 1.0, 0.2, 0.3, 1, 1, x_len, y_len, scale)

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
result = lmfit.minimize(resid, parameters, args=(im,x_len,y_len,scale))
best_fit = draw_galaxy(result.params['flux'].value,
                       result.params['hlr'].value,
                       result.params['e1'].value,
                       result.params['e2'].value,
                       result.params['x0'].value,
                       result.params['y0'].value,x_len,y_len,scale)

lmfit.report_errors(result.params)

#hsm = im.FindAdaptiveMom()
#print hsm.observed_shape.e1
#print hsm.observed_shape.e2

fig = plt.figure()
ax1 = fig.add_subplot(131)
ax1.imshow(im.array,origin='lower')
ax2 = fig.add_subplot(132)
ax2.imshow(best_fit.array,origin='lower')
ax3 = fig.add_subplot(133)
blah = ax3.imshow((im-best_fit).array,origin='lower')
plt.colorbar(blah)
plt.show()