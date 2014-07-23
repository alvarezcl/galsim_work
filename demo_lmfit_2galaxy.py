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
import drawLibrary

# Define independent parameters for image not trying to minimize
x_len = 250
y_len = x_len
scale = 0.2
func = galsim.Gaussian
seed = galsim.BaseDeviate(213456)

# Parameters for the galaxies
flux_a = 5e6; hlr_a = 1; e1_a = 0.0; e2_a = 0.0; x0_a = 0.0; y0_a = 0.0
flux_b = 5e6; hlr_b = 1; e1_b = 0.0; e2_b = 0.0; x0_b = 5.0; y0_b = 0.0


# Draw from some true distribution of points with flux, hlr, e1, e2
# x0, y0, and picture parameters.
im = drawLibrary.drawShoot_galaxy_2(flux_a, hlr_a, e1_a, e2_a, x0_a, y0_a,
                                    flux_b, hlr_b, e1_b, e2_b, x0_b, y0_b,
                                    x_len, y_len, scale, func, func, seed)

# Define some seed that's far from true values and insert into
# lmfit object for galaxy one and two
p0 = (5e6,1.0,0.0,0.0,0,0,
      5e6,1.0,0.0,0.0,2.0,0.0)
parameters = lmfit.Parameters()
parameters.add('flux_1', value=p0[0],min=0.0)
parameters.add('hlr_1', value=p0[1], min=0.0)
parameters.add('e1_1', value=p0[2], min=-1.0, max=1.0)
parameters.add('e2_1', value=p0[3], min=-1.0, max=1.0)
parameters.add('x_center1',value=p0[4])
parameters.add('y_center1',value=p0[5])

parameters.add('flux_2', value=p0[6],min=0.0)
parameters.add('hlr_2', value=p0[7], min=0.0)
parameters.add('e1_2', value=p0[8], min=-1.0, max=1.0)
parameters.add('e2_2', value=p0[9], min=-1.0, max=1.0)
parameters.add('x_center2',value=p0[10])
parameters.add('y_center2',value=p0[11])


# Extract params that minimize the difference of the data from the model.
result = lmfit.minimize(drawLibrary.resid_2, parameters, args=(im,x_len,y_len,scale,func,func))

# Set the differences in an array
flux_ae = result.params['flux_1'].value
hlr_ae = result.params['hlr_1'].value
e1_ae = result.params['e1_1'].value
e2_ae = result.params['e2_1'].value
x0_ae = result.params['x_center1'].value
y0_ae = result.params['y_center1'].value
flux_be = result.params['flux_2'].value
hlr_be = result.params['hlr_2'].value
e1_be = result.params['e1_2'].value
e2_be = result.params['e2_2'].value
x0_be = result.params['x_center2'].value
y0_be = result.params['y_center2'].value

diff = []
diff.append([flux_ae-flux_a,hlr_ae-hlr_a,e1_ae-e1_a,e2_ae-e2_a,x0_ae-x0_a,y0_ae-y0_a,
             flux_be-flux_b,hlr_be-hlr_b,e1_be-e1_b,e2_be-e2_b,x0_be-x0_b,y0_be-y0_b])
diff = np.array(diff)

best_fit = drawLibrary.draw_galaxy_2(flux_ae,hlr_ae,e1_ae,e2_ae,x0_ae,y0_ae,
                                     flux_be,hlr_be,e1_be,e2_be,x0_be,y0_be,
                                     x_len,y_len,scale,func,func)

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