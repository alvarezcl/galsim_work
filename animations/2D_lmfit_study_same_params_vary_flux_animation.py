# -*- coding: utf-8 -*-
"""
Created on Thu Jul 24 14:23:16 2014

@author: luis
"""

from __future__ import division
import galsim
import drawLibrary
import matplotlib.pyplot as plt
import lmfit
import numpy as np
import matplotlib.animation as animation
import matplotlib.gridspec as gridspec

# Function definition to return the original data array, best-fit array,
# residual, and correlation matrix with differences and error on e1 and e2.
def run_2_galaxy_vary_distance(flux_a,HLR_a,e1_a,e2_a,x0_a,y0_a,
                               flux_b,HLR_b,e1_b,e2_b,x0_b,y0_b,
                               size_1,size_2,pixel_scale,func_gauss_1,func_gauss_2,dev_1,dev_2):

    im_1 = drawLibrary.drawShoot_galaxy(flux_a,HLR_a,e1_a,e2_a,x0_a,y0_a,
                                        size,size,pixel_scale,func_gauss,dev_1)
    im_2 = drawLibrary.drawShoot_galaxy(flux_b,HLR_b,e1_b,e2_b,x0_b,y0_b,
                                        size,size,pixel_scale,func_gauss,dev_2)
    
    im = im_1 + im_2
    
    # Obtain the image bounds and domain information
    x_cen,y_cen,x,y,X,Y = drawLibrary.return_domain(im)
    
    # -----------------------------------------------------------------------    
    # Estimate the parameters of the image.
    
    # Define some seed that's far from true values and insert into
    # lmfit object for galaxy one and two
    p0 = 1.0*np.array([flux_a,HLR_a,e1_a,e2_a,x0_a,y0_a,
          flux_b,HLR_b,e1_b,e2_b,x0_b,y0_b])
    parameters = lmfit.Parameters()
    parameters.add('flux_1', value=p0[0])
    parameters.add('hlr_1', value=p0[1], min=0.0)
    parameters.add('e1_1', value=p0[2], min=-0.80, max=0.80)
    parameters.add('e2_1', value=p0[3], min=-0.80, max=0.80)
    parameters.add('x_center1',value=p0[4])
    parameters.add('y_center1',value=p0[5])
    
    parameters.add('flux_2', value=p0[6])
    parameters.add('hlr_2', value=p0[7], min=0.0)
    parameters.add('e1_2', value=p0[8], min=-0.80, max=0.80)
    parameters.add('e2_2', value=p0[9], min=-0.80, max=0.80)
    parameters.add('x_center2',value=p0[10])
    parameters.add('y_center2',value=p0[11])
    
    
    # Extract params that minimize the difference of the data from the model.
    result = lmfit.minimize(drawLibrary.resid_2, parameters, args=(im,size_1,size_2,pixel_scale,func_gauss_1,func_gauss_2))
    best_fit = drawLibrary.draw_galaxy_2(result.params['flux_1'].value,
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
                                       size_1,size_2,pixel_scale,func_gauss_1,func_gauss_2)
    
    if result.covar is None:
        print "Undefined Covariance Matrix\n Rerunning"
        return run_2_galaxy_vary_distance(flux_a,HLR_a,e1_a,e2_a,x0_a,y0_a,
                                          flux_b,HLR_b,e1_b,e2_b,x0_b,y0_b,
                                          size_1,size_2,pixel_scale,
                                          func_gauss_1,func_gauss_2,galsim.BaseDeviate(0),galsim.BaseDeviate(0))
    else:
        print "result.covar is defined."
        
    error_diag = np.sqrt(np.diag(result.covar))
    error_mat = np.outer(error_diag,error_diag)
    correlation_mat = result.covar/error_mat

    diff_e1_a = result.params['e1_1'].value - e1_a
    diff_e2_a = result.params['e2_1'].value - e2_a                               
    diff_e1_b = result.params['e1_1'].value - e1_b
    diff_e2_b = result.params['e2_1'].value - e2_b
    
    diff_e = []
    diff_e.append([diff_e1_a,diff_e2_a,diff_e1_b,diff_e2_b])
    error_e = [error_diag[2],error_diag[3],error_diag[8],error_diag[9]]
    
    diff_e = np.array(diff_e)
    error_e = np.array(error_e)
                                                                      
    return im.array,best_fit.array,(best_fit-im).array,correlation_mat,diff_e,error_e

# Begin with parameters for the image----------------------------------------

# Parameters for object a
flux_a = 5e6          # total counts on the image
HLR_a = 1.            # arcsec
e1_a = 0.0
e2_a = 0.0
x0_a = -HLR_a
y0_a = 0.0

# Parameters for object b
flux_b = 5e6       # total counts on the image
HLR_b = HLR_a         # arcsec
e1_b = 0.0
e2_b = 0.0
x0_b = HLR_b
y0_b = 0.0

d = x0_b - x0_a       # arcsec

# Image properties
psf_sigma = 1e-8      # arcsec
pixel_scale = 1/5     # arcsec / pixel
noise = 1e-8          # standard deviation of the counts in each pixel
size = 50            # pixel

# Galsim function definitions
func_gauss = galsim.Gaussian

# Set the RNG
dev_1 = galsim.UniformDeviate(1)
dev_2 = galsim.UniformDeviate(2)

# Obtain instantiation
im,best_fit,residual,correlation_mat,diff_e,error_e = run_2_galaxy_vary_distance(flux_a,HLR_a,e1_a,e2_a,x0_a,y0_a,
                                                                                 flux_b,HLR_b,e1_b,e2_b,x0_b,y0_b,
                                                                                 size,size,pixel_scale,func_gauss,func_gauss,dev_1,dev_2)

# Provide the fontsize
fonts = 10
# Initialize the geometry of the grid
gs = gridspec.GridSpec(3,5)
# Set the figure size
fig = plt.figure(figsize=(20,11))
plt.suptitle(r'$Image\/Scale:$ $\frac{0.2\rq\rq}{Pixel}$',fontsize=fonts+15)
# Add the first subplot
ax1 = fig.add_subplot(gs[0,0])
plt.title('$Image$',fontsize=fonts+8)
plt.ylabel('$Pixels$',fontsize=fonts+8)
plt.tick_params(axis='both',which='major',labelsize=9)
a = ax1.imshow(im,interpolation='none',origin='lower',vmax=im.max(),vmin=im.min())
plt.colorbar(a,shrink=1)
# Add the second subplot
ax2 = fig.add_subplot(gs[1,0])
plt.ylabel('$Pixels$',fontsize=fonts+8)
plt.tick_params(axis='both',which='major',labelsize=9)
plt.title('$Best Fit$',fontsize=fonts+8)
b = ax2.imshow(best_fit,origin='lower',vmax=im.max(),vmin=im.min())
plt.colorbar(b,shrink=1)
# Add the third subplot
ax3 = fig.add_subplot(gs[2,0])
plt.tick_params(axis='both',which='major',labelsize=9)
plt.title('$Residual$',fontsize=fonts+8)
plt.xlabel('$Pixels$',fontsize=fonts+8)
plt.ylabel('$Pixels$',fontsize=fonts+8)
c = ax3.imshow(residual,interpolation='none',origin='lower')
plt.colorbar(c,shrink=1)
# Add the fourth subplot
ax4 = fig.add_subplot(gs[:,1:])
plt.title('Correlation Coefficient Matrix',fontsize=fonts+8)
d = ax4.imshow(correlation_mat,interpolation='none',origin='lower',vmin=-1,vmax=1)
plt.xticks([0,1,2,3,4,5,6,7,8,9,10,11],['$Flux_a$','$HLR_a$','$e1_a$','$e2_a$','$x0_a$','$y0_a$',
           '$Flux_b$','$HLR_b$','$e1_b$','$e2_b$','$x0_b$','$y0_b$'],fontsize=18)
plt.yticks([0,1,2,3,4,5,6,7,8,9,10,11],['$Flux_a$','$HLR_a$','$e1_a$','$e2_a$','$x0_a$','$y0_a$',
           '$Flux_b$','$HLR_b$','$e1_b$','$e2_b$','$x0_b$','$y0_b$'],fontsize=18)
plt.colorbar(d,shrink=1)
text = ax4.text(0,-1.5,'Flux Ratio: %.2f '%(flux_b/flux_a),fontsize=18)


# Update the plots using the following function
def updategal(*args):
    global flux_b,dev_1,dev_2
    flux_r = float(flux_b/flux_a)
    im, best_fit, residual, correlation_mat,diff_e,error_e = run_2_galaxy_vary_distance(flux_a,HLR_a,e1_a,e2_a,x0_a,y0_a,
                                                                                        flux_b,HLR_b,e1_b,e2_b,x0_b,y0_b,
                                                                                        size,size,pixel_scale,func_gauss,func_gauss,dev_1,dev_2)
    print correlation_mat[1][0]
    print correlation_mat[2][0]
    print correlation_mat[3][0]
    print flux_b                                                                                        
    a.set_array(im)
    b.set_array(best_fit)
    c.set_array(residual)
    d.set_array(correlation_mat)
    text.set_text('Flux Ratio: %.2f '%flux_r)
    flux_b += -50000    
    return a,b,c,d
    
time = 0.5    
ani = animation.FuncAnimation(fig, updategal, frames=95, interval=time, blit=True)
ani.save('flux_text.avi', fps=1, codec='avi')
plt.show()