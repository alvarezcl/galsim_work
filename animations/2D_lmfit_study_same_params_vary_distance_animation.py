# -*- coding: utf-8 -*-
"""
Created on Wed Jul 23 12:29:16 2014

@author: luis
"""

## This file will create an animation of two objects that decrease in
## in separation. 

from __future__ import division
import galsim
import drawLibrary
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
import matplotlib.gridspec as gridspec
import scipy.linalg as scla


# Begin with parameters for the image----------------------------------------

# Parameters for object a
flux_a = 5e6          # total counts on the image
HLR_a = 1.            # arcsec
e1_a = 0.0
e2_a = 0.0

# Parameters for object b
flux_b = flux_a       # total counts on the image
HLR_b = HLR_a         # arcsec
e1_b = 0.0
e2_b = 0.0

# Image properties
psf_sigma = 1e-8      # arcsec
pixel_scale = 1/5     # arcsec / pixel
noise = 1e-8          # standard deviation of the counts in each pixel
size = 100            # pixel

# Set the beginning locations
x0_a = -3
y0_a = 0

x0_b = 3
y0_b = 0.0

distance = x0_b - x0_a # arcsec

# Galsim function definitions
func_gauss = galsim.Gaussian

# Set the RNG
dev_1 = galsim.UniformDeviate(1)
dev_2 = galsim.UniformDeviate(2)

# Noise level
noise_sig = 1e-8

# Obtain instantiation
im,best_fit,result = drawLibrary.run_2_galaxy_full_params(flux_a,HLR_a,e1_a,e2_a,x0_a,y0_a,
                                                        flux_b,HLR_b,e1_b,e2_b,x0_b,y0_b,
                                                        size,size,pixel_scale,func_gauss,func_gauss,dev_1,dev_2)
# Obtain the covariance and correlation matrix.
error_diag = np.sqrt(np.diag(result.covar))
error_mat = np.outer(error_diag,error_diag)
correlation_mat = result.covar/error_mat

# Diagonalize the covariance matrix
EW, EV = scla.eig(result.covar)

# Provide the fontsize
fonts = 10
# Initialize the geometry of the grid
gs = gridspec.GridSpec(3,5)
# Set the figure size
fig = plt.figure(figsize=(20,11))
plt.suptitle(r'$Image\/Scale:$ $\frac{0.2\rq\rq}{Pixel}$',fontsize=fonts+10)
# Add the first subplot
ax1 = fig.add_subplot(gs[0,0])
plt.title('$Image$',fontsize=fonts+8)
plt.ylabel('$Pixels$',fontsize=fonts+8)
plt.tick_params(axis='both',which='major',labelsize=9)
a = ax1.imshow(im.array,interpolation='none',origin='lower',vmax=im.array.max(),vmin=im.array.min())
plt.colorbar(a,shrink=1)
# Add the second subplot
ax2 = fig.add_subplot(gs[1,0])
plt.ylabel('$Pixels$',fontsize=fonts+8)
plt.tick_params(axis='both',which='major',labelsize=9)
plt.title('$Best Fit$',fontsize=fonts+8)
b = ax2.imshow(best_fit.array,origin='lower',vmax=im.array.max(),vmin=im.array.min())
plt.colorbar(b,shrink=1)
# Add the third subplot
ax3 = fig.add_subplot(gs[2,0])
plt.tick_params(axis='both',which='major',labelsize=9)
plt.title('$Residual$',fontsize=fonts+8)
plt.xlabel('$Pixels$',fontsize=fonts+8)
plt.ylabel('$Pixels$',fontsize=fonts+8)
c = ax3.imshow((best_fit-im).array,interpolation='none',origin='lower')
plt.colorbar(c,shrink=1)
# Add the fourth subplot
ax4 = fig.add_subplot(gs[:,1:])
plt.title('Correlation Coefficient Matrix',fontsize=fonts+8)
d = ax4.imshow(correlation_mat,interpolation='none',origin='lower',vmin=-1,vmax=1)
plt.xticks([0,1,2,3,4,5,6,7,8,9,10,11],['$Flux_a$','$HLR_a$','$e1_a$','$e2_a$','$x0_a$','$y0_a$',
           '$Flux_b$','$HLR_b$','$e1_b$','$e2_b$','$x0_b$','$y0_b$'],fontsize=15)
plt.yticks([0,1,2,3,4,5,6,7,8,9,10,11],['$Flux_a$','$HLR_a$','$e1_a$','$e2_a$','$x0_a$','$y0_a$',
           '$Flux_b$','$HLR_b$','$e1_b$','$e2_b$','$x0_b$','$y0_b$'],fontsize=15)
plt.colorbar(d,shrink=1)
text = ax4.text(0,-1.5,'Separation: %.2f arcseconds'%distance,fontsize=15)
plt.show()

# Update the plots using the following function
def updategal(*args):
    # Define globals to update throughout the iteration.
    global x0_a,x0_b
    # Update the distance
    distance = x0_b - x0_a # arcsec
    # If the two objects are on top of each other, separate the two. 
    if distance == 0.0:
        x0_a += 0.5; x0_b += -0.5
        distance = x0_b - x0_a # arcsec
    print "separation = %.2f Arcsec"%distance
    # Obtain the new information from the resulting parameter change.
    im, best_fit,result = drawLibrary.run_2_galaxy_full_params(flux_a,HLR_a,e1_a,e2_a,x0_a,y0_a,
                                                               flux_b,HLR_b,e1_b,e2_b,x0_b,y0_b,
                                                               size,size,pixel_scale,func_gauss,func_gauss,dev_1,dev_2,noise_sig)                                                                                    
    # Check that the covariance matrix of the parameters is empty.  
    while result.covar is None:
        print "Undefined Covariance Matrix\n Rerunning"
        im, best_fit, result = drawLibrary.run_2_galaxy_full_params(flux_a,HLR_a,e1_a,e2_a,x0_a,y0_a,
                                                                  flux_b,HLR_b,e1_b,e2_b,x0_b,y0_b,
                                                                  size,size,pixel_scale,func_gauss,func_gauss,dev_1,dev_2,noise_sig)
    else:
        print "result.covar is defined."    
    
    # Retrive the covariance and correlation matrix.
    error_diag = np.sqrt(np.diag(result.covar))
    error_mat = np.outer(error_diag,error_diag)
    correlation_mat = result.covar/error_mat
    # Set the new data array and text.
    a.set_array(im.array)
    b.set_array(best_fit.array)
    c.set_array((best_fit-im).array)
    d.set_array(correlation_mat)
    text.set_text('Separation is %.2f arcseconds'%distance)
    # Update the parameters. 
    x0_a += 0.5; x0_b += -0.5
    
    return a,b,c,d

## Save the animation.    
#time = 0.2
#frame_num = 8
#fps = 1
#anim = animation.FuncAnimation(fig, updategal, frames=frame_num, interval=time, blit=True)
##anim.save('dist_1e3_noise.avi', codec='avi', fps=fps)
#plt.show()