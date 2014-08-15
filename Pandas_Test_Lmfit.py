# -*- coding: utf-8 -*-
"""
Created on Wed Jul 30 15:16:18 2014

@author: luis
"""

## This is a concurrent Pandas Test file with the goal of implementing an
## lmfit program with reparameterization of two object galaxies.

from __future__ import division
from pandas import Series, DataFrame
import pandas as pd
import numpy as np
import drawLibrary
import galsim
import lmfit
import matplotlib.pyplot as plt
import scipy.linalg as scla
import matplotlib.gridspec as gridspec

# Parameters for object a
flux_a = 5e6          # total counts on the image
hlr_a = 1         # arcsec
e1_a = 0.0
e2_a = 0.0
q_a = (1+(e1_a**2 + e2_a**2))/(1-(e1_a**2 + e2_a**2))

# Parameters for object b
flux_b = 5e6       # total counts on the image
hlr_b = hlr_a         # arcsec
e1_b = 0.0
e2_b = 0.0
q_b = (1+(e1_b**2 + e2_b**2))/(1-(e1_b**2 + e2_b**2))

sum_flux = flux_a + flux_b
dif_flux = flux_a - flux_b

# Image properties
psf_sigma = 1e-8      # arcsec
pixel_scale = 1/5     # arcsec / pixel
noise = 1e-8          # standard deviation of the counts in each pixel
size = 100            # pixel

# Set the centroid location of a
x0_a = -3
y0_a = 0
centroid_a = np.array([x0_a,y0_a])

# Get the angle of the ellipse for a
if e1_a == 0 and e2_a == 0:
    phi_a = 0
elif e1_a == 0 and e2_a != 0:    
    phi_a = np.pi/2
else:    
    phi_a = np.arctan2(e2_a,e1_a)

# Set the centroid location of b
x0_b = 3
y0_b = 0
centroid_b = np.array([x0_b,y0_b])

distance = x0_b - x0_a

# Get the angle of the ellipse for b
if e1_b == 0 and e2_b == 0:
    phi_b = 0
elif e1_b == 0 and e2_b != 0:    
    phi_b = np.pi/2
else:    
    phi_b = np.arctan2(e2_b,e1_b)

# Calculate values for use in correlations
x_sum = (flux_a*x0_a + flux_b*x0_b)/(flux_a + flux_b)
x_dif = (flux_a*x0_a - flux_b*x0_b)/(flux_a + flux_b)

y_sum = (flux_a*y0_a + flux_b*y0_b)/(flux_a + flux_b)
y_dif = (flux_a*y0_a - flux_b*y0_b)/(flux_a + flux_b)

# Calculate the second central moment of object a and b and respective sigmas

Q_a = drawLibrary.Q_Gauss(hlr_a,e1_a,e2_a)
sigma_minus_a = np.linalg.det(Q_a)**(1/4)
sigma_plus_a = ((Q_a[0,0] + Q_a[1,1])/2)**(1/2)
q_a_from_Q_a = (sigma_plus_a/sigma_minus_a)**2

Q_b = drawLibrary.Q_Gauss(hlr_b,e1_b,e2_b)
sigma_minus_b = np.linalg.det(Q_b)**(1/4)
sigma_plus_b = ((Q_b[0,0] + Q_b[1,1])/2)**(1/2)
q_b_from_Q_b = (sigma_plus_b/sigma_minus_b)**2

# Calculate the second central moment of the mixture
Q = drawLibrary.Q_two_mixture(flux_a,hlr_a,e1_a,e2_a,x0_a,y0_a,
                              flux_b,hlr_b,e1_b,e2_b,x0_b,y0_b)
sigma_minus = np.linalg.det(Q)**(1/4)
sigma_plus = ((Q[0,0] + Q[1,1])/2)**(1/2)

# Calculate the angles again and see if they are the same.
if Q_a[0,1] == 0:
    phi_a_from_Q = 0
else:
    phi_a_from_Q = 0.5 * np.arctan2(2.0 * Q_a[0,1],Q_a[0,0] - Q_a[1,1])
    
if Q_b[0,1] == 0:
    phi_b_from_Q = 0
else:
    phi_b_from_Q = 0.5 * np.arctan2(2.0 * Q_a[0,1],Q_a[0,0] - Q_a[1,1])

# Galsim function definitions
func_gauss = galsim.Gaussian

# Set the RNG
dev_1 = galsim.BaseDeviate(1)
dev_2 = galsim.BaseDeviate(2)
dev_3 = galsim.BaseDeviate(3)


# Noise level
sky_level = 0

# Obtain instantiation
im,best_fit,result = drawLibrary.run_2_galaxy(sum_flux,dif_flux,x_sum,x_dif,y_sum,y_dif,
                                              hlr_a,e1_a,e2_a,
                                              hlr_b,e1_b,e2_b,
                                              size,size,pixel_scale,func_gauss,func_gauss,dev_1,dev_2,dev_3)

# Obtain the covariance and correlation matrix.
error_diag = np.sqrt(np.diag(result.covar))
error_mat = np.outer(error_diag,error_diag)
correlation_mat = result.covar/error_mat

# Diagonalize the covariance matrix
EW, EV = scla.eig(result.covar)

# Report errors
lmfit.report_errors(result.params)

# Extract the values from the result
flux_a_est,hlr_a_est,e1_a_est,e2_a_est,x0_a_est,y0_a_est,flux_b_est,hlr_b_est,e1_b_est,e2_b_est,x0_b_est,y0_b_est = drawLibrary.extract_physical_params(result)


# Plots

# Get the domain of the image in pixels
x_cen,y_cen,x,y,X,Y = drawLibrary.return_domain(im)
# Convert all centroids to pixels and with respect to the origin in the image
centroid_a_image = centroid_a/pixel_scale + np.array([x_cen,y_cen]) - np.array([1.5,1.5])
centroid_b_image = centroid_b/pixel_scale + np.array([x_cen,y_cen]) - np.array([1.5,1.5])

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
c = ax3.imshow((best_fit - im).array,interpolation='none',origin='lower')
plt.colorbar(c,shrink=1)
# Add the fourth subplot
ax4 = fig.add_subplot(gs[:,1:])
plt.title('Correlation Coefficient Matrix',fontsize=fonts+8)
d = ax4.imshow(correlation_mat,interpolation='none',vmin=-1,vmax=1)
plt.xticks([0,1,2,3,4,5,6,7,8,9,10,11],['$F_{sum}$','$F_{dif}$','$x_{sum}$','$x_{dif}$','$y_{sum}$','$y_{dif}$','$HLR_a$','$e1_a$','$e2_a$',
           '$HLR_b$','$e1_b$','$e2_b$'],fontsize=15)
plt.yticks([0,1,2,3,4,5,6,7,8,9,10,11],['$F_{sum}$','$F_{dif}$','$x_{sum}$','$x_{dif}$','$y_{sum}$','$y_{dif}$','$HLR_a$','$e1_a$','$e2_a$',
           '$HLR_b$','$e1_b$','$e2_b$'],fontsize=15)
plt.colorbar(d,shrink=1)
for (i,j), val in np.ndenumerate(correlation_mat):
    ax4.annotate('%0.2f'%(val), (j,i), ha='center', va='center',size=12)
text = ax4.text(0,-1.5,'Separation: %.2f arcseconds'%distance,fontsize=15)
plt.show()

