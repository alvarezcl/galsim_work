# -*- coding: utf-8 -*-
"""
Created on Wed Jul  9 19:21:46 2014

@author: luis
"""

## This file draws from a bivariate gaussian distribution 
## and attempts to find the parameters drawn from using
## curve fit. Plots the resultant gaussian contours on the
## drawn points and plots 2D histogram.

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import gauss
import plotutils
from astroML.stats.random import bivariate_normal
from astropy.modeling import models, fitting

# Provide seed for the random draw
np.random.seed(0)

# Define the parameters of the 2D gaussian to draw from.
mean = np.array([0, 0])
sigma_1 = 2
sigma_2 = 1
alpha = np.pi / 4

# Draw from the distribution above and return values using astroML library
np.random.seed(0)
xp, cov = bivariate_normal(mean, sigma_1, sigma_2, alpha, size=1000,
                          return_cov=True)

sigma_x = np.sqrt(cov[0, 0])
sigma_y = np.sqrt(cov[1, 1])
sigma_xy = cov[0, 1]

# Plot the scatter of the first draw
#plt.scatter(xp[:,0],xp[:,1])

# Produce a number of points in x-y from another (or same) 2D distribution.
# This is just a test because of the different parametrizations of
# gaussians and how the function below requires a transpose. 
N = 1000
y,x = np.random.multivariate_normal(mean,cov,N).T

# Plot the scatter of the second draw
p2 = plt.scatter(y,x,marker='x',linewidths=0.5) # Note transpose switches y and x

# Prep bins for histogram
bin_size = 0.25
max_edge = 2*(np.sqrt(cov[0][0])+np.sqrt(cov[1][1])) 
min_edge = -max_edge
bin_num = (max_edge-min_edge)/bin_size
bin_numPlus1 = bin_num + 1
bins = np.linspace(min_edge,max_edge,bin_numPlus1)

# Produce 2D histogram
H,xedges,yedges = np.histogram2d(x,y,bins,normed=False)
m = H # Store for astropy
H = H.ravel()
bin_centers_x = (xedges[:-1]+xedges[1:])/2.0
x_len = len(bin_centers_x)
bin_centers_y = (yedges[:-1]+yedges[1:])/2.0
y_len = len(bin_centers_y)
X,Y = np.meshgrid(bin_centers_x,bin_centers_y)

# Astropy attempt at fitting
p_init = models.Gaussian2D(amplitude=H.max(),x_mean=mean[0],y_mean=mean[1],cov_matrix=cov)
fit_p = fitting.SLSQPFitter()
p = fit_p(p_init, X, Y, m)
plt.contour(X,Y,p(X,Y))

# Initial Guess
p0 = (H.max(),mean[0],mean[1],sigma_x,sigma_y,sigma_xy)

# Curve Fit parameters
coeff, var_matrix = curve_fit(gauss.mult_gaussFun_Fit_Ravel,(X,Y),H,p0=p0)

# Check to see if estimated parameters are close to the original.
sigma_x_est = coeff[3]
sigma_y_est = coeff[4]
sigma_xy_est = coeff[5]
sigma_1_est = np.sqrt((sigma_x_est**2+sigma_y_est**2)/2.0 + np.sqrt(((sigma_x_est**2-sigma_y_est**2)/2.0)**2 + sigma_xy**2))
sigma_2_est = np.sqrt((sigma_x_est**2+sigma_y_est**2)/2.0 - np.sqrt(((sigma_x_est**2-sigma_y_est**2)/2.0)**2 + sigma_xy**2))
alpha_est = 0.5*np.arctan(2*sigma_xy/(sigma_x**2-sigma_y**2))

# Evaulate estimated parameters and produce contour.
Z = gauss.mult_gaussFun_Fit((X,Y),*coeff)
Z = Z.reshape((x_len,y_len))
plt.contour(X,Y,Z)
plt.legend([p2],['Points Drawn'])
plt.title('Scattered Points and Gaussian Contour Fit Using Curve_Fit')
plt.xlabel('$x$'); plt.ylabel('$y$')
plt.show()

# Plot 2D histogram
H = H.reshape((x_len,y_len))
fig = plotutils.hist_3dBar_points(x,y,bins)
ax = fig.gca(projection='3d')
ax.set_xlabel('$x$'); ax.set_ylabel('$y$'); ax.set_zlabel('$Frequency$')
ax.text2D(0.05,0.95,'3D Histogram of Point Density',transform=ax.transAxes)
plt.show()