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
import ipdb

############################################## Create the values for the 2D Gaussian
# Provide seed for the random draw
np.random.seed(0)

# Define the parameters of the 2D gaussian to draw from.
mean = np.array([0, 0])
sigma_1 = 2
sigma_2 = 1
alpha = np.pi / 4

# Draw from the distribution above and return values using astroML library
# np.random.seed(0)
xp, cov = bivariate_normal(mean, sigma_1, sigma_2, alpha, size=1000,
                          return_cov=True)

# Luis believes that the cov returned above is analytic, and not from the 1000 points that are being returned

# Get now the separate entries of the covar mat
sigma_xx = np.sqrt(cov[0, 0])
sigma_yy = np.sqrt(cov[1, 1])
sigma_xy = cov[0, 1]

# Plot the scatter of the first draw
#plt.scatter(xp[:,0],xp[:,1])

# Produce a number of points in x-y from another (or same) 2D distribution.
# This is just a test because of the different parametrizations of
# gaussians and how the function below requires a transpose. 
N = 10000
y,x = np.random.multivariate_normal(mean,cov,N).T # Has to be done this way because python throws an error when not transposed (not fully understood)

# Plot the scatter of the second draw
scatterplot = plt.scatter(y,x,marker='.',linewidths=0.1) # Note transpose switches y and x , linewidth is the dot thickness for this plot


################################ Now make the 2D binned histo (lego plot)
# Prep bins for histogram
bin_size = 0.25 # Good binwidth for the lego plot display
max_edge = 2*(np.sqrt(cov[0][0])+np.sqrt(cov[1][1]))  # Get the width from taking 2* geometric mean of sigma_xx and sigma_yy 
min_edge = -max_edge
numbins = (max_edge-min_edge)/bin_size
numbinsPlus1 = numbins + 1
bins = np.linspace(min_edge,max_edge,numbinsPlus1)

# Produce 2D histogram
twoDHist,xedges,yedges = np.histogram2d(x,y,bins,normed=False)
# m = twoDHist # Store for astropy
twoDHist = twoDHist.ravel()
bin_centers_x = (xedges[:-1]+xedges[1:])/2.0
x_len = len(bin_centers_x)
bin_centers_y = (yedges[:-1]+yedges[1:])/2.0
y_len = len(bin_centers_y)
Xmesh,Ymesh = np.meshgrid(bin_centers_x,bin_centers_y)

# Astropy attempt at fitting
#p_init = models.Gaussian2D(amplitude=twoDHist.max(),x_mean=mean[0],y_mean=mean[1],cov_matrix=cov)
#fit_p = fitting.SLSQPFitter()
#p = fit_p(p_init, Xmesh, Ymesh, m)
#plt.contour(Xmesh,Ymesh,p(Xmesh,Ymesh))

# Initial Guess
initguess = (twoDHist.max(),mean[0],mean[1],sigma_xx,sigma_yy,sigma_xy)

# Curve Fit parameters -- get the functional form we're going to fit to from the gauss utility file, feed it the twoDHist data
outputcoeffs, var_matrix = curve_fit(gauss.mult_gaussFun_Fit_Ravel,(Xmesh,Ymesh), twoDHist, p0=initguess)

# Check to see if estimated parameters are close to the original.
sigma_xx_est = outputcoeffs[3]
sigma_yy_est = outputcoeffs[4]
sigma_xy_est = outputcoeffs[5]

# Getting out sigma1 and 2 from the rotated axes
sigma_1_est = np.sqrt((sigma_xx_est**2+sigma_yy_est**2)/2.0 + np.sqrt(((sigma_xx_est**2-sigma_yy_est**2)/2.0)**2 + sigma_xy**2))
sigma_2_est = np.sqrt((sigma_xx_est**2+sigma_yy_est**2)/2.0 - np.sqrt(((sigma_xx_est**2-sigma_yy_est**2)/2.0)**2 + sigma_xy**2))
alpha_est = 0.5*np.arctan(2*sigma_xy/(sigma_xx**2-sigma_yy**2))

print " sigma_1_est, sigma_2_est = " ,  sigma_1_est, sigma_2_est 


# Evaluate estimated parameters and produce contour.
Z = gauss.mult_gaussFun_Fit((Xmesh,Ymesh),*outputcoeffs)
# Z = Z.reshape((x_len,y_len))
plt.contour(Xmesh,Ymesh,Z)
plt.legend([scatterplot],['Points Drawn'])
plt.title('Scattered Points and Gaussian Contour Fit Using Curve_Fit')
plt.xlabel('$x$'); plt.ylabel('$y$')
plt.show()

# Plot 2D histogram (lego plot)
# twoDHist = twoDHist.reshape((x_len,y_len))
fig = plotutils.hist_3dBar_points(x,y,bins)
ax = fig.gca(projection='3d')
ax.set_xlabel('$x$'); ax.set_ylabel('$y$'); ax.set_zlabel('$Frequency$')
ax.text2D(0.05,0.95,'3D Histogram of Point Density',transform=ax.transAxes)
plt.show()

ipdb.set_trace()
