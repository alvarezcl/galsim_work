galsim_work
===========

Repository contains code used with galsim toolset for the purposes of estimating bias in ellipticity measurements of blended objects.

2D_Gaussian_Estimate.py:
This file draws from a bivariate gaussian distribution 
and attempts to find the parameters drawn from using
curve fit. Plots the resultant gaussian contours on the
drawn points and plots 2D histogram.

demo_lmfit_1galaxy.py
This file creates the functions necessary for usage with
lmfit and obtains the estimates of HLR,Flux,e1,e2,x0, and
y0 for a 1 object image. 

demo_lmfit_2galaxy.py
This file creates the functions necessary for usage with
lmfit and obtains the estimates of HLR,Flux,e1,e2,x0, and
y0 for a 2 object image.

drawLibrary.py
Library containing sets of functions for drawing galaxies and residuals
for usage with fitters.

galsim_First_Test.py
This file looks at the effects of simple galsim objects and methods.

galsim_lmfit_test.py
This file simply uses lmfit on a galsim object in 
an analytic gaussian form in order to obtain
parameters of the best fit.

galsim_lmfit_loop_test.py
This file loops through different sheared galsim objects
in the form of gaussians and calculates the gaussian
corresponding to the data using lmfit. Ellipticity
estimates are then found and compared to the true value.

lmfit_2dGauss.py
This file draws from a bivariate gaussian distribution 
and attempts to find the parameters drawn from using
lm fit. Plots the resultant gaussian contours on the
drawn points and plots 2D histogram.

