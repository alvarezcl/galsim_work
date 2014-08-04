MSSG edits based on LCA file
Cur: 7/31/2014

Repository contains code used with galsim toolset for the purposes of estimating bias in ellipticity measurements of blended objects.

---------------------------------- Under main:

------------- Utility codes:

-------- drawLibrary.py
Library containing sets of functions for drawing galaxies and residuals
for usage with fitters.

-------- plotutils.py

-------- gauss.py

---------------------------------- Under: mssgversion

------ mssg_galsim_Sersic_test.py --- early code that replaced the Gaussians
in Luis' galsim_first_test with Sersic profiles (cur: 8/2/2014)

------ mssg_2D_Sersic_vary_distance.py --- code that replaced the
Gaussians in Luis' 2D_lmfit_study_same_params_vary_distance.py with
Sersic profiles (cur: 8/3/2014)


---------------------------------- Under: galsim_fitting  (mostly from early July 2014)

------ galsim_First_Test.py
This file looks at the effects of simple galsim objects and methods.

------ 2D_Gaussian_Estimate.py:
This file draws from a bivariate gaussian distribution 
and attempts to find the parameters drawn from using
curve fit. Plots the resultant gaussian contours on the
drawn points and plots 2D histogram.  Lego plots also.

-------- demo_lmfit_1galaxy.py
This file creates the functions necessary for usage with
lmfit and obtains the estimates of HLR,Flux,e1,e2,x0, and
y0 for a 1 object image. 

--------- demo_lmfit_2galaxy.py
This file creates the functions necessary for usage with
lmfit and obtains the estimates of HLR,Flux,e1,e2,x0, and
y0 for a 2 object image.

--------- galsim_lmfit_test.py
This file simply uses lmfit on a galsim object in 
an analytic gaussian form in order to obtain
parameters of the best fit.

--------- galsim_lmfit_loop_test.py
This file loops through different sheared galsim objects
in the form of gaussians and calculates the gaussian
corresponding to the data using lmfit. Ellipticity
estimates are then found and compared to the true value.

-------- lmfit_2dGauss.py
This file draws from a bivariate gaussian distribution 
and attempts to find the parameters drawn from using
lm fit. Plots the resultant gaussian contours on the
drawn points and plots 2D histogram.


---------------------------------- Under: variations

--------- 2D_lmfit_study_same_params_vary_angle.py

Varies angles and sees how the fit results vary

--------- 2D_lmfit_study_same_params_vary_distance.py

Varies sep between the two objs and sees how the fit results vary

--------- vary_dist_obtain_6_params.py

Makes the residual plots for all 6 params, varying the distance
between the two objects

--------- vary_flux_obtain_6_params.py

Makes the residual plots for all 6 params, varying the flux ratio
between the two objects

--------- vary_HLR_obtain_6_params.py

Makes the residual plots for all 6 params, varying the HLR ratio
between the two objects

---------------------------------- Under: animations

--------- 2D_lmfit_study_same_params_vary_distance_animation.py
--------- 2D_lmfit_study_same_params_vary_flux_animation.py
--------- 2D_lmfit_study_same_params_vary_hlr_animation.py

---------------------------------- Under: Images

Output stuff









