MSSGversion_code dir
Start: 7/31/2014
Cur: 10/6/2014

Repository contains code used with galsim toolset for the purposes of
estimating bias in ellipticity measurements of blended objects.


----------------------------------------------------------------------- MSSG code

---------------------------------- Under: mssgversion

------ mssg_galsim_Sersic_test.py --- early code that replaced the
Gaussians in Luis' galsim_first_test with Sersic profiles -- creates a
gal with 2 components, adds them, and writes out a FITS file too (cur:
8/2/2014 ) (Didn't trace through in great detail -- 10/9/2014)

------ mssg_2D_Gaussian_Estimate.py:

This file draws from a bivariate gaussian distribution and attempts to
find the parameters drawn from using curve fit. Plots the resultant
gaussian contours on the drawn points and plots 2D histogram.  Lego
plots also.  (Didn't trace through in great detail -- 10/9/2014)



-------- mssg_Sersic_Gaussian_comparison.py (8/2014)

Calls mssg_noiseLibrary to create ind gals of Sersic and Gaussian type w/ or w/o
psf convln, then compares slices and proj's of the two to one another
(Didn't trace through in great detail -- 10/9/2014)


------ mssg_single_Sersic_fit.py  (8/2014)

Simplified code to fit a single comp Sersic profile (Didn't trace
through in great detail -- 10/9/2014) 

------ mssg_double_Sersic_fit.py

Code building up from the single Sersic one to fit 2 Sersic's, or a
Sersic + Gaussian. (Didn't trace through in great detail -- 10/9/2014)



------ mssg_2D_Sersic_vary_distance.py 

Code that replaced the Gaussians in Luis'
2D_lmfit_study_same_params_vary_distance.py with Sersic profiles (cur:
8/3/2014)

This is throwing an error currently -- 10/9/2014.


------- mssg_Sersic_singleComponent_2object_fit.py

Nice and extensive code that varies distance btwn 2 objects and shows
the corr mat, then finally plots the various properties vs. one
another. (Didn't trace through in great detail -- 10/9/2014)



   -----------------------------------------------

-------- mssg_runOverSNRRange.py (9/2014)

Basically shell code to loop over ellips and SNR and call
mssg_singleObjectNoiseStudy.py which will run the fit over multiple
noise rzns (10/2014)

-------- mssg_singleObjectNoiseStudy.py

Main body: 

 - set  e1 and e2,   num_trials , SNR to use 
 - create the 'true' image with the input params
 - Then loop through the num of trials and create a noise image using the above, with dif noise rzns
    - And fit this using the image resid, and the params
    - Extract results of the fit and append resp to 4 vecs: e1, e2, HLR, flux
 - Write out the results of the e1 and e2 vecs to a file that has input e1 and e2 and SNR in its name
 - If plot flag is set, show the results



 (10/2014)

-------- mssg_makeNoiseBiasStudyPlots.py

-------- mssg_plotSlopeVsSNR.py

--------------------------------------- Deblending  (11/2014)

 mssg_deblendloop.py calls
    mssg_deblendingTest.py calls
       mssg_deblend.py



mssg_deblendingTest.py

-  mssg_makeDeblendingStudyPlots.py -- will plot the files on disk

mssg_makeDeblendingHistos.py

mssg_makeRotMatchedDeblendingPlots.py


-------------------------------------------- Mar 2015

- mssg_simpleSingleGalaxyFit.py  -- stripped down code to just clarify for myself how to do fits

- mssg_simultaneousFitAndDeblend.py -- to feed a simfit centroid values into deblend code

- mssg_simpleSimultaneousFitAndDeblend.py -- adding in stuff steadily to make this stripped down code work to fit then deblend

----------------------------------------------------------------------------------
----------------------------------------------------------------------- Luis code
----------------------------------------------------------------------------------
---------------------------------- Under main:

------------- Utility codes:

-------- drawLibrary.py

Library containing sets of functions for drawing galaxies and
residuals for usage with fitters.

-------- plotutils.py

-------- noiseLibrary.py

-------- gauss.py

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









