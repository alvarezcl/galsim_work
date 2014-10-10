# Orig code by JEM (Josh E Meyers)
# Current by MSSG
# 9-11-2014

import numpy as np
import galsim
import lmfit
import matplotlib.pyplot as plt
import ipdb

# Define some params
sky_level = 1.e6
pixel_scale = 0.2
sky_level_pixel = sky_level * pixel_scale**2  # This is the num of photons per pixel?

random_seed = np.random.randint(2**62)  # Why give it such a large seed..?


# Create the orig true image w/o noise
def create_true_image(params):
    e1=params['e1'].value
    e2=params['e2'].value
    hlr=params['hlr'].value
    x=params['x'].value
    y=params['y'].value
    flux=params['flux'].value

    image = galsim.ImageF(32, 32, scale=pixel_scale)
    gal = galsim.Gaussian(half_light_radius=hlr).shear(e1=e1, e2=e2) * flux
    psf = galsim.Moffat(beta=3, fwhm=2.85).shear(e1=0.0, e2=-0.0)
    final = galsim.Convolve(gal,psf)
    final.drawImage(image=image)
    return image

# Create an image with noise in it
def create_noisy_image(image, gal_signal_to_noise, udRandomSeed):
    image2 = galsim.ImageF(32, 32, scale=pixel_scale)*0.0 + image          # Why do we need to have the first part if it's 0, why not just assign image2 = image ?
    noise = galsim.PoissonNoise(udRandomSeed, sky_level=sky_level_pixel)
    image2.addNoiseSNR(noise, gal_signal_to_noise)
    return image2

# Get the resid of the true - noise img
def image_resid(params, noisy_image):
    model = create_true_image(params)
    return (model-noisy_image).array.ravel()


################################# Main body
if __name__ == "__main__":

# Parse command line args
    from argparse import ArgumentParser
    parser = ArgumentParser()

    parser.add_argument("--ntrials", default=5  , type=int, help="number of iterations")
    parser.add_argument("--outfile", default="noisebiasoutput/noise_bias", help="output text filename")
    parser.add_argument("--plot", action="store_true", help="Histogram output values and make plot") # If turned on, this will plot stuff
    parser.add_argument("--snr", default=20, type=int, help="Signal-to-noise ratio")
    parser.add_argument("--e1", default=0, type=float, help="e1 in")
    parser.add_argument("--e2", default=0, type=float, help="e2 in")

    args = parser.parse_args()

# The e1 and e2, SNR to use
    e1val = args.e1
    e2val = args.e2
    snr = args.snr
# Number of trials
    ntrials = args.ntrials # How many iterations to run
    num_trials = ntrials

# The initial params to use in creating the true image 
    params = lmfit.Parameters()
    params.add('e1', value=e1val)
    params.add('e2', value=e2val)
    params.add('flux', value=1.0)
    params.add('x', value=0.2)
    params.add('y', value=0.1)
    params.add('hlr', value=0.43)

# Make the true image first w/o noise
    true_image = create_true_image(params)

# Declare vecs we'll write to disk later
    e1vec = []
    e2vec = []

    HLRvec = []
    fluxvec = []

################################### Now iterate over the num of trials
    for i in range(num_trials):
        print " On iteration ", i
        udRandomSeed = galsim.UniformDeviate(random_seed+i) # Pick the seed
        try:  # The initial guesses
            fit_params = lmfit.Parameters()
            fit_params.add('e1', value=0.3)
            fit_params.add('e2', value=0.0)
            fit_params.add('flux', value=1.)
            fit_params.add('x', value=0.2)
            fit_params.add('y', value=0.1)
            fit_params.add('hlr', value=0.43)
            noisy_image = create_noisy_image(true_image, snr, udRandomSeed) # Send this function the true img, and it will add noise
            fit_params['flux'].value = np.sum(noisy_image.array)            # The true sum is the sum of values of all pixels
            print " Actual flux ", flux

            ml = lmfit.minimize(image_resid, fit_params, args=(noisy_image,)) # Do the fit
            e1vec.append(ml.params['e1'].value)  # Get out e1 val from fit
            e2vec.append(ml.params['e2'].value)  # Get out e2 val

            HLRvec.append(ml.params['hlr'].value)
            fluxvec.append(ml.params['flux'].value)
                
        except:
            continue
            
######################################### End of noise trials

    with open(args.outfile+"_SNRof_"+str(snr)+"_e1In_"+str(e1val)+"_e2In_"+str(e2val)+"e1out.txt", 'w') as f:
        for e1 in e1vec:
            f.write("{}\n".format(e1))
                    
    with open(args.outfile+"_SNRof_"+str(snr)+"_e1In_"+str(e1val)+"_e2In_"+str(e2val)+"e2out.txt", 'w') as f:
        for e2 in e2vec:
            f.write("{}\n".format(e2))
            
# Just to see on screen
    print " np.mean(e1vec) = " , np.mean(e1vec)
    print  " np.std(e1vec) / np.sqrt(len(e1vec)) ", np.std(e1vec) / np.sqrt(len(e1vec))

    print " np.mean(e2vec) = " , np.mean(e2vec)
    print  " np.std(e2vec) / np.sqrt(len(e2vec)) ", np.std(e2vec) / np.sqrt(len(e2vec))

# If plotting is chosen
    if args.plot:
        fig = plt.figure()
    
        fig.suptitle("SNR = " + str(snr))
        
        binsfac = 3
        
        ax = fig.add_subplot(221)
        ax.hist(np.array(e1vec), bins=binsfac * np.sqrt(len(e1vec))/2)
        ax.set_xlabel(r"$e_1$ ")
        ax.set_ylabel("#")
        ax.set_title("e1_in = "+str(e1val))
        
        ax = fig.add_subplot(222)
        ax.hist(np.array(e2vec), bins=binsfac * np.sqrt(len(e2vec))/2)
        ax.set_xlabel(r"$e_2$ ")
        ax.set_ylabel("#")
        ax.set_title("e2_in = "+str(e2val))

        ax = fig.add_subplot(212)
        plt.scatter(100*np.ones(num_trials),resid_vector1)
        
        plt.show()
