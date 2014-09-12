# Orig code by JEM (Josh E Meyers)
# Current by MSSG
# 9-11-2014

import numpy as np
import galsim
import lmfit

sky_level = 1.e6
pixel_scale = 0.2
sky_level_pixel = sky_level * pixel_scale**2  # This is the num of photons per pixel?

random_seed = np.random.randint(2**62)  # Why give it such a large seed..?

# Create the orig image
def create_true_image(params):
    e1=params['e1'].value
    e2=params['e2'].value
    hlr=params['hlr'].value
    x=params['x'].value
    y=params['y'].value
    flux=params['flux'].value

    image = galsim.ImageF(32, 32, scale=pixel_scale)
    gal = galsim.Gaussian(half_light_radius=hlr).shear(e1=e1, e2=e2) * flux
    psf = galsim.Moffat(beta=5, fwhm=0.7).shear(e1=0.01, e2=-0.03)
    final = galsim.Convolve(gal,psf)
    final.drawImage(image=image)
    return image

# Create an image with noise in it
def create_noisy_image(image, gal_signal_to_noise, udRandomSeed):
    image2 = galsim.ImageF(32, 32, scale=pixel_scale)*0.0 + image
    noise = galsim.PoissonNoise(udRandomSeed, sky_level=sky_level_pixel)
    image2.addNoiseSNR(noise, gal_signal_to_noise)
    return image2

# Get the resid of the true - noise img
def image_resid(params, noisy_image):
    model = create_true_image(params)
    return (model-noisy_image).array.ravel()

# Main body
if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()

    parser.add_argument("--niter", default=10  , type=int, help="number of iterations")
    parser.add_argument("--outfile", default="noise_bias.txt", help="output text filename")
    parser.add_argument("--plot", action="store_true", help="Histogram output values and make plot")
    parser.add_argument("--snr", default=20, type=int, help="Signal-to-noise ratio")
    args = parser.parse_args()

    e1val = 0.3
    e2val = 0.0

    params = lmfit.Parameters()
    params.add('e1', value=e1val)
    params.add('e2', value=e2val)
    params.add('flux', value=1.0)
    params.add('x', value=0.2)
    params.add('y', value=0.1)
    params.add('hlr', value=0.43)

    true_image = create_true_image(params)

    snr = args.snr
    e1vec = []
    e2vec = []
    niter = args.niter
    for i in range(niter):
        print " On iteration ", i
        udRandomSeed = galsim.UniformDeviate(random_seed+i)
        try:
            fit_params = lmfit.Parameters()
            fit_params.add('e1', value=0.3)
            fit_params.add('e2', value=0.0)
            fit_params.add('flux', value=1.)
            fit_params.add('x', value=0.2)
            fit_params.add('y', value=0.1)
            fit_params.add('hlr', value=0.43)
            noisy_image = create_noisy_image(true_image, snr, udRandomSeed)
            fit_params['flux'].value = np.sum(noisy_image.array)
            ml = lmfit.minimize(image_resid, fit_params, args=(noisy_image,))
            e1vec.append(ml.params['e1'].value)
            e2vec.append(ml.params['e2'].value)
        except:
            continue

    with open(args.outfile, 'w') as f:
        for e1 in e1vec:
            f.write("{}\n".format(e1))

    print " np.mean(e1vec) = " , np.mean(e1vec)
    print  " np.std(e1vec) / np.sqrt(len(e1vec)) ", np.std(e1vec) / np.sqrt(len(e1vec))

    if args.plot:
        import matplotlib.pyplot as plt
        fig = plt.figure()

        fig.suptitle("SNR = " + str(snr))

        binsfac = 3

        ax = fig.add_subplot(121)
        ax.hist(np.array(e1vec), bins=binsfac * np.sqrt(len(e1vec))/2)
        ax.set_xlabel(r"$e_1$")
        ax.set_ylabel("#")
        ax.set_title("e1 = "+str(e1val))

        ax = fig.add_subplot(122)
        ax.hist(np.array(e2vec), bins=binsfac * np.sqrt(len(e2vec))/2)
        ax.set_xlabel(r"$e_2$")
        ax.set_ylabel("#")
        ax.set_title("e2 = "+str(e2val))

        plt.show()
