# MSSG
# Copied over current JEM code
# 5.28.2015

# Implement the D. Lang et al. LSST all hands meeting slides deblending
# algorithm.  Don't include sky noise for now, since that seems complicated and
# probably requires defining a segmentation map for each blended group.

import numpy as np

def deblend(image, peaks, interpolate=False, force_interpolate=False):
    """ Quick and dirty deblender.

    Args
    ----
    @param image        A numpy array representing an image of a blend.
    @param peaks        A list of tuples representing the peak positions of objects in the blend.
    @param interpolate  If at least one component of rot_center is not a half-integer, use GalSim
                        to rotate the image.  This currently doesn't work very well!!!
    @param force_interpolate   Use GalSim to rotate the image, even if rot_center components are
                               half-integer and rotation via numpy array operations is possible.
                               This currently doesn't work very well!!!
    
    @returns templates, template_fractions, children
    
    """
    work_image = image+1.e-20

    # Step 1: Make symmetric templates
    templates = [np.fmin(work_image, rotate(work_image, peak,
                                            interpolate=interpolate,
                                            force_interpolate=force_interpolate))
                 for peak in peaks]

    # Step 2: Calculate relative contribution of each template
    template_sum = np.sum(templates, axis=0)
    template_fractions = [template/template_sum * (template_sum != 0) for template in templates]

    # Step 3: Calculate deblended children
    children = [t * image for t in template_fractions]

    return templates, template_fractions, children

def rotate(image, rot_center, interpolate=False, force_interpolate=False):
    """ Rotate an image about a point.  Defaults to using numpy array operations, which means that
    the rotation center is rounded to the nearest half-pixel grid point.  Optionally, use GalSim to
    rotate about arbitrary positions, which requires interpolation since the rotated pixel grid
    doesn't generally align with the original pixel grid.

    @param image        Image to rotate
    @param rot_center   Tuple indicating point to be rotated about.  (0,0) indicates rotate about
                        the geometric center of the image (so at the corner of 4 pixels if the
                        image is even-sized, or at the center of a single pixel if the image is
                        odd-sized image.
    @param interpolate  If at least one component of rot_center is not a half-integer, use GalSim
                        to rotate the image.  This currently doesn't work very well!!!
    @param force_interpolate   Use GalSim to rotate the image, even if rot_center components are
                               half-integer and rotation via numpy array operations is possible.
                               This currently doesn't work very well!!!

    @returns   Rotated image.
    """
    height, width = image.shape

    # Round rot_center to nearest half-integer
    rrot_center = [0.5*np.rint(2*p) for p in rot_center]

    if force_interpolate or (interpolate and rrot_center != rot_center):
        try:
            import galsim
        except:
            raise ImportError("can't interpolate w/o galsim")
        imobj = (galsim.InterpolatedImage(galsim.ImageD(image, scale=1),
                                          calculate_stepk=False,
                                          calculate_maxk=False)
                 .shift(-rot_center[0], -rot_center[1])
                 .rotate(180*galsim.degrees)
                 .shift(rot_center[0], rot_center[1]))
        return imobj.drawImage(nx=width, ny=height, scale=1, method='no_pixel').array

    # image_center is 0-indexed and measured from the lower-left corner of the lower-left pixel.
    image_center = (width * 0.5, height * 0.5)
    rot_pix_center = (image_center[0] + rrot_center[0],
                      image_center[1] + rrot_center[1])

    # compute boundary of rotate region
    rot_width = 2.0*min(rot_pix_center[0], width-rot_pix_center[0])
    rot_height = 2.0*min(rot_pix_center[1], height-rot_pix_center[1])
    rot_bounds = [0,width,0,height] # xmin, xmax, ymin, ymax

    # handle edges falling outside original postage stamp
    if rot_pix_center[0] <= image_center[0]:
        rot_bounds[1] = rot_pix_center[0] + rot_width/2
    else:
        rot_bounds[0] = rot_pix_center[0] - rot_width/2
    if rot_pix_center[1] <= image_center[1]:
        rot_bounds[3] = rot_pix_center[1] + rot_height/2
    else:
        rot_bounds[2] = rot_pix_center[1] - rot_height/2
    xmin, xmax, ymin, ymax = rot_bounds

    # and finally, rotate!
    newimage = np.zeros_like(image)
    newimage[ymin:ymax, xmin:xmax] = (image[ymin:ymax, xmin:xmax])[::-1,::-1]
    return newimage

def test_rotate():
    # test odd-size array
    array = np.array([[0,  0,  0,   0,  0],
                      [0, 11, 12,   0,  0],
                      [0,  0, 22,   0,  0],
                      [0,  0,  0,  33,  0],
                      [0,  0,  0,   0, 44]])

    rot = rotate(array, (1,1)) # rotating about the 33 value pixel
    np.testing.assert_array_almost_equal(rot, np.array([[0,  0,  0,   0,  0],
                                                        [0,  0,  0,   0,  0],
                                                        [0,  0, 44,   0,  0],
                                                        [0,  0,  0,  33,  0],
                                                        [0,  0,  0,   0, 22]]),
                                         5, err_msg="incorrect rotate")

    rot = rotate(array, (-1,-1)) # rotating about the 11 value pixel
    np.testing.assert_array_almost_equal(rot, np.array([[22,  0,  0,   0,  0],
                                                        [12, 11,  0,   0,  0],
                                                        [ 0,  0,  0,   0,  0],
                                                        [ 0,  0,  0,   0,  0],
                                                        [ 0,  0,  0,   0,  0]]),
                                         5, err_msg="incorrect rotate")

    rot = rotate(array, (0.5,0.5)) # rotating about point between 22 and 33
    np.testing.assert_array_almost_equal(rot, np.array([[ 0,  0,  0,   0,  0],
                                                        [ 0, 44,  0,   0,  0],
                                                        [ 0,  0, 33,   0,  0],
                                                        [ 0,  0,  0,  22,  0],
                                                        [ 0,  0,  0,  12, 11]]),
                                         5, err_msg="incorrect rotate")

    # test even-size array
    array = np.array([[0,  0,  0,   0],
                      [0, 11, 12,   0],
                      [0,  0, 22,   0],
                      [0,  0,  0,  33]])

    rot = rotate(array, (0.5,0.5)) # rotating about the 22 value pixel
    np.testing.assert_array_almost_equal(rot, np.array([[ 0,  0,  0,   0],
                                                        [ 0, 33,  0,   0],
                                                        [ 0,  0, 22,   0],
                                                        [ 0,  0, 12,  11]]),
                                         5, err_msg="incorrect rotate")

    rot = rotate(array, (0.0,0.0)) # rotating about point between 11 and 22
    np.testing.assert_array_almost_equal(rot, np.array([[33,  0,  0,   0],
                                                        [ 0, 22,  0,   0],
                                                        [ 0, 12, 11,   0],
                                                        [ 0,  0,  0,   0]]),
                                         5, err_msg="incorrect rotate")

    # test non-square array
    array = np.array([[0,  0,  0,   0],
                      [0, 11, 12,   0],
                      [0,  0, 22,   0],
                      [0,  0,  0,  33],
                      [0,  0,  0,  43]])

    rot = rotate(array, (0.0,0.0)) # rotating about point 1/2 unit left of 22
    np.testing.assert_array_almost_equal(rot, np.array([[43,  0,  0,   0],
                                                        [33,  0,  0,   0],
                                                        [ 0, 22,  0,   0],
                                                        [ 0, 12, 11,   0],
                                                        [ 0,  0,  0,   0]]),
                                         5, err_msg="incorrect rotate")

    rot = rotate(array, (0.5,0.5)) # rotating about point 1/2 unit below 22
    np.testing.assert_array_almost_equal(rot, np.array([[ 0,  0,  0,   0],
                                                        [ 0, 43,  0,   0],
                                                        [ 0, 33,  0,   0],
                                                        [ 0,  0, 22,   0],
                                                        [ 0,  0, 12,  11]]),
                                         5, err_msg="incorrect rotate")

    # test that GalSim rotation agrees with numpy rotation when the rotation axis is a
    # half-integer.
    for center in [(0.0, 0.0), (0.5, 0.5), (0.5, -0.5)]:
        numpy_rot = rotate(array, center)
        galsim_rot = rotate(array, center, force_interpolate=True)
        np.testing.assert_array_almost_equal(
            numpy_rot, galsim_rot, 5,
            err_msg="numpy rotation disagrees with galsim rotation at {}".format(center))

def test_deblend():
    try:
        import galsim
    except:
        print "cant test deblend w/o galsim"

    # check that children of symmetric image show same symmetry
    gal1 = galsim.Gaussian(fwhm=5).shift(-5,0)
    gal2 = galsim.Gaussian(fwhm=5).shift(+5,0)
    gals = gal1 + gal2
    img = gals.drawImage(nx=32, ny=24, scale=1.0)
    
    templates, template_fractions, children = deblend(img.array, [(-5, 0), (5, 0)])
    xflip = children[1][:,::-1]
    symdiff = (children[0] - xflip)/img.array
    np.testing.assert_array_almost_equal(children[0], xflip, 10,
                                         "deblend symmetry failed")

    # check again for non-integer shift
    img = galsim.ImageD(32, 24)
    gal1 = galsim.Gaussian(fwhm=5).shift(-5.2,0)
    gal2 = galsim.Gaussian(fwhm=5).shift(+5.2,0)
    gals = gal1 + gal2
    gals.drawImage(image=img, scale=1)
    templates, template_fractions, children = deblend(img.array, [(-5.2, 0), (5.2, 0)])
    xflip = children[1][:,::-1]
    symdiff = (children[0] - xflip)/img.array
    np.testing.assert_array_almost_equal(children[0], xflip, 10,
                                         "deblend symmetry failed")

    # now check that children of transposed image are simliarly transposed
    # use some noise this time.
    gals.drawImage(image=img, method='phot', n_photons=10000)
    _, _, children = deblend(img.array, [(-3, 0), (3, 0)])
    transimage = img.array.transpose()
    _, _, transchildren = deblend(transimage, [(0, -3), (0, 3)])
    np.testing.assert_array_almost_equal(children[0],
                                         transchildren[0].transpose(),
                                         10,
                                         "transposed child of transposed image not equal to child")
    np.testing.assert_array_almost_equal(children[1],
                                         transchildren[1].transpose(),
                                         10,
                                         "transposed child of transposed image not equal to child")

    # compare array operations rotation to Galsim.rotate
    _, _, children2 = deblend(img.array, [(-3, 0), (3, 0)],
                              interpolate=True, force_interpolate=True)
    np.testing.assert_array_almost_equal(children[0],
                                         children2[0],
                                         9,
                                         "array rotate disagrees with galsim.rotate")
    np.testing.assert_array_almost_equal(children[1],
                                         children2[1],
                                         9,
                                         "array rotate disagrees with galsim.rotate")

if __name__ == '__main__':
    test_rotate()
    test_deblend()
