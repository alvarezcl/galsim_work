# MSSG
# Copied over orig JEM code
# 11.18.2014

# Implement the D. Lang et al. LSST AHM slides deblending algorithm
# Don't include sky noise for now, since that seems complicated and
# probably requires defining a segmentation map for each blended
# group.

import numpy as np

def deblend(image, peaks):
    work_image = image+1.e-20
    templates = []
    # Step 1: Make symmetric templates
    for peak in peaks:
        templates.append(np.fmin(work_image, rotate(work_image, peak)))
    # Step 2: Calculate relative contribution of each template
    template_sum = reduce(lambda x,y: x+y, templates, 0)
    template_fractions = []
    for template in templates:
        template_frac = template/template_sum
        template_frac[np.isnan(template_frac)]=0.0
        template_fractions.append(template_frac)
    # Step 3: Calculate deblended children
    children = []
    for template_fraction in template_fractions:
        children.append(template_fraction * image)
    return templates, template_fractions, children

def rotate(image, peak):
    # Assume that origin is in the geometric center of the image (so at the corner of 4 pixels if
    # even-sized image, or at the center of a single pixel if odd-sized image).
    image_height, image_width = image.shape
    # This center is 0-indexed and measured from the lower-left corner of the lower-left pixel.
    image_center = (image_width * 0.5, image_height * 0.5)
    rot_pix_center = (image_center[0] + peak[0],
                      image_center[1] + peak[1])
    rot_width = 2.0*min(rot_pix_center[0], image_width-rot_pix_center[0])
    rot_height = 2.0*min(rot_pix_center[1], image_height-rot_pix_center[1])
    rot_bounds = [0,image_width,0,image_height] # xmin, xmax, ymin, ymax
    if rot_pix_center[0] <= image_center[0]:
        rot_bounds[1] = rot_pix_center[0] + rot_width/2
    else:
        rot_bounds[0] = rot_pix_center[0] - rot_width/2
    if rot_pix_center[1] <= image_center[1]:
        rot_bounds[3] = rot_pix_center[1] + rot_height/2
    else:
        rot_bounds[2] = rot_pix_center[1] - rot_height/2
    xmin, xmax, ymin, ymax = rot_bounds
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

    rot = rotate(array, (0.5,0.5)) # rotating about point 1/2 unint below 22
    np.testing.assert_array_almost_equal(rot, np.array([[ 0,  0,  0,   0],
                                                        [ 0, 43,  0,   0],
                                                        [ 0, 33,  0,   0],
                                                        [ 0,  0, 22,   0],
                                                        [ 0,  0, 12,  11]]),
                                         5, err_msg="incorrect rotate")

def test_deblend():
    try:
        import galsim
    except:
        print "cant test deblend w/o galsim"

    # check that children of symmetric image show same symmetry
    img = galsim.ImageD(32, 24)
    gal1 = galsim.Gaussian(fwhm=5).shift(-5,0)
    gal2 = galsim.Gaussian(fwhm=5).shift(+5,0)
    gals = gal1 + gal2
    gals.drawImage(image=img)
    templates, template_fractions, children = deblend(img.array, [(-3, 0), (3, 0)])
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

if __name__ == '__main__':
    test_rotate()
    test_deblend()
