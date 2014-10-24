# MSSG
# Copied over orig JEM code
# 10/20/2014

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
    image_width, image_height = image.shape
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


if __name__ == '__main__':
    # test rotater
    array = np.zeros((5,5),dtype=np.float)
    array[1,1]=11
    array[1,2]=12
    array[2,2]=22
    array[3,3]=33
    array[4,4]=44
    rot = rotate(array, (1,1))
    print "array"
    print array
    print "rotated about (1,1) wrt center"
    print rot
    print
    rot = rotate(array, (-1,-1))
    print "array"
    print array
    print "rotated about (-1,-1)) wrt center"
    print rot
    print
    rot = rotate(array, (0.5,0.5))
    print "array"
    print array
    print "rotated about (0.5,0.5)) wrt center"
    print rot
    array = array[0:4,0:4]
    print
    rot = rotate(array, (0.5,0.5))
    print "array"
    print array
    print "rotated about (0.5,0.5)) wrt center"
    print rot
    print
    rot = rotate(array, (0.0,0.0))
    print "array"
    print array
    print "rotated about (0,0)) wrt center"
    print rot
