# MSSG, based on JEM and LCA code
# 10/15/2014


# Implement the D. Lang et al. LSST AHM slides deblending algorithm
# Don't include sky noise for now, since that seems complicated and
# probably requires defining a segmentation map for each blended
# group.

import numpy as np

def deblend(image, peaks):
    templates = []
    # Step 1: Make symmetric templates
    for peak in peaks:
        templates.append(np.fmin(image, rotate(image, peak)))
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
    rot_bounds = [0,image_width-1,0,image_height-1] # xmin, xmax, ymin, ymax
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
    newimage[ymin:ymax-1, xmin:xmax-1] = image[ymax-1:ymin:-1, xmax-1:xmin:-1]
    return newimage
