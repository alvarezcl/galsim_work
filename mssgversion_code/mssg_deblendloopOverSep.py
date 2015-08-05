# 10-26-2014
# MSSG

import matplotlib.pyplot as plt
import sys
import os
import numpy as np
import ipdb


e2ain=0.0
e2bin=0.0

e1ain=0.0
e1bin=0.0

############## Run over seps

sepsteps = 10 # How many steps across the pixel to take

for sepstep in xrange(0,sepsteps):
### Run over steps in pixel
    print " ************************************************ We're doing sepstep = ", sepstep/float(sepsteps), "  e1a_in = " , e1ain, "  e2a_in = ", e2ain, " e1b_in = ", e1bin, " e2b_in = ", e2bin
        # Make system call to fitting code
    os.system('python mssg_deblendingTest.py  --e1a ' + str(e1ain) + ' --e2a ' + str(e2ain) + ' --e1b ' + str(e1bin) + ' --e2b ' + str(e2bin) + ' --sepstep ' + str(sepstep/float(sepsteps)) )

#sys.exit()

