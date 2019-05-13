# This script is used to convert feature files 
# into binary numpy files, used by the DN program

import numpy as np
# ID: %#86ccbfb191 
import random
import sys

print "Reading from ",sys.argv[1], ", treating '#' as comments, saving as ", sys.argv[2],  "...";

d = np.genfromtxt(sys.argv[1], comments="#", delimiter=None);
print "Size of read data: ", d.shape

np.save(sys.argv[2], d);
print "Data saved as ", sys.argv[2];
