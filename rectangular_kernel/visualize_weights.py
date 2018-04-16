#!/usr/bin/env python

import numpy as np
import os, argparse, math, sys
import matplotlib.pyplot as plt
from scipy import signal
from scipy.fftpack import fft, fftshift
import numpy.lib.recfunctions as rfn


parser = argparse.ArgumentParser(description="Show weights to determine the last useful lag")

parser.add_argument('-w', '--weightfile', type=str, nargs=1, help="File containing the weights of the ACF ([Lags, Weights])")

args = parser.parse_args()

weightfile = args.weightfile[0]

lags = np.genfromtxt(weightfile)[:,0]
weights = np.genfromtxt(weightfile)[:,1]

plt.plot(lags,weights,'k-')
plt.xlabel("Lags")
plt.ylabel("weights")
plt.show()
