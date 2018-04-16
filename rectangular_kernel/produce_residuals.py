#!/usr/bin/env python

import numpy as np
import os, argparse, math, sys
import matplotlib.pyplot as plt
import numpy.lib.recfunctions as rfn

parser = argparse.ArgumentParser(description="Compute ACF")
parser.add_argument('-f', '--dmfile', type=str, nargs=1,
                    help="File containing DM variations, as observation, date, MJD, DM, edm, angle")
parser.add_argument('-m', '--modelfile', type=str, nargs=1,
                    help="File containing the model, as MJD,DM")

args = parser.parse_args()
dmfile = args.dmfile[0]
mfile = args.modelfile[0]

data = np.genfromtxt(dmfile, dtype= None, names="observation, date, MJD, DM, edm, angle")
model = np.genfromtxt(mfile, dtype= None, names="MJD_m, DM_m")

tomodel = rfn.append_fields(data, ['MJD_m', 'DM_m'],[model['MJD_m'],model['DM_m']], dtypes=None, usemask=False)

tomodel = tomodel[tomodel['edm']<0.0001]

k = np.sum(tomodel['DM']*tomodel['DM_m']/np.power(tomodel['edm'],2))/np.sum(np.power(tomodel['DM_m'],2)/np.power(tomodel['edm'],2))#k that minimises the chi2

plt.errorbar(tomodel['MJD'], tomodel['DM']-k*tomodel['DM_m'],yerr=tomodel['edm'],fmt='k.')
plt.errorbar(tomodel['MJD'], tomodel['DM']-tomodel['DM_m'],yerr=tomodel['edm'],fmt='r.')
plt.show()

data['DM'] -= model['DM_m']*k

np.savetxt("%s_residuals.dat"%(dmfile.split(".")[0]), data, fmt=[b'%s','%s','%.16f','%.16f','%.16f','%.16f']) 

