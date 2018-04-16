
#!/usr/bin/env python

import numpy as np
import os, argparse, math, sys
import matplotlib.pyplot as plt
from acf_module import acf_slotting_gaussian
from acf2ps_module import acf2ps

parser = argparse.ArgumentParser(description="Compute ACF")
parser.add_argument('-f', '--dmfile', type=str, nargs=1,
                    help="File containing DM variations, as observation, date, MJD, DM, edm, angle")
parser.add_argument('--eb_cut', type=float, nargs='?', default=0.0001,
                    help="error bar threshold in the DM plot")
parser.add_argument('--lag_cut', type=int, nargs='?', default=1200,
                    help="Tricky parameter. It defines the ")
parser.add_argument('--a_cut', type=float, nargs='?',
                    help="Window on the Solar angle")

args = parser.parse_args()
dmfile = args.dmfile[0]
if (args.eb_cut):
    cut = args.eb_cut
if (args.a_cut):
    acut = args.a_cut
if (args.lag_cut):
    lc = args.lag_cut

data = np.genfromtxt(dmfile, dtype= None, names="observation, date, MJD, DM, edm, angle")
data = np.sort(data, order='MJD')
data = data[data['edm'] < cut]
if (args.a_cut):
    data = data[data['angle'] > acut]

#acf, acf_n0, bins =acf_slotting_gaussian(data['MJD'],data['DM'],data['edm'])
#
#if (args.a_cut):
#    np.savetxt("%s_angle%s_norm-acf.dat"%(dmfile.split(".")[0],int(acut)),acf)
#    np.savetxt("%s_angle%s_bins.dat"%(dmfile.split(".")[0],int(acut)),bins)
#    np.savetxt("%s_angle%s_norm-acf_n0.dat"%(dmfile.split(".")[0],int(acut)),acf_n0)
#
#else:
#    np.savetxt("%s_norm-acf.dat"%(dmfile.split(".")[0]),acf)
#    np.savetxt("%s_bins.dat"%(dmfile.split(".")[0]),bins)
#    np.savetxt("%s_norm-acf_n0.dat"%(dmfile.split(".")[0]),acf_n0)
#
#
#plt.plot(bins,acf,'r-',label='with 0-lag cp')
#plt.plot(bins,acf_n0,'k-', label='without 0-lag cp')
#plt.xlabel("Lags [d]")
#plt.ylabel("ACF")
#plt.title("ACF for %s"%(dmfile))
#
#if (args.a_cut):
#    plt.savefig("%s_angle%s_ACF.png"%(dmfile.split(".")[0],int(acut)))
#else:
#    plt.savefig("%s_ACF.png"%(dmfile.split(".")[0]))
#
#plt.show()


if (args.a_cut): 
    acf = np.genfromtxt("%s_angle%s_norm-acf.dat"%(dmfile.split(".")[0],int(acut)),names='acf')
    acf_n0 = np.genfromtxt("%s_angle%s_norm-acf_n0.dat"%(dmfile.split(".")[0],int(acut)),names='acf')
    bins = np.genfromtxt("%s_angle%s_bins.dat"%(dmfile.split(".")[0],int(acut)),names='lags')
else:
    acf = np.genfromtxt("%s_norm-acf.dat"%(dmfile.split(".")[0]),names='acf')
    acf_n0 = np.genfromtxt("%s_norm-acf_n0.dat"%(dmfile.split(".")[0]),names='acf')
    bins = np.genfromtxt("%s_bins.dat"%(dmfile.split(".")[0]),names='lags')


plt.plot(bins,acf,'k-')
plt.plot(bins,acf_n0,'r-')
plt.xlabel("Lags")
plt.ylabel("ACF")
plt.show()


corrected_acf,f,ps = acf2ps(bins,acf_n0,lc)

if (args.a_cut):
    np.savetxt("%s_angle%s_corrected-acf.dat"%(dmfile.split(".")[0],int(acut)),acf)
    np.savetxt("%s_angle%s_freq.dat"%(dmfile.split(".")[0],int(acut)),bins)
    np.savetxt("%s_angle%s_ps.dat"%(dmfile.split(".")[0],int(acut)),acf_n0)

else:
    np.savetxt("%s_corrected-acf.dat"%(dmfile.split(".")[0]),acf)
    np.savetxt("%s_freq.dat"%(dmfile.split(".")[0]),bins)
    np.savetxt("%s_ps.dat"%(dmfile.split(".")[0]),acf_n0)



plt.loglog(f,ps,'k-')
plt.axvline(x=1./(356.25),color='r',label=r'1/yr')
plt.xlabel("Frequency [cycles/day]")
plt.ylabel("PS")
plt.title("PS for %s"%(dmfile))
plt.legend(loc='best')
if (args.a_cut):
    plt.savefig("%s_angle%s_PS.png"%(dmfile.split(".")[0],int(acut)))
else:
    plt.savefig("%s_PS.png"%(dmfile.split(".")[0])) 
plt.show()                                                                 


