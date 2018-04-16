#!/usr/bin/env python

import numpy as np
import os, argparse, math, sys
import matplotlib.pyplot as plt
from scipy import signal
from scipy.fftpack import fft, fftshift
import numpy.lib.recfunctions as rfn



def reverse(ar):#reverse the array to feed the FFT algorithm
        rar = np.fft.fftshift(ar)
        if len(ar) % 2 != 0:
	        rar = np.roll(rar,1)
        return rar
	
def pos_freq(ar):#select the positive part of the array 
        if len(ar) % 2 == 0:
	        par = ar[:int(len(ar)/2)]
        else:
                par = ar[:int(len(ar)/2)+1]
        return(par)



parser = argparse.ArgumentParser(description="Compute PS from ACF")

parser.add_argument('-b', '--binfile', type=str, nargs=1, help="File containing the bins of the ACF")
parser.add_argument('-a', '--acffile', type=str, nargs=1, help="File containing the normalized ACF")

parser.add_argument('-lc', '--lag_cut', type=float, nargs=1, help="Last useful lag (i.e., where the ACF goes to 0)")

args = parser.parse_args()

lagfile = args.binfile[0]
acffile = args.acffile[0]

lc = args.lag_cut[0]

lags = np.genfromtxt(lagfile,names='lags')
acf = np.genfromtxt(acffile,names='acf')

data = lags
data = rfn.append_fields(data, ['acf'], [acf['acf']], dtypes=None, usemask=False)

	
binsize=data['lags'][1]-data['lags'][0]
lc = int(lc/binsize)
triangular = np.hstack((1-np.arange(0.,lc)/float(lc),np.zeros(len(data)-lc)))#signal.triang(len(data))
	
corrected_acf = acf['acf'] * triangular

#plt.plot(data['lags'],corrected_acf)
#plt.show()

#np.savetxt("corrected_acf.dat",corrected_acf)

corrected_acf_r =  (np.fliplr([corrected_acf])[0])[:-1]#rearrange the ACF as the FFT algorithm likes. If given in this way, the FFT algorithm will recognize it as a real, symmetric fct and it will give in output a real spectrum. The imaginary part will be minimal
acf = np.hstack((corrected_acf,corrected_acf_r))

ps = np.real(np.fft.fft(acf))*dt*2# multiplication by dt to keep into account the trapezoidal law, and by 2 to keep into account the negative part of the spectrum that we will discard later
posps = pos_freq(ps)

bins = data['lags']
dt = np.mean(np.diff(bins))
fr = np.fft.fftfreq(acf.shape[-1],d=dt)
posfr = pos_freq(fr)
#the step above can be done manually as these two lines. Given that the ACF is real and symmetric, the DFT is also symmetric and real. The negative frequencies don't have any sense for us. So we discard the negative frequencies (this is why the harmonics cover only half of the length of the ACF passed to the fft in the next two lines). To put the harmonics in place, we divide them by the total time spanned by the ACF, given by dt* (2*len(half of the ACF). Note that the last harmonic is the Nyquist frequency 
#fr = np.arange(0,len(corrected_acf)) 
#posfr = fr/(2*len(corrected_acf)*dt)

plt.loglog(posfr,posps,'k-')
plt.axvline(x=1./(356.25),color='r',label=r'1/yr')
plt.legend(loc='best')
plt.xlabel("Frequency [cycles/day]")
plt.ylabel("PS")
plt.title("PS for the corrected ACF")


plt.savefig("PS.png")


plt.show()

#psfinal = np.transpose(np.vstack((posf,posps)))
#np.savetxt("powerspectrum.dat",psfinal)

