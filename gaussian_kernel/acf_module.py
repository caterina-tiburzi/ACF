#!/usr/bin/env python

import numpy as np
import os, argparse, math, sys
import matplotlib.pyplot as plt
from sklearn import preprocessing

def acf_slotting_gaussian(t,y,ye):
    
    plt.errorbar(t,y,yerr=ye'k.')
    plt.xlabel("x")
    plt.ylabel("y")
    plt.show()


    tt_r = np.vstack((t,y))
    tt = np.core.records.fromarrays(tt_r, names='time,value', formats = 'f8,f8')
	    
	
	
	
    #Paragraph 2.3 of Rehfeld et al. 2011, normalization of the time series to have mean 0 (in theory, one should also normalize the variance to 1, but it shouldn't make a difference here. Note that if one normalizes the variance, the s/he should reintroduce it later, when the ACF is calculated).
	
    meandm=np.average(y)
    tt['value'] = tt['value'] - meandm
	
    #Paragraph 2.1 of Rehfeld, the mean sampling time between the elements should be 1 --> IS THIS A PART TO REINTRODUCE LATER??
	
    mean_step = np.ediff1d(tt['time']).mean(axis=0)#note that here I am not following what is mentioned in 2.1 of Rehfeld -- where they perform a rescaling of the time span. This is because I want to end up with the real lags
	
    tau = 3*mean_step #this is 3 times more than what suggested by Rehfeld
    h = tau/4.
	
    print("Bin length:",tau)
	
    #find the boundaries of the ACF bins
    
    deltat, bins = [], []
    
    f = open("%s_weights.dat"%(infile),"w")
    g = open("%s_weights_no0.dat"%(infile),"w")
	
    for i in range(len(tt['value'])):
        for j in range(i,len(tt['value'])):
            deltat.append(tt['time'][j]-tt['time'][i])
	
	
    for k in range(0,int(math.ceil(max(deltat)/tau))+1):
        bins.append(k*tau)
    #Compute the ACF by using the slotting correlation + Gaussian kernel
    print("Number of Lags:",len(bins)-1)
	
    #return deltat,bins
	
    acf, index = [], 0 
    acf_n0, index_n0 = [], 0 #initialize also the ACF without the variance of the time series
	
    for k in xrange(0,int(math.ceil(max(deltat)/tau))+1):
        print "Lag %s of %s"%(k,len(bins)-1)
        b_k, acf_k = [], []
        b_k_n0, acf_k_n0 = [], []
        
        for i in xrange(0,len(tt['value'])):
            for j in xrange(i,len(tt['value'])):
                dt = tt['time'][j]-tt['time'][i]
                if dt != 0:
                    d = np.absolute(dt - k*tau)
                    b_kij = (1./np.sqrt(2.*np.pi*h)) * np.exp(-np.power(d,2)/(2*np.power(h,2)))
                    acf_k.append(tt['value'][i]*tt['value'][j]*b_kij)
                    b_k.append(b_kij)
                    acf_k_n0.append(tt['value'][i]*tt['value'][j]*b_kij)
                    b_k_n0.append(b_kij)
	
                else:
                    d = np.absolute(dt - k*tau)
                    b_kij = (1./np.sqrt(2.*np.pi*h)) * np.exp(-np.power(d,2)/(2*np.power(h,2)))
                    acf_k.append(tt['value'][i]*tt['value'][j]*b_kij)
                    b_k.append(b_kij)
	    
        acf.append(np.sum(acf_k)/np.sum(b_k)) #NORMALIZED ACF
        acf_n0.append(np.sum(acf_k_n0)/np.sum(b_k_n0)) #NORMALIZED ACF, WITHOUT VARIANCE
        f.write("%.3f %.10f\n"%(k*tau,np.sum(b_k)))
        g.write("%.3f %.10f\n"%(k*tau,np.sum(b_k_n0)))
	
	
    f.close()#If one wants the non-normalized ACF, these files can be used to multiply the ACF values
    g.close()
    
    return(np.asarray(acf),np.asarray(acf_n0),bins)
