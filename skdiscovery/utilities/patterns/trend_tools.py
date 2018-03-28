# The MIT License (MIT)
# Copyright (c) 2017 Massachusetts Institute of Technology
#
# Author: Justin D Li
# This software is part of the NASA AIST Project "Computer-Aided Discovery of Earth Surface Deformation Phenomena" , PI: V. Pankratius
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.


# '''
# @package trendTools
#   This module is designed to provide a suite of tools for quick analysis of
#   the linear and sinusoidal (annual, semi-annual, seasonal, and monthly)
#   trends of time-series data (formatted using the pandas format).
#'''


import numpy as np
import pandas as pd
from scipy import optimize
from scipy.ndimage import median_filter


def getTrend(xdata):
    '''
    The getTrend function applies the signal.detrend function

    Returns the trend, given a time index input.
    
    @param xdata: 1D time-series data in a pandas series format
    
    @return the detrended data in pandas series format
    @return the linear trend assuming a 1 day per sample time fit
    @return the parameters for the linear trend
    '''
    
    # returns detrended data and the trend
    # checks the time frequency to create a time vector
    if xdata.index.freq != pd.tseries.offsets.Day():
        print('Frequency of Sampling is not Day, please check that resampling occurred!')
    timeDat = np.arange(0,len(xdata.index))
    dataDat = xdata
    # remove time and data with NaN's for now
    ninds = np.array(pd.isnull(dataDat)==True)
    argtd = (timeDat[ninds==False],dataDat[ninds==False])
    
    fitfunc = lambda p, x: p[0]*x+p[1]          # Target function (linear)
    errfunc = lambda p, x, y: fitfunc(p, x) - y # Distance to the target function
    pAnnual = [1., 5.]                          # Initial guess for the parameters
    pA, success = optimize.leastsq(errfunc, pAnnual[:], args=argtd)
    trend = fitfunc(pA, timeDat)

    # puts the data back into the pandas series format to recover the time index
    xdetr = pd.core.series.Series(data=dataDat-trend,index=xdata.index)
    
    return xdetr, trend, pA


def sinuFits(xdata,fitN=2,rmve=1):
    '''
    The sinuFits function fits annual and semi-annual sinusoid trends

    Other options allow for a monthly and seasonal sinusoid fit. The
    data is expected to be in pandas format

    @param xdata: 1D time-series data in a pandas series format
    @param fitN: the number of sinusoids to fit. 1-annual, 2-semi-annual, 3-seasonal, 4-monthly
    @param rmve: a flag to return sinusoid removed data, or the sinusoids

    @return retrDat: the returned data, either sinusoid removed or the sum of the sinusoids
    '''
    
    # checks the time frequency to create a time vector
    if xdata.index.freq != pd.tseries.offsets.Day():
        print('Frequency of Sampling is not Day, please check that resampling occurred!')
    timeDat = np.arange(0,len(xdata.index))
    # ignores time and data with NaN's for now
    ninds = np.array(pd.isnull(xdata)==True)
    argtd = (timeDat[ninds==False],xdata[ninds==False])
    if fitN >= 1:
        # annual (365 day period) sinusoidal fit
        fitfunc = lambda p, x: p[0]*np.cos(2*np.pi/365*x+p[1]) + p[2]
        errfunc = lambda p, x, y: fitfunc(p, x) - y
        pAnnual = [5., 0., 0.]
        pA, success = optimize.leastsq(errfunc, pAnnual[:], args=argtd)
        retrDat = fitfunc(pA, timeDat)
    if fitN >= 2:
        # semi-annual (182.5 day period) sinusoidal fit
        fitfunc = lambda p, x: p[0]*np.cos(2*np.pi/182.5*x+p[1]) + p[2]
        errfunc = lambda p, x, y: fitfunc(p, x) - y
        pAnnual = [5., 0., 0.]
        pSA, success = optimize.leastsq(errfunc, pAnnual[:], args=argtd)
        retrDat += (fitfunc(pSA, timeDat))
    if fitN >= 3:
        # seasonal / quarterly (91.25 day period) sinusoidal fit
        fitfunc = lambda p, x: p[0]*np.cos(2*np.pi/91.25*x+p[1]) + p[2]
        errfunc = lambda p, x, y: fitfunc(p, x) - y
        pAnnual = [5., 0., 0.]
        pSS, success = optimize.leastsq(errfunc, pAnnual[:], args=argtd)
        retrDat += (fitfunc(pSS, timeDat))
    if fitN >= 4:
        # monthly (30.5 day period) sinusoidal fit
        fitfunc = lambda p, x: p[0]*np.cos(2*np.pi/30.5*x+p[1]) + p[2]
        errfunc = lambda p, x, y: fitfunc(p, x) - y
        pAnnual = [5., 0., 0.]
        pM, success = optimize.leastsq(errfunc, pAnnual[:], args=argtd)
        retrDat += (fitfunc(pM, timeDat))
        
    if rmve==1:
        # if remove flag (rmve) is true, returns the data with the trend removed
        # otherwise, it returns the estimated overall, all components summed, trend
        retrDat = xdata - retrDat
            
    return retrDat
    
    
    
def interpNaN(data):
    '''
    Interpolate data using a linear interpolation

    @param data: 1d numpy or pandas Series with possible NaN's
    @return data after interpolation
    '''

    if isinstance(data, np.ndarray):
        data = pd.Series(data)
        return data.interpolate().as_matrix()

    elif isinstance(data, pd.Series):
        return data.interpolate()



def medianFilter(data, window, interpolate=True):
    '''
    A median filter

    If interpolate is True, data will be interpolated before smoothering.
    Otherwise, all available data within the window will be used
    
    @param data: Input data
    @param window: Size of filter window
    @param interpolate: Interpolate data before smoothing

    @return Smoothed data
    '''

    if interpolate == True:
        data = interpNaN(data)
        result = pd.Series(median_filter(data, size=window), index=data.index)

    else:
        result = data.copy()
        for index, value in data.iteritems():
            if not pd.isnull(value):
                result.loc[index] = np.nanmedian(data[np.logical_and(data.index > index-window/2,
                                                                     data.index < index+window/2)])

    return result

def normalize(in_data):
    in_data = (in_data - np.mean(in_data))
    return in_data / np.std(in_data, ddof=1)
