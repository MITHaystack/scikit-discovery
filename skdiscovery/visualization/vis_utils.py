# The MIT License (MIT)
# Copyright (c) 2017 Massachusetts Institute of Technology
#
# Authors: Evan Wojciechowski
# This software is part of projects sponsored by NASA and NSF (PI: V. Pankratius)

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


# written in python 3.6 by evan wojciechowski
# miscellaneous utilities for use in various time-series decomposition functions
# kept here in an attempt to reduce clutter
#! rename this (and then fix imports in dependencies)

import os, sys, warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# immutable stuff
types = [int, float, complex, np.float32, np.float64, np.int32, np.int64, np.complex64, np.complex128,]
coldict = {0 : 'C0', 1 : 'C1', 2 : 'C2', 3 : 'C3', 4 : 'C4', 5 : 'C5', 6 : 'C6', 7 : 'C7', 8 : '#92C7Ed', 9 : '#FFBB80', 10 : '#9BE49B', 11 : '#EB9393', 12 : '#C0A6D8', 13 : '#D2B3AC', 14 : '#E995D0', 15 : '#BFBFBF'}

# helper functions

def lin_trend(inData, toReturn = 'eval'):
    real = np.isfinite(inData)
    rng  = np.arange(len(inData))

    lin_fit = np.poly1d(np.polyfit(rng[real], inData[real], 1))(rng[real])

    return lin_fit


def index_scale(toScale, endRange = []):
    """
    quick linear adjustment to some input index
    primary: wavelet transforms (pywt)
    @param toScale: input index array-like to be scaled
    @param endRange: optional input list to have (at least) two values that bound output index range; assumes min is [0], max is [1]
    """

    if len(endRange) is 0:
        endRange = [0, len(toScale)]

    iMin = np.min(toScale)
    iMax = np.max(toScale)

    if (iMax - iMin) is not 0:
        delta = iMax - iMin
    else:
        delta = 1 / iMax

    return np.array(((toScale - iMin) / delta) * (endRange[1] - endRange[0]) + endRange[0])

def block_output():
    """
    blocks print output from functions and suppresses warnings
        mostly used to make output look nice, ignoring (often known) minor errors
    primary: emd (pyhht)
    """

    sys.stdout = open(os.devnull, 'w')
    warnings.filterwarnings('ignore')

def enable_output():
    """
    re-enables output, foil of blockOutput
    primary: emd (pyhht)
    """
    sys.stdout = sys.__stdout__
    warnings.filterwarnings('default')


def mod_data(inData, inIndex = None, makeType = None,):
    """
    modifies data for run_spiral so that plotted data is uniform
    @param inData: data values to be used as intensity
    @param inIndex: data values to be used as radial and angular components (once period is applied)
    @param makeType: variable tied to makeIndex which will create different types of index if necessary
    """

    makeIndex = False

    # if a separate index is provided
    if inIndex is not None:
        if len(inIndex) != len(inData):
            # sanity checks on index length
            if len(inIndex) > len(inData):
                # uses a matching amount of indices
                outIndex = inIndex[:len(inData)]
            else:
                # not enough indices, must make some new index
                makeIndex = True

    # data checking; np.ndarray and pd.core.series.Series should be the two most common types of data structure
    if isinstance(inData, np.ndarray):
        # no chance for index, only contains amplitude
        outData   = inData
        makeIndex = True 

    elif isinstance(inData, pd.core.series.Series):
        inIndex = inData.index
        #! need to check index type
        dtype = inIndex.dtype

        if dtype in types:
            # no conversion
            #! still may need to do something
            pass

        elif dtype == 'datetime64[ns]':
            # converts to np.float64, starts at zero (linear removal), unit is days
            inIndex  = inIndex.to_julian_date()
            inIndex -= inIndex.min()

        else:
            print('format not accounted for (FIX)')
            makeIndex = True

        outData  = inData.values.astype('float64')
        outIndex = inIndex

    else:
        print('format unaccounted for (FIX)')
        return

    if inIndex is None:
        # final check to see if any index was assigned
        makeIndex = True

    # creates an index if needed
    #! very rough sketch
    if makeIndex:
        if (makeType == int) or (makeType is None):
            outIndex = np.arange(len(inData))
        elif makeType.lower() == 'fourier':
            outIndex = np.linspace(0, 2*np.pi, len(inData))

    # eliminates potential NaN values from plotted data
    real     = np.isfinite(outData)
    outData  = outData[real]
    outIndex = outIndex[real]

    # returns final product for plotting 
    return outData, outIndex

