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

from .vis_utils import *

import matplotlib.pyplot as plt
import numpy as np

def calc_DFT(t, y):
    '''
    Calculates discrete Fourier transform using np.fft.fft

    @param t: Time array
    @param y: Y (data amplitude) array

    @return Tuple of post-FT frequencies and coefficients
    '''

    tf = np.linspace(t[0], 1/(2*t[-1]/len(t)), len(t)/2)
    yf = np.fft.fft(y)

    return tf, yf

def plot_DFT(tIndex, yData, collage = False, show = True, suptitle = '', hori = True):
    '''
    Plots input data and Fourier transformed coefficients in a subplot grid

    @param tIndex: Input time index for series
    @param yData: Input data amplitude
    @param collage: Boolean that can optionally return certain plot parameters for external usage
    @param show: Boolean to show plot immediately after plot creation   
    @param suptitle: Optional string to add as a plot title
    @param hori: Boolean that optionally changes the orientation of the subplot configuration
    '''

    tfIndex, yfData = calc_DFT(tIndex, yData)

    plotArgs = tIndex, yData, tfIndex, (2/len(tIndex))*np.abs(yfData[:len(tIndex)//2])

    if collage:
        return plotArgs # allows use in some other function

    if hori:
        fig, axes = plt.subplots(1, 2, figsize = (12, 6.75))
    else:
        fig, axes = plt.subplots(2, 1, figsize = (12, 6.75))
    
    if suptitle != '':
        fig.suptitle(suptitle, fontsize = 18)

    axes[0].plot(plotArgs[0], plotArgs[1], color = '#0066FF')
    axes[1].plot(plotArgs[2], plotArgs[3], color = '#FF6600')

    axes[0].set_title('Time series', fontsize = 18)
    axes[0].set_xlabel('Data Index', fontsize = 16)
    axes[0].set_ylabel('Data Amplitude', fontsize = 16)
    axes[1].set_title('Fourier spectrum', fontsize = 16)
    axes[1].set_xlabel('Frequency', fontsize = 16)
    axes[1].set_ylabel('Frequency Amplitude', fontsize = 16)

    for ax in axes:
        ax.tick_params(axis = 'both', labelsize = 12)

    if show:
        plt.show()

#%# wrappers #%#
def run_plotDFT(inData, inIndex = None, collage = False, show = True, suptitle = '', hori = True):
    '''
    Wrapper for plot_DFT

    @param inData: Input data for plotting
    @param inIndex: Possible input index to use in calculating DFT
    @param collage: Boolean that can optionally return certain plot parameters for external usage
    @param show: Boolean to show plot immediately after plot creation   
    @param suptitle: Optional string to add as a plot title
    @param hori: Boolean that optionally changes the orientation of the subplot configuration
    '''


    y, t = mod_data(inData, inIndex, makeType = 'fourier')
    plot_DFT(t, y, collage = collage, show = show, suptitle = suptitle, hori = hori)

