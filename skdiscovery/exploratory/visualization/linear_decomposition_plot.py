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


# linear interpolation functions and plotting
# mostly for use in tandem with other functions

import matplotlib.pyplot as plt
import numpy as np

from .vis_utils import coldict


def lin_trend(inData):
    '''
    Calculates a linear polynomial fit and evaluates

    @param inData: Input data to fit

    @return Array of evaluated points for the linear fit
    '''

    real = np.isfinite(inData)
    rng  = np.arange(len(inData))

    lin_fit = np.poly1d(np.polyfit(rng[real], inData[real], 1))(rng[real])

    return lin_fit

def calc_lin_interp(inData, iterStep = 100):
    '''
    Calculates a piecewise linear interpolated fit for some data

    @param inData: Input data to fit
    @param iterStep: Number of data points per interpolation step

    @return Array of interpolated values
    '''

    if iterStep > len(inData):
        iterStep = input('Step size too large (need < {}), enter new: '.format(len(inData)))

    intIndex = np.arange(len(inData))
    interp = np.interp(intIndex[::iterStep], intIndex, inData) 
    return interp

def plot_lin_trend(inData, plotIndex = None, show = True):
    '''
    Plots a linear linear trend against its source data

    @param inData: Input data to fit and plot
    @param plotIndex: Optional index array to pass for plotting
    @param show: Boolean to show plot immediately after plot creation
    '''

    real = np.isfinite(inData)

    if plotIndex is None:
        try:
            plotIndex = inData[real].index
        except:
            plotIndex = np.arange(len(inData[real]))

    lin_fit = lin_trend(inData)

    plt.plot(plotIndex, inData, color = coldict[0])
    plt.plot(plotIndex, lin_fit, color = coldict[0], linestyle = '--', alpha = .4)

    if show:
        plt.show()

def plot_lin_interp(inData, interps = None, plotIndex = None, iterSteps = [100], pRange = [], mainTitle = 'Piecewise Decomposition', plotReal = True, show = True):
    '''
    Plots linear interpolation against its source data

    @param inData: Input data to fit and plot
    @param interps: Optional interpolated data to be plotted, will be made if not given
    @param plotIndex: Optional index array to pass for plotting
    @param iterSteps: List of iterStep values to calculate/plot
    @param pRange: Range over which to plot, defaults to start and end of original data
    @param mainTitle: Optional string plot title
    @param plotReal: Boolean variable to optionally disinclude source data
    @param show: Boolean to show plot immediately after plot creation

    @return Multidimensional array of interpreted data values
    '''

    col = 0

    fig = plt.figure(figsize = (12, 6.75))

    if interps is None:
        interps = [calc_lin_interp(inData, step) for step in iterSteps]

    if len(pRange) != 2:
        pRange = [0, len(inData)]

    if plotIndex is None:
        try:
            plotIndex = inData.index
            if ((plotIndex.dtype == np.dtype('<M8[ns]')) or (plotIndex.dtype == np.dtype('datetime64[ns]'))):
                xlabel = 'Time'
            else:
                xlabel = 'Index'
        except:
            plotIndex = np.arange(len(inData))
            xlabel    = 'Index'
    else:
        if ((plotIndex.dtype == np.dtype('<M8[ns]')) or (plotIndex.dtype == np.dtype('datetime64[ns]'))):
             xlabel = 'Time'
        else:
             xlabel = 'Index'

    if plotReal:
        plt.plot(plotIndex[pRange[0] : pRange[1]], inData[pRange[0] : pRange[1]], color = '#CC0000', label = 'Raw data')

    for num, step in enumerate(iterSteps):
        plt.plot(plotIndex[pRange[0] : pRange[1] : step], interps[num][pRange[0] : pRange[1]], color = coldict[col % 16], label = 'Interpolated ({} points/segment)'.format(step))
        col += 1

    plt.legend(fontsize = 14)
    plt.xlabel(xlabel, fontsize = 18)
    plt.ylabel('Amplitude', fontsize = 18)
    plt.title(mainTitle, fontsize = 18)
    plt.xticks(fontsize = 18)
    plt.yticks(fontsize = 18)

    if show:
        plt.show()

    return interps


def plot_lin_slope(inData, interps = None, plotIndex = None, mainTitle = 'Piecewise Decomposition and Slopes', iterSteps = [100], pRange = [], plotReal = True, show = True):
    '''
    Plots raw data, linear interpolated data, and interpolated slope

    @param inData: Input data to fit and plot
    @param interps: Optional interpolated data to be plotted, will be made if not given
    @param plotIndex: Optional index array to pass for plotting
    @param iterSteps: List of iterStep values to calculate/plot
    @param pRange: Range over which to plot, defaults to start and end of original data
    @param mainTitle: Optional string plot title
    @param plotReal: Boolean variable to optionally disinclude source data
    @param show: Boolean to show plot immediately after plot creation   

    @return Tuple of interpolated values array and corresponding gradient array
    '''

    coli = 0; colg = 0

    if interps is None:
        interps = [calc_lin_interp(inData, step) for step in iterSteps]

    if len(pRange) != 2:
        pRange = [0, len(inData)]

    if plotIndex is None:
        try:
            plotIndex = inData.index
            if ((plotIndex.dtype == np.dtype('<M8[ns]')) or (plotIndex.dtype == np.dtype('datetime64[ns]'))):
                xlabel = 'Time'
            else:
                xlabel = 'Index'
        except:
            plotIndex = np.arange(len(inData))
            xlabel    = 'Index'
    else:
        if ((plotIndex.dtype == np.dtype('<M8[ns]')) or (plotIndex.dtype == np.dtype('datetime64[ns]'))):
             xlabel = 'Time'
        else:
             xlabel = 'Index'

    gradients = [np.gradient(interp) for interp in interps]

    fig, axes = plt.subplots(2 + plotReal, 1, sharex = True, figsize = (12, 6.75))
    
    fig.add_subplot(111, frameon = False)
    plt.tick_params(labelcolor = 'none', top = 'off', bottom = 'off', left = 'off', right = 'off')
    plt.xlabel(xlabel, fontsize = 18)
    plt.ylabel('Amplitude', fontsize = 18)

    if plotReal:
        axes[0].plot(plotIndex, inData, color = '#CC0000')
        axes[0].set_title(mainTitle, fontsize = 18)

    for num, interp in enumerate(interps):
        axes[0 + plotReal].plot(plotIndex[::iterSteps[num]], interp, color = coldict[coli % 16], label = '{} points/segment'.format(iterSteps[num]))
        axes[0 + plotReal].set_title('Interpolated', fontsize = 18)
        coli += 1
    axes[0 + plotReal].legend(fontsize = 12)

    for num, grad in enumerate(gradients):
        axes[1 + plotReal].plot(plotIndex[::iterSteps[num]], grad, color = coldict[colg % 16], label = '{} points/segment'.format(iterSteps[num]))
        axes[1 + plotReal].set_title('Interpolated slope', fontsize = 18)
        colg += 1
    axes[1 + plotReal].legend(fontsize = 12)

    for ax in axes:
        ax.tick_params(axis = 'both', labelsize = 14)
    
    if show:
        plt.show()

    return interps, gradients

"""
def calc_gradient_extrema(gradient, numExt = 3):
    '''
    NOT DONE
    '''

    #! make it calculate where on original index extrema lie
    
    arr = gradient.copy()
    maxs = {}
    mins = {}
    
    for _ in np.arange(numExt):
        maxs[np.argmax(arr)] = np.max(arr)
        mins[np.argmin(arr)] = np.min(arr)
        arr = np.delete(arr, (np.argmax(arr), np.argmin(arr)))

    extrema = np.concatenate((list(maxs.values()), list(mins.values())))
    indices = np.concatenate((list(maxs.keys()), list(mins.keys())))

    return indices, extrema
"""
