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


# emd (empirical mode decomposition) and hht (hilbert huang transform) functions
#! robust up until ~1500 data points
#! rename this file?
#! need to look over wrappers, streamline

import os, sys, warnings

from matplotlib.ticker import FormatStrFormatter
from scipy import stats
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyhht

from .vis_utils import *
    # block_output(), enable_output(),


# core functionality

def calc_imfs(rawData, nbsym = False):
    '''
    IMF calculation function, streamlined and quieted

    @param rawData: Input data for EMD calculation
    @param nbsym: Boolean that would add extra data points near boundaries when calculating; breaks some datasets unless False

    @return 2D numpy.ndarray of IMFs
    '''

    block_output()
    imfs = pyhht.emd.EMD(rawData, nbsym = nbsym).decompose()
    enable_output()
    return imfs

def calc_imfs_sum(imfs, highNum = 2, high = True, residual = False):
    '''
    IMF summation helper function

    @param imfs: Input array of IMFs to be summed
    @param highNum: Number of high frequency IMFs to sum, starting from IMF1 (indexed at 0)
    @param high: Boolean that determines which class of frequency to sum (default True to sum HF)
    @param residual: Boolean that optionally includes the residual function when summing low frequency IMFs (default False to disinclude residual)

    @return 1D numpy.ndarray of summed IMFs
    '''

    if high:
        summed = np.sum(imfs[: highNum], axis = 0)
    else:
        summed = np.sum(imfs[highNum : len(imfs) - residual], axis = 0)

    return summed

def plot_imfs(rawData, imfs, toPlot = [], mainTitle = 'IMFs', show = True, figsize=(12,10)):
    '''
    Plots raw data and IMFs in a subplot grid (n Imfs [rows] x 1 [col])
    
    @param rawData: Input data for plotting
    @param imfs: Input array of IMFs for plotting
    @param toPlot: List of which IMFs to plot (default is all)
    @param mainTitle: Main title of the plot
    @param show: Boolean to show plot immediately after plot creation
    @param figsize: Size of figure
    '''

    try:
        plotIndex = rawData.index
        if ((plotIndex.dtype == np.dtype('<M8[ns]')) or (plotIndex.dtype == np.dtype('datetime64[ns]'))):
            xlabel = 'Time'
        else:
            xlabel = 'Index'

    except:
        plotIndex = np.arange(len(rawData))
        xlabel    = 'Index'

    if len(toPlot) == 0:
        toPlot = np.arange(len(imfs))
    else:
        toPlot = [i - 1 for i in toPlot]

    fig, axes = plt.subplots(len(toPlot) + 1, 1, sharex = True, figsize = figsize)

    naxis = fig.add_subplot(111, frameon = False)
    plt.tick_params(labelcolor = 'none', top = 'off', bottom = 'off', left = 'off', right = 'off')

    if len(imfs) <= 10:
        lSize = 14
    else:
        lSize = 12


    axes[0].plot(plotIndex, rawData, color = '#0066FF')
    axes[0].set_title(mainTitle, fontsize = 18)
    axes[0].set_ylabel('Raw', fontsize = 18)

    axes[-1].set_xlabel(xlabel, fontsize = 18)

    axnum = 1

    for i in toPlot:
        if i != (len(imfs) - 1):
            axes[axnum].plot(plotIndex, imfs[i], color = '#009900')
            axes[axnum].set_ylabel('IMF{}'.format(i + 1), fontsize = lSize)
            axnum += 1

        else:
            axes[axnum].plot(plotIndex, imfs[-1], color = '#000000')
            axes[axnum].set_ylabel('Resid', fontsize = lSize)
            axnum += 1

    for ax in axes:
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

        ax.tick_params(axis = 'x', labelsize = 18)
        ax.tick_params(axis = 'y', labelsize = lSize)

    if show:
        plt.show()

def plot_imfs_split(rawData, imfs, highNum = 2, residual = False, mainTitle = 'Raw data', collage = False, show = True):
    '''
    Plots raw data and summed IMFs based on HF/LF cut, can optionally plot the residual separately from LF

    @param rawData: Input data for plotting
    @param imfs: Input array of IMFs for summing and then plotting
    @param highNum: Number of high frequency IMFs to sum, starting from IMF1 (indexed at 0)
    @param residual: Boolean that optionally includes the residual function when summing low frequency IMFs (default False to disinclude residual)
    @param mainTitle: Title string of plot
    @param collage: Boolean that can optionally return certain plot parameters for external usage
    @param show: Boolean to show plot immediately after plot creation

    @return Tuple of HF summed data array and LF summed data array
    '''
    try:
        plotIndex = rawData.index
        if ((plotIndex.dtype == np.dtype('<M8[ns]')) or (plotIndex.dtype == np.dtype('datetime64[ns]'))):
            xlabel = 'Time'
        else:
            xlabel = 'Index'

    except:
        plotIndex = np.arange(len(rawData))
        xlabel    = 'Index'
    
    highData = calc_imfs_sum(imfs, highNum = highNum, residual = residual)
    lowData  = calc_imfs_sum(imfs, highNum = highNum, residual = residual, high = False)

    plotArgs = rawData, highData, lowData, imfs[-1],
    if collage:
        return plotArgs

    fig, axes = plt.subplots(3 + residual, 1, sharex = True, figsize = (12, 6.75))

    fig.add_subplot(111, frameon = False)
    plt.tick_params(labelcolor = 'none', top = 'off', bottom = 'off', left = 'off', right = 'off')
    plt.xlabel(xlabel, fontsize = 18)
    plt.ylabel('Amplitude', fontsize = 18)

    axes[0].plot(plotIndex, plotArgs[0], color = '#0066FF')
    axes[1].plot(plotIndex, plotArgs[1], color = '#CC33FF')
    axes[2].plot(plotIndex, plotArgs[2], color = '#CC0000')

    axes[0].set_title(mainTitle, fontsize = 18)
    axes[1].set_title('{} highest frequency IMFs'.format(highNum), fontsize = 18)
    if highNum == 1:
        axes[1].set_title('Highest frequency IMF', fontsize = 18)
    axes[2].set_title('{} lowest frequency IMFs and residual'.format(len(imfs) - highNum -1), fontsize = 18)
    if (len(imfs) - highNum) == 2:
            axes[2].set_title('Lowest frequency IMF and residual', fontsize = 18)

    if residual:
        axes[3].plot(plotIndex, plotArgs[3], '#000000')
            
        axes[2].set_title('{} lowest frequency IMFs'.format(len(imfs) - highNum - 1), fontsize = 18)
        if (len(imfs) - highNum) == 2:
            axes[2].set_title('Lowest frequency IMF', fontsize = 18)
        axes[3].set_title('Residual', fontsize = 18)
    
    for ax in axes:
        ax.tick_params(axis = 'x', labelsize = 12)
        ax.tick_params(axis = 'y', labelsize = 12)

    if show:
        plt.show()

    return highData, lowData


def plot_imfs_split_comp(rawData, imfs, highNums = [2, 3], residual = False, plotRaw = True, mainTitle = 'Raw data', collage = False, show = True):
    '''
    Like plot_imfs_split, plots raw data and summed IMFs based on two HF/LF cuts

    @param rawData: Input data for plotting
    @param imfs: Input array of IMFs for summing and then plotting
    @param highNums: Number of high frequency IMFs to sum and compare
    @param residual: Boolean that optionally includes the residual function when summing low frequency IMFs (default False to disinclude residual)
    @param plotRaw: Boolean to optionally disinclude raw data plot above IMF summation comparison
    @param mainTitle: Title string of plot
    @param collage: Boolean that can optionally return certain plot parameters for external usage
    @param show: Boolean to show plot immediately after plot creation

    @return Tuple of both HF summed data arrays and borh LF summed data arrays
    '''

    try:
        plotIndex = rawData.index
        if ((plotIndex.dtype == np.dtype('<M8[ns]')) or (plotIndex.dtype == np.dtype('datetime64[ns]'))):
            xlabel = 'Time'
        else:
            xlabel = 'Index'

    except:
        plotIndex = np.arange(len(rawData))
        xlabel    = 'Index'

    highData1 = calc_imfs_sum(imfs, highNum = highNums[0], residual = residual)
    highData2 = calc_imfs_sum(imfs, highNum = highNums[1], residual = residual)
    lowData1  = calc_imfs_sum(imfs, highNum = highNums[0], residual = residual, high = False)
    lowData2  = calc_imfs_sum(imfs, highNum = highNums[1], residual = residual, high = False)

    plotArgs = plotIndex, rawData, highData1, highData2, lowData1, lowData2, imfs[-1]
    if collage:
        return plotArgs

    fig = plt.figure(figsize = (12, 6.75))
    
    fig.add_subplot(111, frameon = False)
    plt.tick_params(labelcolor = 'none', top = 'off', bottom = 'off', left = 'off', right = 'off')
    plt.xlabel(xlabel, fontsize = 18)
    plt.ylabel('Amplitude', fontsize = 18)


    if plotRaw:
        ax_main = fig.add_subplot(3 + residual, 1 , 1)
        ax_main.set_title(mainTitle, fontsize = 18)
        ax_main.plot(plotArgs[0], plotArgs[1], color = '#0066FF')

    axhigh1 = fig.add_subplot(2 + plotRaw + residual, 2, 1 + 2*plotRaw)
    axhigh2 = fig.add_subplot(2 + plotRaw + residual, 2, 2 + 2*plotRaw, sharey = axhigh1)
    axlow1  = fig.add_subplot(2 + plotRaw + residual, 2, 3 + 2*plotRaw, sharex = axhigh1)
    axlow2  = fig.add_subplot(2 + plotRaw + residual, 2, 4 + 2*plotRaw, sharex = axhigh2, sharey = axlow1)

    axhigh1.set_title('{} highest frequency IMFs'.format(highNums[0]), fontsize = 16)
    axhigh2.set_title('{} highest frequency IMFs'.format(highNums[1]), fontsize = 16)
    axlow1.set_title('{} lowest frequency IMFs and residual'.format(len(imfs) - highNums[0] - 1), fontsize = 16)
    axlow2.set_title('{} lowest frequency IMFs and residual'.format(len(imfs) - highNums[1] - 1), fontsize = 16)

    axhigh1.plot(plotArgs[0], plotArgs[2], color = '#CC33FF')
    axhigh2.plot(plotArgs[0], plotArgs[3], color = '#CC33FF')
    axlow1.plot(plotArgs[0], plotArgs[4], color = '#CC0000')
    axlow2.plot(plotArgs[0], plotArgs[5], color = '#CC0000')

    axL = [ax_main, axhigh1, axhigh2, axlow1, axlow2]

    if residual:
        axresid = fig.add_subplot(2 + plotRaw + residual, 2, (7,8))

        axlow1.set_title('{} lowest frequency IMFs'.format(len(imfs) - highNums[0] - 1), fontsize = 16)
        axlow2.set_title('{} lowest frequency IMFs'.format(len(imfs) - highNums[1] - 1), fontsize = 16)
        axresid.set_title('Residual', fontsize = 16)

        axresid.plot(plotArgs[0], plotArgs[6], color = '#000000')
        axL.append(axresid)
    
    for ax in axL:
        ax.tick_params(axis = 'x', labelsize = 12)
        ax.tick_params(axis = 'y', labelsize = 12)

    plt.tight_layout()

    if show:
        plt.show()

    return highData1, highData2, lowData1, lowData2


def plot_imfs_noise(imfs, guessType = 'high', noiseNum = 2, collage = False, show = True):
    '''
    Plots assumed noise from IMF summation in a histogram, with overlaid graphs of fit probability distributions to check if assumption can be validated
    
    @param imfs: Input array of IMFs to be summed
    @param guessType: String of noise guess type ('high' or 'low' are possibilities)
    @param noiseNum: Number of IMFs to sum
    @param collage: Boolean that can optionally return certain plot parameters for external usage
    @param show: Boolean to show plot immediately after plot creation

    @return Array of plotted noise
    '''

    #! Need to fix K-S test
    
    if type(guessType) == str:
        if guessType.lower() == 'high':
            highFreq = True
        elif guessType.lower() == 'low':
            highFreq = False
        else:
            print('Invalid string argument for guessType')
            return
    else:
        ratio = guessType / len(imfs[0])
        if ratio >= 1. / 10:
            highFreq = True
        else:
            highFreq = False

    if highFreq:
        imfNoise = np.sum(imfs[: noiseNum], axis = 0)
    else:
        imfNoise = np.sum(imfs[-(noiseNum +1) : -1], axis = 0)


    mu  = np.mean(imfNoise)
    sig = np.std(imfNoise)
    

    plotArgs = imfNoise,
    if collage:
        return plotArgs

    plt.figure(figsize = (12, 6.75))
    n, bins, patches = plt.hist(imfNoise, bins = 250, normed = 1, color = '#CC33FF', label = 'Data bins')

    #? if p large (~ 1) cannot reject null; if p << 1, most likely can reject null (different distributions)
    #? kstest does not seem to be working as expected; alternative: i still dont understand what it is trying to do

    gauss   = stats.norm.pdf(bins, loc = mu, scale = sig)
    gauss_p = stats.kstest(imfNoise, 'norm', args = (mu, sig))

    tmp = '{:.2f}'.format(mu)
    if (tmp == '0.00') or (tmp == '-0.00'):
        muFmt = '{:.2}'.format(mu)
    else:
        muFmt = tmp

    tmp = '{:.2f}'.format(sig)
    if (tmp == '0.00') or (tmp == '-0.00'):
        sigFmt = '{:.2}'.format(sig)
    else:
        sigFmt = tmp

    plt.plot(bins, gauss, color = '#0066FF', linestyle = '--', label = r'Normal (fit): $\mu$ = {0}, $\sigma$ = {1}'.format(muFmt, sigFmt, gauss_p[1]))

    cGamma    = np.max(stats.cauchy.pdf(bins) / np.max(n))
    cauchyS   = stats.cauchy.pdf(bins, loc = 0, scale = cGamma)
    cauchyS_p = stats.kstest(imfNoise, 'cauchy', args = (0, cGamma))

    tmp = '{:.2f}'.format(cGamma)
    if (tmp == '0.00') or (tmp == '-0.00'):
        gamFmt = '{:.2}'.format(cGamma)
    else:
        gamFmt = tmp

    plt.plot(bins, cauchyS, color = '#FF6600', linestyle = '--', label = r'Cauchy (fit): $\gamma$ = {0}'.format(gamFmt, cauchyS_p[1]))

    plt.xlabel('Noise values')
    plt.ylabel('Probability (density)')
    plt.title('IMF noise ({})'.format(guessType))
    plt.legend(fontsize = 8)
    
    if show:
        plt.show()

    return imfNoise
    

#%# wrappers #%#

def run_plotImfs(inData, imfs = None, nbsym = False, toPlot = [], mainTitle = 'IMFs', show = True, figsize=(12,10)):
    '''
    Wrapper for plot_imfs

    @param inData: Input data for plotting
    @param imfs: Input array of IMFs for plotting
    @param nbsym: Boolean that would add extra data points near boundaries when calculating; breaks some datasets unless False
    @param toPlot: List of which IMFs to plot (default is all)
    @param mainTitle: Main title of plot
    @param figsize: Tuple containing the figure size
    @param show: Boolean to show plot immediately after plot creation

    @return Intrinsic mode functions
    '''
    real = np.isfinite(inData)
    if imfs is None:
        imfs = calc_imfs(inData[real], nbsym = nbsym)

    plot_imfs(inData[real], imfs, toPlot = toPlot, mainTitle = mainTitle, show = show, figsize = figsize)
    return imfs

def run_plotImfsSplit(inData, imfs = None, nbsym = False, highNum = 2, residual = False, mainTitle = 'Raw data', collage = False, show = True):
    # simple subplot-grid imf frequency split (high versus low) plot
    '''
    Wrapper for plot_imfs_split

    @param inData: Input data for plotting
    @param imfs: Input array of IMFs for summing and then plotting
    @param nbsym: Boolean that would add extra data points near boundaries when calculating; breaks some datasets unless False
    @param highNum: Number of high frequency IMFs to sum, starting from IMF1 (indexed at 0)
    @param residual: Boolean that optionally includes the residual function when summing low frequency IMFs (default False to disinclude residual)
    @param mainTitle: Title string of plot
    @param collage: Boolean that can optionally return certain plot parameters for external usage
    @param show: Boolean to show plot immediately after plot creation

    @return Tuple of HF summed data array and LF summed data array
    '''
    real = np.isfinite(inData)
    if imfs is None:
        imfs = calc_imfs(inData[real], nbsym = nbsym)

    outSplit = plot_imfs_split(inData[real], imfs, highNum = highNum, residual = residual, mainTitle = mainTitle, collage = collage, show = show)
    return outSplit


def run_plotImfsSplitComp(inData, imfs = None, nbsym = False, highNums = [2, 3], residual = False, plotRaw = True, mainTitle = 'Raw data', collage = False, show = True):
    '''
    Wrappper for plot_imfs_split_comp
    @param inData: Input data for plotting
    @param imfs: Input array of IMFs for summing and then plotting
    @param nbsym: Boolean that would add extra data points near boundaries when calculating; breaks some datasets unless False
    @param highNums: Number of high frequency IMFs to sum and compare
    @param residual: Boolean that optionally includes the residual function when summing low frequency IMFs (default False to disinclude residual)
    @param plotRaw: Boolean to optionally disinclude raw data plot above IMF summation comparison
    @param mainTitle: Title string of plot
    @param collage: Boolean that can optionally return certain plot parameters for external usage
    @param show: Boolean to show plot immediately after plot creation

    @return Tuple of both HF summed data arrays and borh LF summed data arrays
    '''
    real = np.isfinite(inData)
    if imfs is None:
        imfs = calc_imfs(inData[real], nbsym = nbsym)
        
    outSplit = plot_imfs_split_comp(inData[real], imfs, highNums = highNums, residual = residual, plotRaw = plotRaw, mainTitle = mainTitle, collage = collage, show = show)
    return outSplit

def run_plotImfsNoise(inData, imfs = None, nbsym = False, noiseNum = 2, guessType = 'high', show = True):
    '''
    Wrapper for plot_imfs_noise
    @param inData: Input data for plotting
    @param imfs: Input array of IMFs to be summed
    @param nbsym: Boolean that would add extra data points near boundaries when calculating; breaks some datasets unless False
    @param guessType: String of noise guess type ('high' or 'low' are possibilities)
    @param noiseNum: Number of IMFs to sum
    @param collage: Boolean that can optionally return certain plot parameters for external usage
    @param show: Boolean to show plot immediately after plot creation

    @return Array of noise data values
    '''
    real = np.isfinite(inData)
    if imfs is None:
        imfs = calc_imfs(inData[real], nbsym = nbsym)

    outNoise = plot_imfs_noise(imfs, guessType = guessType, noiseNum = noiseNum, show = show)
    return outNoise

def run_plotImfsSplitNoise(inData, imfs = None, nbsym = False, highNum = 2, residual = False, mainTitle = 'Raw data', noiseNum = 2, guessType = 'high', show = False):
    '''
    Wrapper for both plot_imfs_split and plot_imfs_noise

    @param inData: Input data for plotting
    @param imfs: Input array of IMFs for summing and then plotting
    @param nbsym: Boolean that would add extra data points near boundaries when calculating; breaks some datasets unless False
    @param highNum: Number of high frequency IMFs to sum
    @param residual: Boolean that optionally includes the residual function when summing low frequency IMFs (default False to disinclude residual)
    @param mainTitle: Title string of plot
    @param guessType: String of noise guess type ('high' or 'low' are possibilities)
    @param noiseNum: Number of IMFs to sum
    @param show: Boolean to show plot immediately after plot creation

    @return Tuple of split tuple and noise array

    '''
    real = np.isfinite(inData)
    if imfs is None:
        imfs = calc_imfs(inData[real], nbsym = nbsym)

    outSplit = plot_imfs_split(inData[real], imfs, highNum = highNum, residual = residual, mainTitle = mainTitle, show = show)
    outNoise = plot_imfs_noise(imfs, guessType = guessType, noiseNum = noiseNum, show = show)

    plt.show()

    return outSplit, outNoise
