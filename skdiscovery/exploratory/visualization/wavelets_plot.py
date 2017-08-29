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


import matplotlib.pyplot as plt
import numpy as np
import pywt

from .vis_utils import *
    # block_output(), enable_output(), index_scale(),

# core functions

def calc_wp_deconstruct(calcData, wavelet = None):
    """
    simple function to calculate a wavelet deconstruction
    """

    if wavelet is None:
        wavelet = input('Wavelet not given, please enter a name (ex. haar): ').lower()

    #block_output() # necessary?
    deconPacket = pywt.WaveletPacket(data = calcData, wavelet = wavelet)
    #enable_output()

    return deconPacket

def calc_wp_reconstruct(deconPacket = None, calcData = None, wavelet = None, reconNodes = []):

    if deconPacket is None:
        if calcData is not None:
            deconPacket = calc_wp_deconstruct(calcData = calcData, wavelet = wavelet)
        else:
            print('Need one of deconPacket OR calcData')
            return

    reconPacket = pywt.WaveletPacket(data = None, wavelet = deconPacket.wavelet)

    if len(reconNodes) == 0:
        reconNodes = np.arange(1, deconPacket.maxlevel + 1)

    for n in reconNodes:
        reconPacket['a' * n] = deconPacket['a' * n]
        reconPacket['d' * n] = deconPacket['d' * n]

    reconPacket.reconstruct()

    return reconPacket


def plot_wp_deconstruct(deconPacket, deconNodes = [], mainTitle = 'Wavelet Deconstruction', plotRaw = True, show = True):

    col = 0

    if len(deconNodes) == 0:
        deconNodes = np.arange(1, deconPacket.maxlevel + 1)

    plt.figure(figsize = (12, 6.75))

    if plotRaw:
        plt.plot(deconPacket.data, color = coldict[col % 16], label = 'Raw')
        col += 1

    for n in deconNodes:
        indScaled = index_scale(toScale = np.arange(len(deconPacket['a' * n].data)), endRange = [0, len(deconPacket.data)])
        plt.plot(indScaled, deconPacket['a' * n].data / (np.sqrt(2) ** (n + 1)), color = coldict[col % 16], label = 'a.{}'.format(n))
        col += 1

    plt.legend(loc = 'best', fontsize = 12)
    plt.xticks(fontsize = 18)
    plt.yticks(fontsize = 18)
    plt.xlabel('Integer index from Raw (nodes scaled)', fontsize = 18)
    plt.ylabel('Amplitude', fontsize = 18)
    plt.title(mainTitle, fontsize = 18)

    if show:
        plt.show()

def plot_wp_showall(deconPacket, deconNodes = [], mainTitle = 'Individual Wavelet Nodes', plotRaw = True, show = True):

    col = 0

    if len(deconNodes) == 0:
        deconNodes = np.arange(1, deconPacket.maxlevel + 1)

    fig, axes = plt.subplots(len(deconNodes) + plotRaw, 1, sharex = True, figsize = (12, 6.75))

    if plotRaw:
        axes[0].plot(deconPacket.data, color = coldict[col % 16])
        axes[0].set_ylabel('Raw', fontsize = 18)
        col += 1

    for n in deconNodes:
        indScaled = index_scale(toScale = np.arange(len(deconPacket['a' * n].data)), endRange = [0, len(deconPacket.data)])
        axes[n - (not plotRaw)].plot(indScaled, deconPacket['a' * n].data / (np.sqrt(2) ** (n + 1)), color = coldict[col % 16])
        axes[n - (not plotRaw)].set_ylabel('a.{}'.format(n), fontsize = 18)
        col += 1

    for ax in axes:
        ax.tick_params(axis = 'x', labelsize = 18)
        ax.tick_params(axis = 'y', labelsize = 14)

    plt.xlabel('Integer index from Raw (nodes scaled)', fontsize = 18)
    axes[0].set_title(mainTitle, fontsize = 18)
    #plt.suptitle(mainTitle)

    if show:
        plt.show()
    

def plot_wp_reconstruct(reconPacket, calcData, mainTitle = 'Wavelet Reconstruction', plotRaw = True, show = True):

    col = 0

    plt.figure(figsize = (12, 6.75))

    plotArgs = calcData, reconPacket.data
    
    if plotRaw:
        plt.plot(calcData, color = coldict[col], label = 'Raw')
        col += 1

    plt.plot(reconPacket.data, color = coldict[col], label = 'Reconstruction ({} nodes)'.format(reconPacket.maxlevel))

    plt.title(mainTitle)
    plt.legend()

    if show:
        plt.show()


# wrappers

def run_plotWPDecon(inData, wavelet = None, deconNodes = [], mainTitle = 'Wavelet Deconstruction', plotRaw = True, show = True,):
    outData, outIndex = mod_data(inData)

    deconPacket = calc_wp_deconstruct(outData, wavelet = wavelet)

    plot_wp_deconstruct(deconPacket, deconNodes = deconNodes, mainTitle = mainTitle, plotRaw = plotRaw, show = show,)


def run_plotWPRecon(inData, wavelet = None, reconNodes = [], mainTitle = 'Individual Wavelet Nodes', plotRaw = True, show = True,):
    outData, outIndex = mod_data(inData)

    deconPacket = calc_wp_deconstruct(outData, wavelet = wavelet)
    reconPacket = calc_wp_reconstruct(deconPacket, calcData = inData, wavelet = wavelet, reconNodes = reconNodes)

    plot_wp_reconstruct(reconPacket, calcData = outData, mainTitle = mainTitle, plotRaw = plotRaw, show = show,)

def run_plotWPShowall(inData, wavelet = None, deconNodes = [], mainTitle = 'Wavelet Reconstruction', plotRaw = True, show = True):
    outData, outIndex = mod_data(inData)

    deconPacket = calc_wp_deconstruct(outData, wavelet = wavelet)

    plot_wp_showall(deconPacket, deconNodes = deconNodes, mainTitle = mainTitle, plotRaw = plotRaw, show = show)
