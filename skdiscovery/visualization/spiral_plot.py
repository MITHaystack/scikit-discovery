# The MIT License (MIT)
# Copyright (c) 2017 Massachusetts Institute of Technology
#
# Authors: Evan Wojciechowski
# This software is part of the NSF DIBBS Project "An Infrastructure for
# Computer Aided Discovery in Geoscience" (PI: V. Pankratius) and 
# NASA AIST Project "Computer-Aided Discovery of Earth Surface 
# Deformation Phenomena" (PI: V. Pankratius)
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

import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .vis_utils import *
    # mod_data() for a (more) thorough checking routine than in other functions, keeping plot_spiral from breaking

def plot_spiral(plotData, plotIndex, T, mainTitle = 'Spiral plot', barLabel = 'Amplitude', plotTS = False, show = True):
    '''
    Plots data in a spiral pattern via a polar plot

    @param plotData: Input data values/amplitudes
    @param plotIndex: Input index (series time coordinates)
    @param T: Period value with which to wrap data around the plot
    @param mainTitle: Title for plot
    @param barLabel: Colorbar label
    @param plotTS: Optional flag to plot the time series of the data in a separate window
    @param show: Boolean to show plot immediately after plot creation
    '''
    r     = plotIndex - plotIndex.min()
    theta = 2 * np.pi * r / T
    z     = plotData

    fig = plt.figure()
    ax  = fig.add_axes([.1, .1, .8, .8], projection = 'polar', frameon = False)
    ax.set_title(mainTitle, fontsize = 18)

    heatMap = ax.scatter(theta, r, c = z, s = 10, edgecolors = 'none', cmap = 'plasma')
    cbar = fig.colorbar(heatMap)
    cbar.set_label(barLabel, fontsize = 18)
    cbar.ax.tick_params(labelsize = 12)

    ax.set_rmax(max(r))
    ax.set_rlabel_position(-22.5)
    ax.set_xticks(np.pi/180. * np.linspace(0, 360, 8, endpoint = False))
    ax.grid(True)

    ax.text(0, 0, 'Period = {:.4f}'.format(T), fontsize = 14, transform = ax.transAxes)

    if plotTS: # only works for kepler, basically broken for others
        plt.figure()
        plt.plot(plotIndex % T, plotData, '.')
    
    if show:
        plt.show()

#%# wrappers #%#
def run_spiral(inData, period, inIndex = None, mainTitle = 'Spiral plot', barLabel = 'Amplitude', plotTS = False, show = True):
    '''
    Wrapper for plot_spiral

    @param inData: Input data to use in plot
    @param period: Period value with which to wrap data around the plot
    @param inIndex: Input index (series time coordinates)
    @param mainTitle: Title for plot
    @param barLabel: Colorbar label
    @param plotTS: Optional flag to plot the time series of the data in a separate window
    @param show: Boolean to show plot immediately after plot creation
    '''

    plotData, plotIndex = mod_data(inData, inIndex, makeType = int)
    plot_spiral(plotData, plotIndex, T = period, mainTitle = mainTitle, barLabel = barLabel, plotTS = plotTS, show = show)

def run_spiralInteractive(inData, period, pParams = [], inIndex = None, mainTitle = 'Spiral plot', barLabel = 'Amplitude', plotTS = False):
    '''
    Wrapper for plot_spiral that is interactive when used in Jupyter notebooks

    @param inData: Input data to use in plot
    @param period: Period value with which to wrap data around the plot
    @param pParams: List of plot's period parameters [min, max, step] necessary for interactive
    @param inIndex: Input index (series time coordinates)
    @param mainTitle: Title for plot
    @param barLabel: Colorbar label
    @param plotTS: Optional flag to plot the time series of the data in a separate window
    '''

    from ipywidgets import fixed, interact, interactive
    from IPython.display import clear_output, display, HTML
    
    plotData, plotIndex = mod_data(inData, inIndex)

    if len(pParams) != 3:
        pParams = [period - .1, period + .1, .01]

    inter = interactive(plot_spiral, plotData = fixed(plotData), plotIndex = fixed(plotIndex), mainTitle = fixed(mainTitle), barLabel = fixed(barLabel), plotTS = fixed(plotTS), T = (pParams[0], pParams[1], pParams[2]))
    display(inter)
