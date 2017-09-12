# The MIT License (MIT)
# Copyright (c) 2017 Massachusetts Institute of Technology
#
# Authors: Victor Pankratius, Justin Li, Cody Rude
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

from skdiscovery.data_structure.framework import PipelineItem
import numpy as np
import matplotlib.pyplot as plt
import math

class Plotter(PipelineItem):
    '''
    Make a plot of series data
    '''
    def __init__(self, str_description, num_columns = 3, errorbars=False, width=13, height=4, **kwargs):
        '''
        Initialize Plotter
        
        @param str_description: String describing accumulator
        @param num_columns: Number of columns to use when plotting data
        @param errorbars: Flag indicating if errorbars should be used
        @param width: Total width of all columns combined
        @param height: Height of single row of plots
        @param **kwargs: Any additional keyword arguments are passed on to matplotlib
        '''
        self.kwargs = kwargs
        self.num_columns = num_columns
        self.errorbars = errorbars
        self.height = height
        self.width = width
        super(Plotter, self).__init__(str_description, [])
    
    def process(self, obj_data):
        '''
        Plot each column in obj_data

        @param obj_data: Data Wrapper
        '''

        width = self.width
        height = self.height
        
        labels = []
        data = []
        errors = []
        for label, series, err in obj_data.getIterator():
            data.append(series)
            errors.append(err)
            labels.append(label)

        # if len(data) > self.num_columns:
        #     width *= self.num_columns
        # else:
        #     width *= len(data)

        if len(labels) > 0:

            rows = math.ceil(len(labels) / self.num_columns)
            height *= rows

            figure = plt.figure()
            figure.set_size_inches(width, height, True)


            for label, series, err, num in zip(labels, data, errors, range(1,len(labels)+1)):
                plt.subplot(rows, self.num_columns, num)
                plt.title(series.name)
                plt.ylabel(label)

                plt.xticks(rotation=45)

                if self.errorbars:
                    plt.errorbar(np.array(series.index),np.array(series), yerr=np.array(err), **self.kwargs)
                else:
                    plt.plot(series, **self.kwargs)

            plt.tight_layout()

            if(obj_data.run_id > -1):
                figure.suptitle( "Run: " + str(obj_data.run_id), y=1.02)

