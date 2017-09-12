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
    Make a plot of table data
    '''
    def __init__(self, str_description, column_names=None, error_column_names = None, num_columns = 3, width=13, height=4, columns_together=False,
                 annotate_column = None, annotate_data = None, xlim = None, ylim = None, **kwargs):
        '''
        Initialize Plotter
        
        @param str_description: String describing accumulator
        @param column_names: Columns to be plot
        @param error_column_names Columns containing uncertainties to be plot, no errorbars if None
        @param num_columns: Number of columns to use when plotting data
        @param width: Total width of all columns combined
        @param height: Height of single row of plots
        @param columns_together: If true, plot the columns on the same graph
        @param annotate_column: Column of annotation data to use for annotation
        @param annotate_data: Annotation data
        @param xlim: The x limit
        @param ylim: The y limit
        @param **kwargs: Any additional keyword arguments are passed on to matplotlib
        '''
        self.xlim = xlim
        self.ylim = ylim
        self.kwargs = kwargs
        self.num_columns = num_columns
        self.height = height
        self.width = width
        self.column_names = column_names
        self.annotate_column = annotate_column
        self.annotate_data = annotate_data
        
        self.error_column_names = error_column_names
        self.columns_together = columns_together
        super(Plotter, self).__init__(str_description, [])
    
    def process(self, obj_data):
        '''
        Plot each column in obj_data

        @param obj_data: Data Wrapper
        '''

        if self.column_names == None:
            column_names = obj_data.getDefaultColumns()
        else:
            column_names = self.column_names

        width = self.width
        height = self.height

        # Determine total number of figures needed
        if self.columns_together == False:
            num_figures = obj_data.getLength() * len(column_names)
        else:
            num_figures =  obj_data.getLength()
        if num_figures > 0:

            # Determine number of rows and height needed to plot all figures
            rows = math.ceil( num_figures / self.num_columns)
            height *= rows

            figure = plt.figure()
            figure.set_size_inches(width, height, True)

            if self.xlim != None:
                plt.xlim(*self.xlim)

            if self.ylim != None:
                plt.ylim(*self.ylim)
            
            num = 0
            # Main loop that iterates over all data
            for label, data in obj_data.getIterator():
                if self.columns_together == True:
                    num += 1
                    
                # Plotting with errorbars
                if self.error_column_names != None:
                    for column, err_column in zip(column_names, self.error_column_names):
                        if self.columns_together == False:
                            num += 1

                        plt.subplot(rows, self.num_columns, num)
                        plt.title(label)
                        plt.ylabel(column)

                        plt.xticks(rotation=45)
                        plt.errorbar(np.array(data.index),np.array(data[column]), yerr=np.array(data[err_column]), **self.kwargs)

                        if self.annotate_column is not None:
                            try:
                                for vline in self.annotate_data[label][self.annotate_column]:
                                    plt.axvline(vline,color='black',linewidth=3,alpha=0.5)
                            except KeyError:
                                pass
                                # print('cannot find info')
                        
                        elif self.annotate_data is not None:
                            try:
                                for vline in self.annotate_data[label]:
                                    plt.axvline(vline,color='black',linewidth=3,alpha=0.5)
                            except KeyError:
                                pass
                                # print('cannot find info')

                # Plotting without errorbars
                else:
                    for column in column_names:
                        if self.columns_together == False:
                            num += 1
                        
                        plt.subplot(rows, self.num_columns, num)
                        plt.title(label)
                        plt.ylabel(column)

                        plt.xticks(rotation=45)

                        plt.plot(data[column], **self.kwargs)
                        if self.annotate_column is not None:
                            try:
                                for vline in self.annotate_data[label][self.annotate_column]:
                                    plt.axvline(vline,color='black',linewidth=3,alpha=0.5)
                            except KeyError:
                                pass

                        elif self.annotate_data is not None:
                            try:
                                for vline in self.annotate_data[label]:
                                    plt.axvline(vline,color='black',linewidth=3,alpha=0.5)
                            except KeyError:
                                pass
                                
                        

            # Tight layout usually dispalys nicer
            plt.tight_layout()

            # If run_id is > -1, display run number on figure
            if(obj_data.run_id > -1):
                figure.suptitle( "Run: " + str(obj_data.run_id), y=1.02)
