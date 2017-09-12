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

from skdiscovery.utilities.patterns import trend_tools


class MedianFilter(PipelineItem):
    '''
    A Median filter for table data
    '''

    def __init__(self, str_description, ap_paramList, interpolate=True,
                 subtract = False,regular_period=True, min_periods=1):
        '''
        Initialize MedianFilter

        @param str_description: String describing filter
        @param ap_paramList[ap_window]: median filter window width
        @param interpolate: Interpolate data points before filtering
        @param subtract: Subtract filtered result from original
        @param regular_period: Assume the data is regularly sampled
        @param min_periods: Minimum required number of data points in window
        '''

        self.interpolate = interpolate
        self.subtract = subtract
        self.ap_paramNames = ['windowSize']
        self.regular_period = regular_period
        self.min_periods = min_periods
        super(MedianFilter, self).__init__(str_description, ap_paramList)

    
    def process(self, obj_data):
        ''' 
        Apply median filter to data set
        
        @param obj_data: Input panda's data series. Changes are made in place.
        '''

        ap_window = self.ap_paramList[0]()

        column_names = obj_data.getDefaultColumns()

        for label, data in obj_data.getIterator():
            for column in column_names:
                if self.interpolate == True or self.regular_period == False:
                    result = trend_tools.medianFilter(data[column], ap_window, self.interpolate)
                else:
                    result = data[column].rolling(ap_window,min_periods=self.min_periods, center=True).median()
                    
                if self.subtract == True:
                    obj_data.updateData(label, data.index, column, data[column] - result)
                else:
                    obj_data.updateData(label, data.index, column, result)
                    
