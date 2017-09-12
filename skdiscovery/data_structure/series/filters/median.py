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
    A Median filter for series data
    '''

    def __init__(self, str_description, ap_paramList, interpolate=True, subtract = False):
        '''
        Initialize MedianFilter
        
        @param str_description: String describing filter
        @param ap_paramList[ap_window]: median filter window width
        @param interpolate: Flag to interpolate data points before filtering
        @param subtract: Flag to subtract filtered result from original
        '''

        self.interpolate = interpolate
        self.subtract = subtract
        self.ap_paramNames = ['windowSize']
        super(MedianFilter, self).__init__(str_description, ap_paramList)

    
    def process(self, obj_data):
        ''' 
        Apply median filter to data set
        
        @param obj_data: Input DataWrapper. Changes are made in place.
        '''

        ap_window = self.ap_paramList[0]()

        for label, data, err in obj_data.getIterator():

            result = trend_tools.medianFilter(data, ap_window, self.interpolate)
            if self.subtract == True:
                data.update(data - result)
            else:
                data.update(result)
