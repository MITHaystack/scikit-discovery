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


class TrendFilter(PipelineItem):
    '''
    Trend Filter that removes linear and sinusoidal (annual, semi-annual) trends on series data.
    '''
    def __init__(self, str_description, ap_paramList):
        ''' 
        Initialize Trend Filter

        @param str_description: String describing filter
        @param ap_paramList[list_trendTypes]: List of trend types. List can contain 
                           any mix of "linear", "annual", or "semiannual". The default
                           is to remove the linear, annual, and semiannual trends
        '''
        super(TrendFilter, self).__init__(str_description, ap_paramList)
        self.ap_paramNames = ['trend_list']
        
    def process(self, obj_data):
        ''' 
        Apply trend filter to data set.
        
        @param obj_data: Input data. Changes are made in place.
        '''

        filter_list = None
        if len(self.ap_paramList) != 0:
            filter_list = self.ap_paramList[0].val()
            
        for label, data, err in obj_data.getIterator():
            if filter_list == None or 'linear' in filter_list:
                data.update(trend_tools.getTrend(data)[0])

            if filter_list == None or 'semiannual' in filter_list:
                data.update(trend_tools.sinuFits(data))
            elif 'annual' in filter_list:
                data.update(trend_tools.sinuFits(data, fitN=1))
