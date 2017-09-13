# The MIT License (MIT)
# Copyright (c) 2017 Massachusetts Institute of Technology
#
# Authors: Victor Pankratius, Justin Li, Cody Rude
# This software has been created in projects supported by the US National
# Science Foundation and NASA (PI: Pankratius)
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

from skdiscovery.data_structure.framework.base import PipelineItem
import numpy as np
import pandas as pd
from statsmodels.robust import mad


class MIDAS(PipelineItem):
    '''
    *In Development* A basic MIDAS trend estimator

    See http://onlinelibrary.wiley.com/doi/10.1002/2015JB012552/full
    '''
        
    def __init__(self, str_description,column_names = None):
        '''
        Initiatlize the MIDAS filtering item

        @param str_description: String description of filter
        @param column_names: List of column names to analyze
        '''
        
        super(MIDAS, self).__init__(str_description, [])
        self.column_names = column_names

    def process(self, obj_data):
        '''
        Apply the MIDAS estimator to generate velocity estimates

        Adds the result to the data wrapper

        @param obj_data: Data wrapper
        '''
        
        if self.column_names == None:
            column_names = obj_data.getDefaultColumns()
        else:
            column_names = self.column_names

        time_diff = pd.to_timedelta('365d')
        results = dict()
        for label, data in obj_data.getIterator():
            start_date = data.index[0]
            end_date = data.index[-1]
            for column in column_names:
                start_data = data.loc[start_date:(end_date-time_diff), column]
                end_data = data.loc[start_date+time_diff:end_date, column]
                
                offsets = end_data.values - start_data.values
                offsets = offsets[~np.isnan(offsets)]
                med_off = np.median(offsets)
                mad_off = mad(offsets)
                
                cut_offsets = offsets[np.logical_and(offsets < med_off + 2*mad_off, 
                                                     offsets > med_off - 2*mad_off)]
                final_vel = np.median(cut_offsets)
                final_unc = np.sqrt(np.pi/2) * mad(cut_offsets) / np.sqrt(len(cut_offsets))

                results[label] = pd.DataFrame([final_vel,final_unc], ['velocity', 'uncertainty'] ,[column])
                
        obj_data.addResult(self.str_description, pd.Panel.fromDict(results,orient='minor'))
        

