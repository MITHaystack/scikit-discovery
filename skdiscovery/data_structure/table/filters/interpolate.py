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
import pandas as pd
import numpy as np


class InterpolateFilter(PipelineItem):
    ''' Interpolate missing values on table data '''
    
    def process(self, obj_data):
        ''' 
        Interpolate missing data in obj_data DataWrapper 

        @param obj_data: Input DataWrapper, will be modified in place
        '''
        
        for label, data in obj_data.getIterator():
            for column in obj_data.getDefaultColumns():
                index = pd.Series(np.arange(len(data)))
                good_index = index.loc[np.array(pd.notnull(data[column]))]
                nan_index = index.loc[np.array(pd.isnull(data[column]))]
                interpolated_values = np.interp(nan_index, good_index, data.iloc[good_index].loc[:,column])

                obj_data.updateData(label, data.index[nan_index], column, interpolated_values)
                


