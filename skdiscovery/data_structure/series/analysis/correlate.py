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


class Correlate(PipelineItem):
    '''
    Computes the correlation for series data

    Stores the result as a matrix
    '''

    def __init__(self, str_description, labels=None, column_names=None):
        '''
        Initialize Correlate analysis item

        @param str_description: String describing analysis item
        @param labels: List of labels used to select data
        @param column_names: List of column names used to select data
        '''
        
        super(Correlate, self).__init__(str_description,[])
        self.labels = labels
        self.column_names = column_names
    
    def process(self, obj_data):
        ''' 
        Computes the correlation between all the time series

        The results are stored in obj_data

        @param obj_data: Data wrapper for correlating
        '''
        data = []
        index = []

        for label, series, err in obj_data.getIterator():
            if (self.labels is None or label in self.labels) and \
               (self.column_names is None or series.name in self.column_names):

                data.append(series)
                index.append(series.name)

        index = np.array(index)
        

        result = []
        for s1 in data:
            row = []
            for s2 in data:
                row.append(s1.corr(s2))
            result.append(row)

        obj_data.addResult(self.str_description, pd.DataFrame(result, index=index, columns=index))
