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
    Computes the correlation for table data and stores the result as a matrix.
    '''

    def __init__(self, str_description, column_names = None, local_match = False, correlation_type = 'pearson'):
        '''
        Initialize Correlate analysis item for use on tables

        @param str_description: String describing analysis item
        @param column_names: List of column names to correlate
        @param local_match: Only correlate data on the same frames
        @param correlation_type: Type of correlation to be passed to pandas ('pearson', 'kendall', 'spearman')
        '''
        super(Correlate, self).__init__(str_description,[])
        self.column_names = column_names
        self.local_match = local_match
        self.corr_type = correlation_type
    
    
    def process(self, obj_data):
        ''' 
        Computes the correlation between columns and stores the results in obj_data

        @param obj_data: Data wrapper
        '''

        if self.column_names == None:
            column_names = obj_data.getDefaultColumns()
        else:
            column_names = self.column_names
        
        if self.local_match == False:
            data = []
            index = []

            for label, data_in in obj_data.getIterator():
                for column in column_names:
                    data.append(data_in[column])
                    index.append(label + '.' + column)

            index = np.array(index)


            result = []
            for s1 in data:
                row = []
                for s2 in data:
                    row.append(s1.corr(s2, method=self.corr_type))
                result.append(row)

            obj_data.addResult(self.str_description, pd.DataFrame(result, index=index, columns=index))

        else:
            full_results = dict()
            for label, data_in in obj_data.getIterator():
                data = []
                index = []
                for column in column_names:
                    data.append(data_in[column])
                    index.append(column)

                result = []
                for s1 in data:
                    row = []
                    for s2 in data:
                        row.append(s1.corr(s2, method=self.corr_type))
                    result.append(row)

                full_results[label] = pd.DataFrame(result,index=index,columns=index)

            obj_data.addResult(self.str_description, pd.Panel.from_dict(full_results))
