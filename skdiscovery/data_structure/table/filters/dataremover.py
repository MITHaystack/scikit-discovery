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

class DataRemover(PipelineItem):
    ''' Sets specified table data to NaN '''    
    def __init__(self, str_description, column_names, start=None, end=None, labels=None):
        ''' 
        Initialize DataRemover
        
        @param str_description: String describing filter
        @param column_names: List of column names to select data to be removed (using None will apply to all columns)
        @param start: Starting index value
        @param end: Ending index value (inclusive)
        @param labels: List of labels used to select data to be removed (using None will apply to all labels)
        '''
        
        self.labels = labels
        self.column_names = column_names
        self.start = start
        self.end = end

        super(DataRemover, self).__init__(str_description, [])

    def process(self, obj_data):
        ''' 
        NaN's data from DataWrapper

        @param obj_data: Input DataWrapper, will be modified in place
        '''

        labels = self.labels
        
        for label, data, in obj_data.getIterator():
            if (labels is None or label in labels):
                for column in self.column_names:
                    index = data.index

                    if self.start is None:
                        start = index[0]
                    else:
                        start = self.start

                    if self.end is None:
                        end = index[-1]
                    else:
                        end = self.end

                    index=data.loc[start:end].index

                    obj_data.updateData(label, data.loc[start:end].index,column,np.nan)
