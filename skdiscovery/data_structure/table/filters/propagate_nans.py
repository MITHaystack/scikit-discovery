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

# Library imports
import pandas as pd
import numpy as np

# skdiscovery imports
from skdiscovery.data_structure.framework.base import PipelineItem


class PropagateNaNs(PipelineItem):
    ''' Propagates NaN's from one column to other columns '''
    
    def __init__(self,str_description,nan_column,target_columns):
        '''
        Initialize PropagateNaNs Filter

        @param str_description: String describing filter
        @param nan_column: Column used to select which rows should be NaN's
        @param target_columns: Rows in these column will be set to NaN's based on nan_column
        '''

        super(PropagateNaNs,self).__init__(str_description,[])
        self.nan_column = nan_column
        self.target_columns = target_columns
        
    def process(self, obj_data):
        '''
        PropagateNaNs on table data wrapper

        @param obj_data: Input table data wrapper
        '''
        
        for label, data in obj_data.getIterator():
            nan_index = data[pd.isnull(data[self.nan_column])].index
            
            for column in self.target_columns:
                new_data = data[column].copy()
                new_data[nan_index] = np.nan
                
                obj_data.updateData(label,new_data.index,column,new_data)
