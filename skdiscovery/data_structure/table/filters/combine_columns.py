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



# Framework import
from skdiscovery.data_structure.framework.base import PipelineItem

# 3rd party libraries import
import pandas as pd


class CombineColumns(PipelineItem):
    '''
    Create a new column by selecting data from a column

    Fills in any missing values using a second column
    '''
    def __init__(self, str_description, column_1, column_2, new_column_name):
        '''
        Initialize a CombineColumns object

        @param str_description: String describing filter
        @param column_1: Name of primary column
        @param column_2: Name of secondary column to be used
        when data from the primary column is not avaiable
        @param new_column_name: Name of resulting column
        '''
        self.column_1 = column_1
        self.column_2 = column_2
        self.new_column_name = new_column_name
        
        super(CombineColumns,self).__init__(str_description)
    def process(self, obj_data):
        ''' 
        Apply combine column filter to data set, operating on the data_obj
        
        @param obj_data: Table data wrapper.
        '''
        for label, data in obj_data.getIterator():
            if self.column_1 in data.columns and self.column_2 in data.columns:
                 # replacing all null median data with mean data
                col1_null_index = pd.isnull(data.loc[:,self.column_1])


                data.loc[:,self.new_column_name] = data.loc[:,self.column_1]

                # Check if there is any replacement data available
                if (~pd.isnull(data.loc[col1_null_index, self.column_2])).sum() > 0:
                    data.loc[col1_null_index, self.new_column_name] = data.loc[col1_null_index, self.column_2]
                
                
            elif self.column_2 in data.columns and self.column_1 not in data.columns:
                data.loc[:,self.new_column_name] = data.loc[:,self.column_2]

            elif self.column_2 not in data.columns and self.column_1 in data.columns:
                data.loc[:,self.new_column_name] = data.loc[:,self.column_1]

            else:
                raise KeyError('data needs either "' + self.column_2 + '" or "' + self.column_1 + '" or both')
