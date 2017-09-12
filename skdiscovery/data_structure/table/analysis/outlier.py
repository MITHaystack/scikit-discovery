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
from statsmodels.robust import mad


class Outlier(PipelineItem):

    ''' 
    Computes (data / mad(data)) for outlier detection

    Creates a new column for the result
    '''

    def __init__(self, str_description, columns=None, name_prefix = 'MAD_Scale_'):
        '''
        Initalize Outlier Item

        @param str_description: Name of Item
        @param columns: List of of column names
        @param name_prefix: Prefix of newly created column
        '''

        self.columns = columns
        self.name_prefix = name_prefix
        super(Outlier, self).__init__(str_description)

    def process(self, obj_data):
        '''
        Process the data object to add a column with the outlier scores

        @param obj_data: Input table data wrapper
        '''

        if self.columns == None:
                column_names = obj_data.getDefaultColumns()
        else:
                column_names = self.columns

        for label, data in obj_data.getIterator():
            for column in column_names:

                clean_data = data[column].dropna()
                
                mad_scale = clean_data / mad(clean_data)
                obj_data.addColumn(label, self.name_prefix + column, mad_scale)


