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
from skdaccess.utilities.pbo_util import removeAntennaOffset


class AntennaOffset(PipelineItem):
    '''
    Applies corrections to fix offsets in PBO GPS data induced by antenna changes
    '''

    def __init__(self, str_description, antenna_data, min_diff = 0.0, column_list = None):
        '''
        Initialize AntennaOffset function

        @param str_description: String describing the filter
        @param antenna_data: Data containing the log of antenna changes
        @param min_diff: Difference in position needed to be considered an offset
        @param column_list: Names of the columns to apply the function to
        '''

        super(AntennaOffset, self).__init__(str_description, [])
        self.antenna_data = antenna_data
        self.column_list = column_list
        self.min_diff = min_diff


    def process(self, obj_data):
        '''
        Applies the function to the data, updating in place

        @param obj_data: Table data wrapper
        '''
        
        if self.column_list == None:
            column_names = obj_data.getDefaultColumns()

        else:
            column_names = self.column_list


        for label, data in obj_data.getIterator():
            for column in column_names:
                try:
                    updated_data = removeAntennaOffset(self.antenna_data[label], data[column], min_diff = self.min_diff)
                    index = data[column].index
                    obj_data.updateData(label,index,column,updated_data)
                    
                except KeyError:
                    pass
                    
