# The MIT License (MIT)
# Copyright (c) 2016 Massachusetts Institute of Technology
#
# Author: Cody
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


from skdiscovery.framework import PipelineItem
from skdiscovery.utilities import VariantDBScan


class VDBScan(PipelineItem):
    '''
    Runs Variant DBscan on table data

    Adds cluster information columns to data
    '''

    def __init__(self, str_description, variants, column_names):
        '''
        Initialize VDBScan pipelne item

        @param str_description: Description of item
        @param variants: Dataframe containing column of epsilon values and column of min points
        @param column_names: List of column names to use
        '''


        self.__column_names = column_names
        self.__variants = variants
        super(VDBScan, self).__init__(str_description, [])
    
    def process(self, obj_data):
        ''' 
        Run VDBScan on data

        @param obj_data: Data wrapper to process
        '''

        for label, data in obj_data.getIterator():
            vdb = VariantDBScan(self.__variants, data, self.__column_names)
            results = vdb.run()

            obj_data.addColumn(label, results.columns, results)


        
        
