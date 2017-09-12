# The MIT License (MIT)
# Copyright (c) 2017 Massachusetts Institute of Technology
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


from skdiscovery.data_structure.framework import PipelineItem
from sklearn.cluster import DBSCAN

class DBScan(PipelineItem):
    '''
    Runs DBScan on table data

    Adds cluster information column to data
    '''

    def __init__(self, str_description, ap_paramList, column_names):
        '''
        Initialize DBScan pipelne item

        @param str_description: Description of item
        @param ap_paramList[epsilon]: Distance between two nodes for them to be considered connected
        @param ap_paramList[min_points]: Minimum number of points for a cluster
        @param column_names: List of column names to use
        '''
        self.column_names = column_names
        super(DBScan, self).__init__(str_description, ap_paramList)
    
    def process(self, obj_data):
        ''' 
        Run DBScan on data. Stores result in data wrapper

        @param obj_data: Data wrapper to be processed
        '''

        epsilon = self.ap_paramList[0]()
        min_points = self.ap_paramList[1]()

        results = dict()
        

        for label, data in obj_data.getIterator():
            results[label] = DBSCAN(eps=epsilon, min_samples = min_points).fit_predict(data.loc[:,self.column_names])

        obj_data.addResult(self.str_description, results)
