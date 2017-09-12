# The MIT License (MIT)
# Copyright (c) 2017 Massachusetts Institute of Technology
#
# Authors: Cody Rude
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


from skdiscovery.data_structure.framework.base import PipelineItem
from skdiscovery.data_structure.framework.stagecontainers import *
from skdaccess.framework.param_class import *
from skdaccess.geo.imsdnhs import DataFetcher as SDF
from collections import OrderedDict


class SnowFusion(PipelineItem):
    ''' 
    Adds snow time series data to table based on geographic coordinates 

    Works on table data (original data from http://nsidc.org/data/g02156)
    '''
    def __init__(self, str_description, metadata, column_data_name = 'Snow'):
        '''
        Initialize Snow Fusion item

        @param str_description: String describing item
        @param metadata: Metadata that contains lat,lon coordinates based on data labels
        @param column_data_name: Name of column for Snow data
        '''
        
        super(SnowFusion, self).__init__(str_description, [])
        self.metadata = metadata
        self.column_data_name = column_data_name



    def process(self, obj_data): 
        ''' 
        Adds column for snow (g02156) data

        @param obj_data: Input DataWrapper, will be modified in place
        '''
        coordinate_dict = OrderedDict()
        for label, data in obj_data.getIterator():
            coordinate_dict[label] = (self.metadata[label]['Lat'], self.metadata[label]['Lon'])
            start_date = data.index[0]
            end_date = data.index[-1]
            sdf = SDF(coordinate_dict, start_date, end_date)
            snow = sdf.output().get()[label].loc[:,'Snow']

            obj_data.addColumn(label, self.column_data_name, snow)
