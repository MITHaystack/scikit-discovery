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

from skdiscovery.data_structure.framework.base import PipelineItem
import pandas as pd


class Resample(PipelineItem):
    '''
    Resample data
    '''
    
    def __init__(self, str_description, start_date = None, end_date = None, frequency='D'):
        '''
        Initialize Resample filter

        @param str_description: String describing filter
        @param start_date: Starting date
        @param end_date: Ending date
        @param period: New sampling rate
        @param frequency: Frequency of the resmampled data
                          (see Pandas DataFrame reindex for more options)
        '''
        self.start_date = start_date
        self.end_date = end_date
        self.frequency = frequency
        super(Resample, self).__init__(str_description)
        
    def process(self, obj_data):
        '''
        Calibrates GRACE, updating in place

        @param obj_data: Table data wrapper
        '''
        label_list = []
        frame_list = []
        for label, data in obj_data.getIterator():

            if self.start_date == None:
                start_date = data.index[0]
            else:
                start_date = self.start_date

            if self.end_date == None:
                end_date = data.index[-1]
            else:
                end_date = self.end_date

            # Bug in Data frame occasionally prevents reindexing data frames will all nans
            if (~pd.isnull(data)).sum().sum() != 0:
                frame_list.append(data.reindex(pd.date_range(start_date, end_date, freq=self.frequency)))
            else:
                frame_list.append(pd.DataFrame(columns = data.columns, index = pd.date_range(start_date, end_date, freq=self.frequency)))

            label_list.append(label)
            
        obj_data.updateFrames(label_list,frame_list)
