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
from skdaccess.utilities.grace_util import computeEWD
import pandas as pd


class CalibrateGRACE(PipelineItem):
    '''
    Calibrate Grace Data

    Averages the three solutions and applies a scale factor
    '''
    
    def __init__(self, str_description, ewd_column_name='EWD', round_dates = True, apply_scale_factor = True):
        '''
        Initialize GRACE calibration filter

        @param str_description: String describing filter
        @param ewd_column_name: Name of new column for the calibrated GRACE data
        @param round_dates: Option for rounding to dates to the nearest day
        @param apply_scale_factor: Boolean indicating whether or not a corrective scale factor should be applied
        '''
        self.ewd_column_name = ewd_column_name
        self.round_dates = round_dates
        self.apply_scale_factor = apply_scale_factor
        super(CalibrateGRACE, self).__init__(str_description)
        
    def process(self, obj_data):
        '''
        Calibrates GRACE, updating in place

        @param obj_data: Table data wrapper
        '''
        label_list = []
        frame_list = []
        for label, data in obj_data.getIterator():

            if self.apply_scale_factor == True:
                scale_factor = obj_data.info(label)['scale_factor']
            else:
                scale_factor = 1

            frame_list.append(pd.DataFrame(computeEWD(data, scale_factor, round_nearest_day=self.round_dates),
                                           columns=[self.ewd_column_name]))

            frame_list[-1].loc[:,self.ewd_column_name + '_Error'] = obj_data.info(label)['measurement_error']
            label_list.append(label)
            
        obj_data.updateFrames(label_list,frame_list)
        obj_data.default_columns = [self.ewd_column_name]
