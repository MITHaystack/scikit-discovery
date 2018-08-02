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


class CalibrateGRACEMascon(PipelineItem):
    '''
    Calibrate Grace Data

    This can apply a scale factor and round dates to the nearest day
    '''

    def __init__(self, str_description, round_dates=True, apply_scale_factor = True):
        '''
        Initialize GRACE Mascon calibration filter

        @param str_description: String describing filter
        @param round_dates: Option for rounding to dates to the nearest day
        @param apply_scale_factor: Boolean indicating whether or not a corrective scale factor should be applied
        '''
        self.round_dates = round_dates
        self.apply_scale_factor = apply_scale_factor
        super(CalibrateGRACEMascon, self).__init__(str_description)

    def process(self, obj_data):
        '''
        Calibrates GRACE, updating in place

        @param obj_data: Table data wrapper
        '''
        label_list = []
        frame_list = []
        for label, data in obj_data.getIterator():
            label_list.append(label)

            new_data = data
            if self.apply_scale_factor == True:
                new_data.loc[:,'EWD'] = data.loc[:,'EWD'] * obj_data.info(label)['scale_factor']

            if self.round_dates == True:
                new_data.index = new_data.index.round('D')

            frame_list.append(new_data)

        obj_data.updateFrames(label_list,frame_list)
