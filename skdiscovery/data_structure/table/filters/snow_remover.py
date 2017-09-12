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
from skdaccess.framework.param_class import AutoParam
import numpy as np

class SnowRemover(PipelineItem):
    ''' Removes data with snow errors '''

    def __init__(self,str_description,ap_paramList = [AutoParam(1.5)],column_name='dN', snow_column='Snow'):
        '''
        Initialize snow remover for use on table data.

        @param str_description: String describing filter
        @param ap_paramList[sigma_clip]: remove station if the stddev of snowdays is sigma_clip 
                                         times greater than non-snow days, default 1.5
        @param column_name: Name of column to check
        @param snow_column: Name of snow column to determine snowdays/non snow days
        '''
        self.column_name = column_name
        self.snow_column = snow_column


        super(SnowRemover, self).__init__(str_description, ap_paramList)


    def process(self, obj_data):
        ''' 
        Removes table data with large snow errors

        @param obj_data: Input DataWrapper, will be modified in place
        '''

        bad_stations = []
        sigma_multiplier = self.ap_paramList[0]()

        for label, data in obj_data.getIterator():

            if len(data[data[self.snow_column]==4]) > 0 and len(data[data[self.snow_column]==2]) > 0:
                snow = data[data[self.snow_column]==4].loc[:,self.column_name]
                no_snow = data[data[self.snow_column]==2].loc[:,self.column_name]

                non_snow_std = np.nanstd(no_snow)
                snow_std = np.nanstd(snow)

                if snow_std > sigma_multiplier * non_snow_std:
                    bad_stations.append(label)

        if len(bad_stations) > 0:
            obj_data.removeFrames(bad_stations)

