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

import numpy as np
import pandas as pd
from sklearn.tree import DecisionTreeRegressor

class OffsetDetrend(PipelineItem):
    '''
    Trend filter that fits a stepwise function to linearly detrended table data

    On detrended data this filter fits a stepwise function (number of
    steps provided by the user) to correct the linear fit by
    accounting for discontinuous offsets, such as due to a change in
    the antenna or from an earthquake. The final linear fit handles
    each portion of the offset independently. If the number of
    discontinuities is not provided as an autoparam, the filter
    assumes a single discontinuity.
    '''

    def __init__(self, str_description, column_names, ap_paramList = [], labels=None, time_point=None, time_interval=None):
        ''' 
        Initialize OffsetDetrend filter for use on table data
        
        @param str_description: String describing filter
        @param column_names: List of column names to select data to be removed (using None will apply to all columns)
        @param ap_paramList[step_count]: Number of steps to remove from data (Default: 1)
        @param labels: List of labels used to select data to be removed (using None will apply to all labels)
        @param time_point: Time of offset
        @param time_interval: Interval within which the offset occurs
        '''
        
        self.labels = labels
        self.column_names = column_names
        self.time_point = time_point
        if time_interval == None:
            self.time_interval = [-500,500]
        else:
            if type(time_interval) == int:
                self.time_interval = [-time_interval,time_interval]        
            else:
                self.time_interval = time_interval
        self.ap_paramNames = ['step_count']
            
        super(OffsetDetrend, self).__init__(str_description, ap_paramList)

    
    def process(self, obj_data):
        ''' 
        Apply offset estimation and detrending filter to data set.
        
        @param obj_data: Input data. Changes are made in place.
        '''

        labels = self.labels
        column_names = self.column_names

        # user provided number of steps/offsets in the data
        step_count = 1
        if len(self.ap_paramList) != 0:
            step_count = self.ap_paramList[0]()
            
        for label, data in obj_data.getIterator():
            for column in column_names:
                if (labels is None or label in labels):
                    # keep track of the time index and the location of nan's                
                    tindex = data.index
                    reind = np.array(np.isnan(data))
                    # a temporary time index and data array without nan's
                    nts = np.arange(len(data))
                    nts = np.delete(nts,nts[reind])
                    nys = data[reind==False]

                    # Decision Tree Regressor for finding the discontinuities
                    regr_1 = DecisionTreeRegressor(max_depth=step_count)
                    if self.time_point == None:
                        regr_1.fit(nts[:,np.newaxis], nys)
                    else:
                        # make time_point (a string) into an index
                        time_point = np.where(tindex==self.time_point)[0][0]
                        regr_1.fit(nts[(time_point+self.time_interval[0]):(time_point+self.time_interval[1]),np.newaxis],
                                   nys[(time_point+self.time_interval[0]):(time_point+self.time_interval[1])])
                    r1 = regr_1.predict(nts[:,np.newaxis])

                    # offset the discontinuity to be continous and fit a single line
                    # (using median of 5 points on either side of discontinuity)
                    nys[r1==r1[-1]] += np.median(nys[r1==r1[0]][-5:-1]) - np.median(nys[r1==r1[-1]][0:5])
                    z3 = np.polyfit(nts, nys, 1)

                    # make the data into a pd series and correctly index
                    x3 = pd.Series(data=nys-(z3[0]*nts+z3[1]),index=tindex[reind==False])
                    x3 = x3.reindex(tindex)
                    # and then use that to update in place
                    obj_data.updateData(label, x3.index, column, x3)
            
