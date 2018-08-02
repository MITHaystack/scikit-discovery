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
import numpy as np


class WeightedAverage(PipelineItem):
    ''' This filter performs a rolling weighted average using standard deviations as weight '''

    def __init__(self, str_description, ap_paramList, column_names, std_dev_column_names=None, propagate_uncertainties=False):
        '''
        Initializes a WeightedAverage object

        @param str_description: String describing filter
        @param ap_paramList[window]: Window to use for computing rolling weighted average
        @param column_names: Names of columns to apply the weighted average
        @param std_dev_column_names: Names of columns of the standard deviations. If none a regular mean is computed.
        @param propagate_uncertainties: Propagate uncertainties assuming uncorrelated errors
        '''

        super(WeightedAverage,self).__init__(str_description, ap_paramList)
        self.column_names = column_names
        self.std_dev_column_names = std_dev_column_names
        self.propagate_uncertainties = propagate_uncertainties


    def process(self, obj_data):
        '''
        Apply the moving (weighted) average filter to a table data wrapper.n

        Changes are made in place.

        @param obj_data: Input table data wrapper
        '''

        window = self.ap_paramList[0]()

        for label, data in obj_data.getIterator():

            if self.std_dev_column_names != None:
                for column, std_dev_column in zip(self.column_names,
                                                  self.std_dev_column_names):

                    weights = 1 / data[std_dev_column]**2
                    weighted_data = data[column] * weights

                    scale = weights.rolling(window=window,center=True, min_periods=1).sum()

                    weighted_average = weighted_data.rolling(window=window, center=True, min_periods=1).sum() / scale

                    obj_data.updateData(label, weighted_average.index, column,weighted_average)

                    if self.propagate_uncertainties:
                        # Uncertainty determined using the standard error propagation technique
                        # Assumes data is uncorrelated
                        uncertainty =  1 / np.sqrt(scale)
                        obj_data.updateData(label, uncertainty.index, std_dev_column, uncertainty)


            else:
                for column in self.column_names:

                    weighted_average = data[column].rolling(window=window, center=True, min_periods=1).mean()

                    obj_data.updateData(label,weighted_average.index,column,weighted_average)
