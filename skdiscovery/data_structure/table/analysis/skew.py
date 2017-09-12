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


from collections import defaultdict

from skdiscovery.data_structure.framework import PipelineItem
from scipy.stats import skew
import numpy as np

class Skew(PipelineItem):
    ''' Calculates the skew of table data '''
    def process(self, obj_data):
        '''
        Apply Skew analysis with results added to the data wrapper

        @param obj_data: Data wrapper
        '''

        column_names = obj_data.getDefaultColumns()
        
        results = defaultdict(dict)
        # for label, frame in tqdm(obj_data.getIterator()):
        for label, frame in obj_data.getIterator():
            for column in column_names:
                # dropping missing data in order to remove top and bottom 2%
                data = frame[column].dropna()
                # Remove top and bottom 2%
                rem_num = round(len(data)*0.02)
                res = skew(data.sort_values(ascending=True)[rem_num:-rem_num])
                if isinstance(res, np.ma.masked_array):
                    res = np.float(res.data)
                results[label][column] = res
                obj_data.addResult(self.str_description, results)
