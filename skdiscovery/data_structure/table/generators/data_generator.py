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

from skdaccess.framework.data_class import DataFetcherBase
from skdaccess.framework.data_class import TableWrapper
import pandas as pd
import numpy as np

class DataGenerator(DataFetcherBase):

    ''' In Development: Class for generating random data '''
    
    def __init__(self, length, *args, seed = None, final_function = None):
        '''
        Initialize Random data generator


        @param length: Number of rows to generate
        @param *args: Dictionaries containing entries: 'name', 'start', 'end', and optionally 'func'
        @param seed: Seed to use when generating random data
        @param final_function: Final function to call on random data
        '''

        self.length = length
        self.seed = seed
        self.args = args
        self.final_function = final_function

    def output(self):
        '''
        Generate data

        @return Table data wrapper of generated data
        '''
        if self.seed is not None:
            np.random.seed(self.seed)
            

        new_data = dict()
        name_list = []
        for arg in self.args:
            new_data[arg['name']] = np.random.rand(self.length) * (arg['end'] - arg['start']) + arg['start']

            name_list.append(arg['name'])

            if 'func' in arg:
                new_data[arg['name']] = arg['func'](new_data[arg['name']])


        if self.final_function is not None:
            new_data = pd.DataFrame.from_dict(new_data)
            new_data, updated_column_names = self.final_function(new_data)
            if updated_column_names is not None:
                default_columns = updated_column_names


        data = {'generated_data' : new_data}
        data_wrapper = TableWrapper(data, default_columns = name_list)
        
        return data_wrapper
