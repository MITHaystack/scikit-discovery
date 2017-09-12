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
from skdiscovery.utilities.patterns import pbo_tools
from scipy.optimize import curve_fit
from collections import defaultdict
import pandas as pd
import numpy as np


class HTanFilter(PipelineItem):
    '''
    Filter to subtract arctan fit from data

    [DEPRECATED] [will be removed]
    '''
    def __init__(self, str_description, t0, amplitude=5, timescale=1., offset=0,
                 slope=0, labels=None, column_names=None, start_time_limit = None,
                 end_time_limit = None, start = None, end = None):
        '''
        Fit and remove hyperbolic tangent function from data.

        @param str_description: String description of data
        @param t0: Initial time offset of arctangent
        @param amplitude: Initial amplitude of arctangent
        @param timescale: Timescale of fit
        @param offset: Initial Y offset of arctangent
        @param slope: Slope of the data 
        @param labels: Labels to apply arctangent function to
        @param column_names: Column names to apply arctanget function to
        @param start_time_limit: Starting time bound for fit to arctan (default: no bound)
        @param end_time_limit: Ending time bound for fit to arctan (default: no bound)
        @param start: Index of the first data point to fit (default: index of first data point)
        @param end: Index of the last data point to fit (default: index of last data point)
        '''
        self.a = amplitude
        self.t0 = t0
        self.c = timescale
        self.slope = slope
        self.offset = offset
        self.labels = labels
        self.column_names = column_names
        self.start_time_limit = start_time_limit
        self.end_time_limit = end_time_limit
        self.start = start
        self.end = end

        super(HTanFilter, self).__init__(str_description, [])


    def process(self, obj_data):
        ''' 
        Apply Arctangent filter to data param.

        @param obj_data: Input data. Changes are made in place.
        '''
        
        parameters = [self.a, pbo_tools.datetimeToNumber(self.t0), self.c, self.offset, self.slope]

        def fitfunc(t,a,t0,c,b,m):

            # return np.piecewise(t, [t <= t0, t > t0], [ lambda t: a*np.arctan((t-t0)/c)+b1+m1*t,
            #                                             lambda t: a*np.arctan((t-t0)/c)+b2+m2*t ] )

            return a*np.tanh((t-t0)/c)+b+m*t

        labels = self.labels
        column_names = self.column_names

        results = defaultdict(pd.DataFrame)


        if self.start_time_limit is not None:
            start_time_limit = pbo_tools.datetimeToNumber(self.start_time_limit)
        else:
            start_time_limit = -np.inf

        if self.end_time_limit is not None:
            end_time_limit = pbo_tools.datetimeToNumber(self.end_time_limit)
        else:
            end_time_limit = np.inf


        bounds = ( (-np.inf, start_time_limit, -np.inf, -np.inf, -np.inf ),
                   ( np.inf, end_time_limit,    np.inf,  np.inf,  np.inf) )



        for label, data, err in obj_data.getIterator():

            if (labels is None or label in labels) and \
               (column_names is None or data.name in column_names):


                index = data.index

                if self.start is None:
                    start = index[0]
                else:
                    start = self.start

                if self.end is None:
                    end = index[-1]
                else:
                    end = self.end

                # new_index = pbo_tools.datetimeToNumber(data.index)
                new_index = pbo_tools.datetimeToNumber(data.loc[start:end].index)
                full_new_index = pbo_tools.datetimeToNumber(data.index)
                
                
                pA, success = curve_fit(fitfunc, new_index, data.loc[start:end], p0=parameters, sigma = err.loc[start:end],
                                        absolute_sigma=True, bounds=bounds)

                data.update(data - fitfunc(full_new_index, *pA))

                results[label][data.name] = pd.Series(fitfunc(full_new_index, *pA), index=index)

        obj_data.addResult(self.str_description, pd.Panel(results))
