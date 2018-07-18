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

import numpy as np
import pandas as pd
from skdiscovery.data_structure.framework import PipelineItem

from skdiscovery.utilities.patterns import kalman_smoother


class KalmanFilter(PipelineItem):
    '''
    Runs a forward and backward Kalman Smoother with a FOGM state on series data

    For more information see: Ji, K. H. 2011, PhD thesis, MIT, and
    Fraser, D. C., and Potter, J. E. 1969, IEEE Trans. Automat. Contr., Acl4, 4, 387-390
    '''

    def __init__(self, str_description, ap_paramList, uncertainty_clip=5):
        '''
        Initialize Kalman Smoother

        @param str_description: String describing filter        
        @param ap_paramList[ap_tau]: the correlation time
        @param ap_paramList[ap_sigmaSq]: the data noise
        @param ap_paramList[ap_R]: the process noise
        @param uncertainty_clip: Clip data with uncertainties greater than uncertainty_clip * median uncertainty
        '''

        super(KalmanFilter, self).__init__(str_description, ap_paramList)
        self.uncertainty_clip = uncertainty_clip
        self.ap_paramNames = ['Tau','SigmaSq','R']
    
    def process(self, obj_data):
        ''' 
        Apply kalman smoother to data set
        
        @param obj_data: Input DataWrapper. Changes are made in place.
        ''' 

        uncertainty_clip = self.uncertainty_clip

        ap_tau = self.ap_paramList[0]()
        ap_sigmaSq = self.ap_paramList[1]()
        ap_R = self.ap_paramList[2]()

        for label, data, err in obj_data.getIterator():


            # if label == 'AV37' and data.name == 'dN':
            #     import matplotlib.pyplot as plt
            #     plt.figure()
            #     plt.title(label)
            #     plt.plot(err, 'o')
            #     plt.plot(data, 'o')

            # Clip data with high uncertainties
            data.loc[np.logical_and(~pd.isnull(err), err > np.nanmedian(err) * uncertainty_clip)] = np.nan

            # clip = np.nanmedian(err) * uncertainty_clip
            err.loc[np.logical_and(~pd.isnull(err), err > np.nanmedian(err) * uncertainty_clip)] = np.nan

            # If the beginning is missing data, the smoother will diverge
            if np.sum(~np.isnan(data.iloc[:20])) == 0:
                data.iloc[:2] = np.nanmedian(data)

            if ap_R == 'formal':
                R = err
            else:
                R = ap_R
            
            # Smooth the data
            smoothed, variance, t, sigma_sq, R = kalman_smoother.KalmanSmoother(data,
                                                                                t = ap_tau,
                                                                                sigma_sq = ap_sigmaSq,
                                                                                R = R)

            # Set uncertainties for missing data to those estimated from
            # the filter.
            err.loc[pd.isnull(err)] = variance[pd.isnull(err)]

            # Calculate the sample variance
            T = len(data)
            r = np.exp(-1 / t)
            sample_var = sigma_sq * (T / (T - 1)) * ( 1 - ((1+r) / (T * (1-r))) + ((2*r*(1-r**T)) / (T**2 * (1-r)**2)))
            err.update(np.sqrt(err**2 + sample_var))
            data.update(smoothed)

