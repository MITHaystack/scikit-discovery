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
    Runs a forward and backward Kalman Smoother with a FOGM state on table data

    For more information see: Ji, K. H. 2011, PhD thesis, MIT, and
    Fraser, D. C., and Potter, J. E. 1969, IEEE Trans. Automat. Contr., Acl4, 4, 387-390
    '''

    def __init__(self, str_description, ap_paramList, uncertainty_clip=5, column_names=None,
                 error_column_names = None, fillna=True):
        '''
        Initialize Kalman Smoother

        @param str_description: String describing filter        
        @param ap_paramList[ap_tau]: the correlation time
        @param ap_paramList[ap_sigmaSq]: the data noise
        @param ap_paramList[ap_R]: the process noise
        @param uncertainty_clip: Clip data with uncertainties greater than uncertainty_clip * median uncertainty
        @param column_names: List of column names to smooth (using None will apply to all columns)
        @param error_column_names: List of error column names to smooth (using None will use default error columns)
        @param fillna: Fill in missing values
        '''
        
        super(KalmanFilter, self).__init__(str_description, ap_paramList)
        self.uncertainty_clip = uncertainty_clip
        self.ap_paramNames = ['Tau','SigmaSq','R']
        self.column_names = column_names
        self.error_column_names = error_column_names
        self.fillna = fillna
    
    def process(self, obj_data):
        ''' 
        Apply kalman smoother to data set
        
        @param obj_data: Input data. Changes are made in place.
        ''' 

        uncertainty_clip = self.uncertainty_clip

        ap_tau = self.ap_paramList[0]()
        ap_sigmaSq = self.ap_paramList[1]()
        ap_R = self.ap_paramList[2]()

        if self.column_names is None:
            column_names = obj_data.getDefaultColumns()
        else:
            column_names = self.column_names

        if self.error_column_names is None:
            error_column_names = obj_data.getDefaultErrorColumns()
        else:
            error_column_names = self.error_column_names


        for label, dataframe in obj_data.getIterator():
            for column, error_column in zip(column_names, error_column_names):
                data = dataframe.loc[:,column].copy()
                err = dataframe.loc[:,error_column].copy()
                
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

                if self.fillna == True:
                    obj_data.updateData(label, data.index, error_column, np.sqrt(err**2 + sample_var))
                    obj_data.updateData(label, data.index, column, smoothed)

                else:
                    obj_data.updateData(label, dataframe.loc[:,error_column].dropna().index, error_column, np.sqrt(err**2 + sample_var))
                    obj_data.updateData(label, dataframe.loc[:,column].dropna().index, column, smoothed)

                
    def _applyKalman(self,label_dataframe,obj_data,run_Params):
        column_names = run_Params[0]
        error_column_names = run_Params[1]
        uncertainty_clip = run_Params[2]
        ap_tau = run_Params[3]
        ap_sigmaSq = run_Params[4]
        ap_R = run_Params[5]
        
        label = label_dataframe[0]
        dataframe = label_dataframe[1]
        result = {label:dict()}
        for column, error_column in zip(column_names, error_column_names):
            data = dataframe.loc[:,column].copy()
            err = dataframe.loc[:,error_column].copy()
            
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
            
            if obj_data != None:
                obj_data.updateData(label, data.index, error_column, np.sqrt(err**2 + sample_var))
                obj_data.updateData(label, data.index, column, smoothed)
                
            if obj_data == None:
                result[label]['index'] = data.index
                result[label][error_column] = np.sqrt(err**2 + sample_var)
                result[label][column] = smoothed
        if obj_data == None:
            return result, label
