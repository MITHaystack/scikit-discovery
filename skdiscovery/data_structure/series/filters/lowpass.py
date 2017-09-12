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
from scipy import signal


class LowPassFilter(PipelineItem):
    '''
    A FIR Remez (Parks-McLellan) designed low pass filter for series data
    '''
    
    def __init__(self, str_description, ap_paramList):
        '''
        Initialize LowPassFilter
        @param str_description: String describing filter
        @param ap_paramList[ntaps]: Number of filter taps
        @param ap_paramList[fpassf_per]: Frequency passband ratio/percentage
        @param ap_paramList[fstopf_per]: Frequency stopband ratio/percentage
        @param ap_paramList[wghts]: Band importance weights
        @param ap_paramList[miter]: Maximum number of iterations for generating the filter
        '''

        super(LowPassFilter, self).__init__(str_description, ap_paramList)
        self.ap_paramNames = ['ntaps','fpassf_per','fstopf_per','weights','miter']
        
    def process(self, obj_data):
        ''' 
        Apply lowpass filter to data set, with changes applied in place
        
        @param obj_data: Input data with data
        ''' 
        ntaps = self.ap_paramList[0]()
        fpassf_per = self.ap_paramList[1]()
        fstopf_per = self.ap_paramList[2]()
        wghts = self.ap_paramList[3]()
        miter = self.ap_paramList[4]()
        b_filt=signal.fir_filter_design.remez(numtaps=ntaps,
                                              bands=[0,fpassf_per,fstopf_per,.5],
                                              desired=[1,0],weight=wghts,maxiter=miter)
        
        for label, data, err in obj_data.getIterator():
            data.update(signal.filtfilt(b_filt,1,data))
                    
