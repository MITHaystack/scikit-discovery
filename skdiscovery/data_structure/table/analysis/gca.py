# The MIT License (MIT)
# Copyright (c) 2015 Massachusetts Institute of Technology
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
from sklearn.decomposition import PCA
from sklearn.decomposition import FastICA


class General_Component_Analysis(PipelineItem):
    ''' 
    Performs a general component analysis on table data.

    Currently, the two built-in types of analysis are either ICA or PCA.
    '''

    def __init__(self, str_description, ap_paramList, n_components, column_names, **kwargs):
        '''
        Initialize Analysis object

        @param str_description: String description of analysis
        @param ap_paramList[component_type]: Type of CA; either PCA or ICA
        @param ap_paramList[start_time]: Starting time for CA
        @param ap_paramList[end_time]: ending time for CA
        @param n_components: Number of components to compute
        @param column_names: Columns names to use 
        @param kwargs: Extra keyword arguments to pass on to ICA (ignored for PCA)
        '''
        
        self.str_description = str_description
        self.ap_paramList = ap_paramList
        self.ap_paramNames = ['component_type','start_time','end_time']
        self.n_components = n_components
        self.column_names = column_names
        self.kwargs = kwargs
        self.results = dict()

    
    def process(self, obj_data):
        ''' 
        Perform component analysis on data

        Results are added to the data wrapper as a dictionary with
        results['CA'] = Eigenvenctors
        results['Projection'] = Projection on to the eigenvectors
        
        @param obj_data: Data wrapper
        '''

        component_type = self.ap_paramList[0]()
        start_time = self.ap_paramList[1]()
        end_time = self.ap_paramList[2]()

        num_components = self.n_components

        results = dict()
        results['start_date'] = start_time
        results['end_date'] = end_time

        cut_data = []
        label_list = []
        for label, data  in obj_data.getIterator():
            for column in self.column_names:
                cut_data.append(data.loc[start_time:end_time, column])
                label_list.append(label)

        cut_data = np.array(cut_data)

        if len(cut_data) > 0:
            if component_type == 'ICA' :
                ca = FastICA(n_components = num_components, **self.kwargs)
            else:
                ca = PCA(n_components = num_components)
            
            time_projection = ca.fit_transform(cut_data.T)
            results['CA'] = ca
            results['Projection'] = time_projection

        else:
            results['CA'] = None
            results['Projection'] = None

        results['labels'] = label_list

        obj_data.addResult(self.str_description, results)
