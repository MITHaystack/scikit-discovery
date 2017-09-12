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
import scipy as sp
import matplotlib.pyplot as plt


class HCluster(PipelineItem):
    '''
    Hierarchical Clustering function that produces a cluster map of the distance matrix
    '''
    def __init__(self, str_description, obj_name):
        '''
        Initialize HCluster

        @param str_description: String describing accumulator
        @param obj_name: Name of distance matrix parameter in the obj_data results
        '''
        
        super(HCluster, self).__init__(str_description)
        self.obj_name = obj_name
        
    def process(self, obj_data):
        '''
        Produces a cluster map and stores the linkage results.
    
        @param obj_data: Data wrapper
        '''

        import seaborn as sns
        
        data = obj_data.getResults()[self.obj_name]

        linkage = sp.cluster.hierarchy.linkage(data, method='average')

        plt.figure()
        
        g = sns.clustermap(data, col_linkage = linkage, row_linkage=linkage)
        for item in g.ax_heatmap.get_yticklabels():
            item.set_rotation(0)
        

        plt.figure()

        sp.cluster.hierarchy.dendrogram(linkage)

        obj_data.addResult(self.str_description, linkage)
