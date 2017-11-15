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


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy as sp

from skdaccess.utilities import pbo_util
from skdiscovery.data_structure.series.analysis.mogi import MogiVectors

from skdiscovery.utilities.patterns import pbo_tools


def calc_distance_map(pipeline, ap_name, ca_name, ca_type ,plotFlag = True, histIdx = False, fontsize=10):
    """ 
    Calculates distances/similarities between pipeline runs

    Optionally visualizes the result as a seaborn clustermap for PBO
    pipelines (requires multiple stations)

    Calculates the square root of the summed squared differences between eigenvectors.
    Only works, because of internal assumptions, on pipelines with multiple stations
    Returns the distances as a pandas dataframe
    
    @param pipeline: Pipeline to analyze.
    @param ap_name: Name of the pipeline item that is being perturbed
    @param ca_name: Name of the pipeline item used as the comparison metric for calculating the distance
    @param ca_type: Type of comparison metric [PCA for PCA, MogiSource of Mogi Source, MogiVector for Mogi vectors]
    @param plotFlag: Boolean flag for plotting the clustermap of distances
    @param histIdx: Flag for returning the perturbed pipeline item parameters
    @param fontsize: Fontsize adjustments

    @return cg: The generated clustermap of the calculated distances/similarities
    @return dist_mat: A matrix of the calculated distances/similarities
    @return history: The record of the perturbed pipeline item parameters
    """
    # a history of the perturbed pipeline item
    history = []
    for runInfo in pipeline.getMetadataHistory():
        for stageItem in runInfo:
            if ap_name in stageItem:
                history.append(stageItem.rsplit(ap_name)[1].strip('[]'))
        
    # number of runs
    num_results = len(pipeline.RA_results)  
    dist_mat = np.zeros([num_results, num_results])
        
    # compute distances between all pairs of runs
    for i in range(num_results):
        for j in range(i, num_results):       
            # Check the ca_name type to properly format the type of comparison
            if ca_type == 'PCA':
                ctitle='PCA Vector Distance Similarity'
                summation = 0
                rstation_list = pipeline.RA_results[i][ca_name]['labels']
                num_stations = len(rstation_list)
                # some redundancy because of the way the modules/functions were structured
                # eigenvectors (lat and lon) for the one run
                coord_list = pbo_util.getStationCoords(pipeline.data_fetcher.meta_data, rstation_list)
                _,_,EV1_lat,EV1_lon,_ = pbo_tools.dirEigenvectors(coord_list, pipeline.RA_results[i][ca_name]['CA'].components_[0])
                # and for the other run
                _,_,EV2_lat,EV2_lon,_ = pbo_tools.dirEigenvectors(coord_list, pipeline.RA_results[j][ca_name]['CA'].components_[0])
                
                for k in range(num_stations):
                    ev1 = np.hstack((EV1_lat[k], EV1_lon[k]))
                    ev2 = np.hstack((EV2_lat[k], EV2_lon[k]))
                    # calculate the euclidean distance difference at each station
                    summation += sp.spatial.distance.euclidean(ev1, ev2)**2
                
            elif ca_type == 'MogiSource':
                ctitle='Mogi Source Distance [deg+km]'
                # for now, just uses lat, lon, and depth for Mogi comparison
                ev1 = np.array([pipeline.RA_results[i][ca_name]['lat'], pipeline.RA_results[i][ca_name]['lon'], pipeline.RA_results[i][ca_name]['depth']])
                ev2 = np.array([pipeline.RA_results[j][ca_name]['lat'], pipeline.RA_results[j][ca_name]['lon'], pipeline.RA_results[j][ca_name]['depth']])
                summation = sp.spatial.distance.euclidean(ev1, ev2)**2
            elif ca_type == 'MogiVector':
                ctitle='Summed Mogi Vector Distance [mm]'
                # do the same comparison as eigenvectors for the mogi modeled vectors
                summation = 0
                rstation_list = pipeline.RA_results[i][ca_name]['labels']
                num_stations = len(rstation_list)
                # got mogi vectors for the two runs
                coord_list = np.array(pbo_util.getStationCoords(pipeline.data_fetcher.meta_data, rstation_list))
                mogi_x_1, mogi_y_1 = MogiVectors(pipeline.RA_results[i][ca_name],coord_list[:,0],coord_list[:,1])
                mogi_x_2, mogi_y_2 = MogiVectors(pipeline.RA_results[j][ca_name],coord_list[:,0],coord_list[:,1])
                
                for k in range(num_stations):
                    ev1 = np.hstack((mogi_x_1[k], mogi_y_1[k]))
                    ev2 = np.hstack((mogi_x_2[k], mogi_y_2[k]))
                    # calculate the euclidean distance difference at each station
                    summation += sp.spatial.distance.euclidean(ev1, ev2)**2
                # as all pca amplitudes are the same, to scale Mogi to difference in mm
                summation *= (pipeline.RA_results[0][ca_name]['pca_amplitude'])**2
                
            dist_mat[i][j] = np.sqrt(summation)

    dist_mat += dist_mat.transpose()
    if histIdx:
        dist_mat = pd.DataFrame(dist_mat, index=['Configuration '+str(ii).zfill(2) for ii in range(len(history))], columns=['Configuration '+str(ii).zfill(2) for ii in range(len(history))])
    else:
        dist_mat = pd.DataFrame(dist_mat, index=history, columns=history)
    
    if plotFlag:
        # Importing seaborn here in order to
        # return matplotlib to it's original style
        import seaborn as sns
        sns.set(font='serif')

        cg = sns.clustermap(dist_mat)
        plt.setp(cg.ax_heatmap.yaxis.get_majorticklabels(), rotation=0, fontsize=fontsize);
        plt.setp(cg.ax_heatmap.xaxis.get_majorticklabels(), rotation=90, fontsize=fontsize);
        cg.cax.set_title(ctitle,fontsize=fontsize);
        sns.reset_orig()

        if histIdx:
            return cg, history
        else:
            return cg
    else:
        if histIdx:
            return dist_mat, history
        else:
            return dist_mat
