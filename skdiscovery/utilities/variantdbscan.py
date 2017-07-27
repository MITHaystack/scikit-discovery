# The MIT License (MIT)
# Copyright (c) 2016 Massachusetts Institute of Technology
#
# Authors: Michael Gowanlock, Cody Rude
# This software is part of the NASA AIST Project "Computer-Aided Discovery of
# Earth Surface Deformation Phenomena" and the NSF DIBBS Project "An
# Infrastructure for Computer Aided Discovery in Geoscience", PI: V. Pankratius
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


#!/usr/bin/env python3
import os
import tempfile
import subprocess
import sys
import math
import pandas as pd
import numpy as np
from functools import total_ordering


class VariantDBScan(object):
    ''' Wrapper for VariantDBScan '''
    def __init__(self, variants, data, column_names ):
        '''
        Initialize DBScan  pipeline item

        @param variants: DataFrame of epsilon (label column 'eps') and minpoints (label column 'mp')
        @param data: Data Pandas DataFrame to be clustered 
        @param column_names: List of column names in DataFrame to cluster (Can be 2 or 3 columns)
        '''

        self.__data = data
        self.__column_names = column_names

        if type(variants) is pd.DataFrame:
            self.__variants = variants
        else:
            self.__variants = pd.DataFrame(variants,columns=('eps','mp'))

        if len(self.__column_names) == 2:
            self.__program = '/bin/jupyter/vdbscan/main'
            self.__data_file = 'output_clusters_CLUSTER_REUSE_list_'
        else:
            self.__program = '/bin/jupyter/vdbscan/VDBSCAN_3D'
            self.__data_file = 'output_clusters_MSB_list_'            

    def run(self, verbose=False):

        '''
        Runs VariantDBScan on data

        @param verbose: Print additional information about run

        @return a dataframe with a column for each dbscan run which contains the cluster id for each object.
                Note: A value of 0 indicates object is a noise point
        '''

        cwd = os.getcwd()
        
        tmp_dir = tempfile.TemporaryDirectory(dir='./')
        os.chdir(tmp_dir.name)


        temp_data = tempfile.NamedTemporaryFile(dir='./', prefix='tmp_data')
        temp_variants = tempfile.NamedTemporaryFile(dir='./', prefix='tmp_variant')

        self.__data.loc[:, self.__column_names].to_csv(temp_data.name, header=False, index=False, float_format='%.5f')
        self.__variants.to_csv(temp_variants.name, header=False, index=False)
        args=(self.__program, temp_data.name, temp_variants.name)

        popen = subprocess.Popen(args, stdout=subprocess.PIPE, shell=False)

        popen.wait()
        output = popen.stdout.read()
        if verbose:
            print(output.decode('ascii'))

        data_order = pd.read_csv('variants_after_sort.txt', skiprows=1, header=None, names=('eps','mp'))

        result_list = []
        result_names = []
        for index,row in self.__variants.iterrows():
            # Find releavent matches
            matches = np.logical_and(np.isclose(row['eps'],data_order['eps']), np.isclose(row['mp'],data_order['mp']))

            # Get the indices of the matching parametrs
            indices = np.where(matches)[0]

            # Check if there is more than one match
            if len(indices) > 1:
                print('Parameter run with repeated entries!')

            # If there are not matching data sets something has gone awry...
            elif len(indices) == 0:
                print('No matching parameter set found, aborting!')
                continue

            # The index of the matching data set:
            output_index = indices[0]
                
            results = pd.read_csv(self.__data_file + str(output_index) + '.txt', skiprows=2, names=self.__column_names +  ['cluster_id'])
            result_list.append(results['cluster_id'])
            result_names.append('cluster_' + str(index))
            

        index = pd.read_csv('data_enumerated_after_sort.txt', skiprows = 1, header=None, names= self.__column_names + ['index'])['index']

        # Cleanup
        os.chdir(cwd)
        temp_data.close()
        temp_variants.close()
        tmp_dir.cleanup()


        return pd.DataFrame(np.array(result_list).T, columns=(result_names), index=index).sort_index()
