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

# 3rd part imports
import numpy as np
import pandas as pd
from scipy.optimize import brute
from fastdtw import fastdtw

# scikit discovery imports
from skdiscovery.data_structure.framework import PipelineItem
from skdiscovery.utilities.patterns import trend_tools as tt

# Standard library imports
from collections import OrderedDict


class RotatePCA(PipelineItem):
    """
    *** In Development *** Class for rotating PCA to seperate superimposed signals
    """

    def __init__(self, str_description, ap_paramList, pca_name, model, norm=None, num_components=3):
        '''

        @param str_description: String description of this item
        @param ap_paramList[fit_type]: Fitness test to use (either 'dtw' or 'remove')
        @param ap_paramList[resolution]: Fitting resolution when using brute force
        @param pca_name: Name of pca results
        @param model: Model to compare to (used in dtw)
        @param norm: Normalization to use when comparing data and model (if None, absolute differences are used)
        @param num_components: Number of pca components to use
        '''
        self._pca_name = pca_name
        self._model = tt.normalize(model)
        self.norm = norm

        if num_components not in (3,4):
            raise NotImplementedError('Only 3 or 4 components implemented')

        self.num_components = num_components
        super(RotatePCA, self).__init__(str_description, ap_paramList)



    def _rotate(self, col_vector, az, ay, ax):
        '''
        Rotate column vectors in three dimensions

        Rx * Ry * Rz * col_vectors

        @param col_vector: Data as a column vector
        @param az: Z angle
        @param ay: Y angle
        @param ax: X angle

        @return rotated column vectors
        '''

        rz = np.array([[np.cos(az), -np.sin(az), 0], [np.sin(az), np.cos(az), 0], [0, 0, 1]])
        ry = np.array([[np.cos(ay), 0, np.sin(ay)], [0, 1, 0], [-np.sin(ay), 0, np.cos(ay)]])
        rx = np.array([[ 1, 0, 0], [0, np.cos(ax), -np.sin(ax)], [0, np.sin(ax), np.cos(ax)]])

        rot = rx @ ry @ rz

        return rot @ col_vector

    def _rotate4d(self, col_vector, rot_angles):
        '''
        Rotate column vectors in four dimensions

        @param col_vector: Data as a column vector
        @param rot_angles: Rotation angles ('xy', 'yz', 'zx', 'xw', 'yw', 'zw')
        @return rotated column vectors
        '''

        index_list = []
        index_list.append([0,1])
        index_list.append([1,2])
        index_list.append([0,2])
        index_list.append([0,3])
        index_list.append([1,3])
        index_list.append([2,3])

        # Two different types:
        # left sine is negative: Type 0
        # right sine is negative: Type 1
        type_list = [0, 0, 1, 0, 1, 1]

        rotation_dict = OrderedDict()

        # The order of the rotation matrix is as follows:
        # (see https://hollasch.github.io/ray4/Four-Space_Visualization_of_4D_Objects.html#s2.2)
        label_list = ['xy', 'yz', 'zx', 'xw', 'yw', 'zw']
        for angle, label, index, negative_type in zip(rot_angles, label_list, index_list, type_list):
            ct = np.cos(angle)
            st = np.sin(angle)

            rotation_matrix = np.eye(4)
            rotation_matrix[index[0], index[0]] = ct
            rotation_matrix[index[1], index[1]] = ct
            rotation_matrix[index[0], index[1]] = st
            rotation_matrix[index[1], index[0]] = st

            if negative_type == 0:
                rotation_matrix[index[1], index[0]] *= -1
            elif negative_type == 1:
                rotation_matrix[index[0], index[1]] *= -1
            else:
                raise RuntimeError('Invalid value of negative_type')

            rotation_dict[label]=rotation_matrix


        rot_matrix = np.eye(4)
        for label, matrix in rotation_dict.items():
            rot_matrix = rot_matrix @ matrix

        return rot_matrix @ col_vector


    def _rollFastDTW(self, data, centered_tiled_model, model_size):
        '''
        Compute minimum fastdtw distance for a model to match length of real data at all possible phases


        @param data: Real input data
        @param centered_tiled_model: Model after being tiled to appropriate length and normalized (mean removed an scaled by standard devation)
        @param model_size: Size of the original model (before tiling)
        @return Index of minimum distance, minimum distance
        '''

        centered_data = tt.normalize(data)

        fitness_values = [fastdtw(centered_data, np.roll(centered_tiled_model, i), dist=self.norm)[0] for i in range(model_size)]
        min_index = np.argmin(fitness_values)

        return min_index, fitness_values[min_index]


    def _tileModel(self, in_model, new_size):
        '''
        Tile a model to increase its length

        @param in_model: Input model
        @param new_size: Size of tiled model
        @return Tiled model
        '''

        num_models = int(np.ceil(new_size / len(in_model)))
        return np.tile(in_model, num_models)[:new_size]

    def _fitness(self, z, data, model, fit_type = 'dtw', num_components=3):
        '''
        Compute fitness of data given a model and rotation

        @param z: Rotation angles
        @param data: Input data
        @param model: Input model
        @param fit_type: Choose fitness computation between dynamic time warping ('dtw') or 
                         by comparing to an seasonal and linear signal ('remove')
        @param num_components: Number of pca components to use. Can be 3 or 4 for fit_type='dtw' 
                               or 3 for fit_type='remove'
        @return fitness value

        '''
        if num_components == 3:
            new_data = self._rotate(data.as_matrix().T, *z)
        elif num_components == 4:
            new_data = self._rotate4d(data.as_matrix().T, z)


        if fit_type == 'dtw':
            return self._fitnessDTW(new_data, model, num_components)
        elif fit_type == 'remove' and num_components == 3:
            return self._fitnessRemove(pd.DataFrame(new_data.T, columns=['PC1','PC2','PC3'],
                                                    index=data.index))
        elif fit_type == 'remove':
            raise NotImplementedError("The 'remove' fitness type only works with 3 components")
        else:
            raise NotImplementedError('Only "dtw" and "remove" fitness types implemented')


    def _fitnessDTW(self, new_data, model, num_components=3):
        '''
        Compute fitness value using dynamic time warping

        @param new_data: Input data
        @param model: Input model
        @param: Number of pca components to use (3 or 4)

        @return fitness value using dynamic time warping
        '''


        tiled_model = tt.normalize(self._tileModel(model, new_data.shape[1]))

        roll, primary_results = self._rollFastDTW(new_data[num_components-1,:], tiled_model, len(model))

        # pc1_results = np.min([fastdtw(tt.normalize(new_data[0,:]), np.roll(tiled_model, roll))[0],
        #                       fastdtw(-tt.normalize(-new_data[0,:]), np.roll(tiled_model, roll))[0]])

        # pc2_results = np.min([fastdtw(tt.normalize(new_data[1,:]), np.roll(tiled_model, roll))[0],
        #                       fastdtw(tt.normalize(-new_data[1,:]), np.roll(tiled_model, roll))[0]])

        other_pc_results = 0
        for i in range(num_components-1):
            other_pc_results += self._rollFastDTW(new_data[i,:], tiled_model, len(model))[1]

        return primary_results - other_pc_results

    def _fitnessRemove(self, new_data):
        '''
        fitness value determined by how well seasonal and linear signals can be removed frm first two components

        @param new_data: Input data
        @return fitness value determined by comparison of first two components to seasonal and linear signals
        '''

        linear_removed = tt.getTrend(new_data['PC1'].asfreq('D'))[0]
        annual_removed = tt.sinuFits(new_data['PC2'].asfreq('D'), 1, 1)

        return linear_removed.var() + annual_removed.var()


    def process(self, obj_data):
        '''
        Compute rotation angles for PCA

        @param obj_data: Input table data wrapper
        '''

        fit_type = self.ap_paramList[0]()
        resolution = self.ap_paramList[1]()

        pca_results = obj_data.getResults()[self._pca_name]

        date_range = pd.date_range(pca_results['start_date'],  pca_results['end_date'])
        column_names = ['PC' + str(i+1) for i in range(pca_results['CA'].n_components)]
        pca = pd.DataFrame(data = pca_results['Projection'], index = date_range, columns=column_names)
        pca.index.name='Date'

        pca = pca.loc[:,['PC' + str(i+1) for i in range(self.num_components)]]


        end_point = 360 - (360/resolution)

        if self.num_components == 3:
            num_ranges = 3
        elif self.num_components == 4:
            num_ranges = 4
        else:
            raise ValueError('Wrong number of components')

        ranges = []
        for i in range(num_ranges):
            ranges.append((0, np.deg2rad(end_point)))

        new_angles = brute(func=self._fitness,
                           ranges=ranges,
                           Ns=resolution,
                           args=(pca, self._model, fit_type, self.num_components))


        final_score = self._fitness(new_angles, pca, self._model, fit_type, self.num_components)

        rotated_pcs = pd.DataFrame(self._rotate(pca.T, *new_angles).T, index=pca.index, columns = pca.columns)

        results = OrderedDict()

        results['rotation_angles'] = new_angles
        results['rotated_pcs'] = rotated_pcs
        results['final_score'] = final_score
        results['rotated_components'] =  self._rotate(pca_results['CA'].components_, *new_angles)

        obj_data.addResult(self.str_description, results)
