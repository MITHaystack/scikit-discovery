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

# scikit discovery imports
from skdiscovery.data_structure.framework import PipelineItem
from skdiscovery.utilities.patterns.trend_tools import normalize

# Standard library imports
from collections import OrderedDict


class RotatePCA(PipelineItem):

    def __init__(self, str_description, pca_name, model, resolution = 20):

        self._pca_name = pca_name
        self._model = normalize(model)
        self._resolution = resolution
        super(RotatePCA, self).__init__(str_description)



    def _rotate(self, row_vector, az, ay, ax):
        rz = np.array([[np.cos(az), -np.sin(az), 0], [np.sin(az), np.cos(az), 0], [0, 0, 1]])
        ry = np.array([[np.cos(ay), 0, np.sin(ay)], [0, 1, 0], [-np.sin(ay), 0, np.cos(ay)]])
        rx = np.array([[ 1, 0, 0], [0, np.cos(ax), -np.sin(ax)], [0, np.sin(ax), np.cos(ax)]])

        rot = rx @ ry @ rz

        return rot @ row_vector


    def _shape_fitness(self, z, data):
        az, ay, ax = z
        new_data = self._rotate(data.as_matrix().T, az, ay, ax)

        normalized_pc1 = normalize(new_data[0,:])
        normalized_pc2 = normalize(new_data[1,:])
        normalized_pc3 = normalize(new_data[2,:])

        #linear_removed = tt.getTrend(new_data['PC1'].asfreq('D'))[0]
        #annual_removed = tt.sinuFits(new_data['PC2'].asfreq('D'), 1, 1)

        return np.std(np.correlate(normalized_pc1, self._model),ddof=1) + np.std(np.correlate(normalized_pc2, self._model),ddof=1) - \
               np.std(np.correlate(normalized_pc3, self._model),ddof=1)


    def process(self, obj_data):
        pca_results = obj_data.getResults()[self._pca_name]

        date_range = pd.date_range(pca_results['start_date'],  pca_results['end_date'])
        column_names = ['PC' + str(i+1) for i in range(pca_results['CA'].n_components)]
        pca = pd.DataFrame(data = pca_results['Projection'], index = date_range, columns=column_names)
        pca.index.name='Date'


        end_point = 360 - (360/self._resolution)
        new_angles = brute(self._shape_fitness,
                           ranges = ((0, np.deg2rad(end_point)),
                                     (0, np.deg2rad(end_point)),
                                     (0, np.deg2rad(end_point))),
                           Ns=self._resolution,
                           args=(pca,))


        final_score = self._shape_fitness(new_angles, pca)


        rotated_pcs = pd.DataFrame(self._rotate(pca.T, *new_angles).T, index=pca.index, columns = pca.columns)

        results = OrderedDict()

        results['rotation_angles'] = new_angles
        results['rotated_pcs'] = rotated_pcs
        results['final_score'] = final_score

        obj_data.addResult(self.str_description, results)
