# The MIT License (MIT)
# Copyright (c) 2018 Massachusetts Institute of Technology
#
# Authors: Cody Rude
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

# 3rd party imports
import numpy as np
import pandas as pd


def getPCAComponents(pca_results):
    '''
    Retrieve PCA components from PCA results

    @param pca_results: PCA results from a pipeline run

    @return Pandas DataFrame containing the pca components
    '''

    date_range = pd.date_range(pca_results['start_date'],  pca_results['end_date'])
    column_names = ['PC' + str(i+1) for i in range(pca_results['CA'].n_components)]
    pca = pd.DataFrame(data = pca_results['Projection'], index = date_range, columns=column_names)
    pca.index.name='Date'

    return pca

def rotate(col_vectors, az, ay, ax):
    '''
    Rotate col vectors in three dimensions

    Rx * Ry * Rz * row_vectors

    @param col_vectors: Three dimensional Column vectors
    @param az: Z angle
    @param ay: Y angle
    @param ax: X angle

    @return rotated col vectors
    '''
    rz = np.array([[np.cos(az), -np.sin(az), 0], [np.sin(az), np.cos(az), 0], [0, 0, 1]])
    ry = np.array([[np.cos(ay), 0, np.sin(ay)], [0, 1, 0], [-np.sin(ay), 0, np.cos(ay)]])
    rx = np.array([[ 1, 0, 0], [0, np.cos(ax), -np.sin(ax)], [0, np.sin(ax), np.cos(ax)]])

    rot = rx @ ry @ rz

    return rot @ col_vector

def translate(col_vectors, delta_x, delta_y, delta_z):
    '''
    Translate col vectors by x, y, and z

    @param col_vectors: Row vectors of positions

    @param delta_x: Amount to translate in the x direction
    @param delta_y: Amount to translate in the y direction
    @param delta_z: Amount to translate in the y direction
    '''

    col_vectors = col_vectors.copy()
    col_vectors[0,:] += delta_x
    col_vectors[1,:] += delta_y
    col_vectors[2,:] += delta_z

    return col_vectors


def formatColorbarLabels(colorbar, pad=29):
    """
    Adjust the labels on a colorbar so they are right aligned

    @param colorbar: Input matplotlib colorbar
    @param pad: Amount of padding to use
    """
    for t in colorbar.ax.get_yticklabels():
        t.set_horizontalalignment('right')
        colorbar.ax.yaxis.set_tick_params(pad=pad)
