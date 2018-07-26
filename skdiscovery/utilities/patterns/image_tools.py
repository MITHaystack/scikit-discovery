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


import statsmodels.api as sm
import numpy as np
import imreg_dft as ird
import shapely
import scipy as sp

def buildMatchedPoints(in_matches, query_kp, train_kp):
    '''
    Get postions of matched points

    @param in_matches: Input matches
    @param query_kp: Query key points
    @param train_kp: Training key points

    @return Tuple containing the matched query and training positions
    '''
    query_index = [match.queryIdx for match in in_matches]
    train_index = [match.trainIdx for match in in_matches]

    sorted_query_kp = [query_kp[i] for i in query_index]
    sorted_train_kp = [train_kp[i] for i in train_index]


    query_positions = [[kp.pt[0], kp.pt[1]] for kp in sorted_query_kp]
    train_positions = [[kp.pt[0], kp.pt[1]] for kp in sorted_train_kp]

    return query_positions, train_positions

def scaleImage(input_data, vmin=None, vmax=None):
    '''
    Scale image values to be within 0 and 255

    @param input_data: Input data
    @param vmin: Minimum value for scaled data, where smaller values are clipped, defaults to Median - stddev as determined by mad
    @param vmax: Maximum value for scaled data, where larger values are clipped, defaults to Median - stddev as determined by mad)
    @return input_data scaled to be within 0 and 255 as an 8 bit integer
    '''
    if vmin==None or vmax==None:
        stddev = sm.robust.mad(input_data.ravel())
        middle = np.median(input_data.ravel())

    if vmin == None:
        vmin = middle - 1*stddev

    if vmax == None:
        vmax = middle + 1*stddev

    input_data = input_data.astype(np.float)
    input_data[input_data<vmin] = vmin
    input_data[input_data>vmax] = vmax

    input_data = np.round((input_data - vmin) * 255 / (vmax-vmin)).astype(np.uint8)

    return input_data


def divideIntoSquares(image, size, stride):
    """
    Create many patches from an image

    Will drop any patches that contain NaN's

    @param image: Source image
    @param size: Size of one side of the square patch
    @param stride: Spacing between patches (must be an integer greater than 0)

    @return Array containing the extent [x_start, x_end, y_start, y_end] of each patch and an array of the patches
    """
    def compute_len(size, stride):
        return (size-1) // stride + 1


    num_x = compute_len(image.shape[-1]-size, stride)
    num_y = compute_len(image.shape[-2]-size, stride)

    if image.ndim == 2:
        array_data = np.zeros((num_x * num_y, size, size), dtype = image.dtype)

    elif image.ndim == 3:
        array_data = np.zeros((num_x * num_y, image.shape[0], size, size), dtype = image.dtype)


    extent_data = np.zeros((num_x * num_y, 4), dtype = np.int)

    index = 0

    for x in range(0, image.shape[-1]-size, stride):
        for y in range(0, image.shape[-2]-size, stride):
            if image.ndim == 2:
                cut_box = image[y:y+size, x:x+size]
            elif image.ndim == 3:
                cut_box = image[:, y:y+size, x:x+size]

            array_data[index, ...] = cut_box
            extent_data[index, :] = np.array([x, x+size, y, y+size])
            index += 1

    if image.ndim==2:
        valid_index = ~np.any(np.isnan(array_data), axis=(1,2))
    else:
        valid_index = ~np.any(np.isnan(array_data), axis=(1,2,3))

    return extent_data[valid_index], array_data[valid_index]


def generateSquaresAroundPoly(poly, size=100, stride=20):
    '''
    Generate that may touch a shapely polygon

    @param poly: Shapely polygon
    @param size: Size of boxes to create
    @param stride: Distance between squares

    @return list of Shapely squares that may touch input polygon
    '''
    x_start, x_end = np.min(poly.bounds[0]-size).astype(np.int), np.max(poly.bounds[2]+size).astype(np.int)
    y_start, y_end = np.min(poly.bounds[1]-size).astype(np.int), np.max(poly.bounds[3]+size).astype(np.int)

    x_coords = np.arange(x_start, x_end+1, stride)
    y_coords = np.arange(y_start, y_end+1, stride)

    x_mesh, y_mesh = np.meshgrid(x_coords, y_coords)

    return [shapely.geometry.box(x, y, x+size, y+size) for x, y in zip(x_mesh.ravel(), y_mesh.ravel())]
