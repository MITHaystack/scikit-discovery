# The MIT License (MIT)
# Copyright (c) 2018 Massachusetts Institute of Technology
#
# Author: Guillaume Rongier
# We acknowledge support from NSF ACI1442997 (PI: V. Pankratius)
#                         and NASA AISTNNX15AG84G (PI: V. Pankratius)
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

################################################################################
# Package imports
################################################################################

import math
import numpy as np

from numba import jit

from skdiscovery.utilities.planetary.geographical_computation import haversine_distance_math

################################################################################
# Function definitions
################################################################################

# Border management

@jit(nopython = True)
def add_symmetric_border(raster_array, border_size = 1):
    raster_width = raster_array.shape[1]
    raster_height = raster_array.shape[0]
    bordered_raster_array = np.full((raster_height + 2*border_size, 
                                     raster_width + 2*border_size), 
                                    np.nan)
    bordered_raster_array[border_size:-border_size, border_size:-border_size] = raster_array
    
    bordered_raster_array[border_size:-border_size, -border_size:] = raster_array[:, -1:-border_size - 1:-1]
    bordered_raster_array[border_size:-border_size, border_size::-1] = raster_array[:, 0:border_size]
    bordered_raster_array[0:border_size, :] = bordered_raster_array[2*border_size - 1:border_size - 1:-1, :]
    bordered_raster_array[-border_size:, :] = bordered_raster_array[-border_size - 1:-2*border_size - 1:-1, :]
    
    return bordered_raster_array

@jit(nopython = True)
def add_planet_border(raster_array, border_size = 1):
    raster_width = raster_array.shape[1]
    raster_height = raster_array.shape[0]
    bordered_raster_array = np.full((raster_height + 2*border_size, 
                                     raster_width + 2*border_size), 
                                    np.nan)
    bordered_raster_array[border_size:-border_size, border_size:-border_size] = raster_array
    
    bordered_raster_array[border_size:-border_size, 0:border_size] = raster_array[:, -border_size:]
    bordered_raster_array[border_size:-border_size, -border_size:] = raster_array[:, 0:border_size]
    bordered_raster_array[0:border_size, :] = bordered_raster_array[2*border_size - 1:border_size - 1:-1, ::-1]
    bordered_raster_array[-border_size:, :] = bordered_raster_array[-border_size - 1:-2*border_size - 1:-1, ::-1]
    
    return bordered_raster_array

# Computation of topographic features

@jit(nopython = True)
def compute_gradient(j, i, raster_array, longitude_array, latitude_array, planet_radius, axis = 1):
    
    cell_1 = (j + 1, i + 1)
    cell_2 = (j + 1, i - 1)
    cell_3 = (j, i + 1)
    cell_4 = (j, i - 1)
    cell_5 = (j - 1, i + 1)
    cell_6 = (j - 1, i - 1)
    if axis == 0:
        cell_2 = (j - 1, i + 1)
        cell_3 = (j + 1, i)
        cell_4 = (j - 1, i)
        cell_5 = (j + 1, i - 1)

    distance_p1 = haversine_distance_math(longitude_array[cell_1], 
                                          latitude_array[cell_1], 
                                          longitude_array[cell_2], 
                                          latitude_array[cell_2], 
                                          planet_radius)
    distance_c = haversine_distance_math(longitude_array[cell_3], 
                                         latitude_array[cell_3], 
                                         longitude_array[cell_4], 
                                         latitude_array[cell_4], 
                                         planet_radius)
    distance_m1 = haversine_distance_math(longitude_array[cell_5], 
                                          latitude_array[cell_5], 
                                          longitude_array[cell_6], 
                                          latitude_array[cell_6], 
                                          planet_radius)
    return ((raster_array[cell_1] - raster_array[cell_2])/distance_p1 +
            2*(raster_array[cell_3] - raster_array[cell_4])/distance_c +
            (raster_array[cell_5] - raster_array[cell_6])/distance_m1)/4.

@jit(nopython = True)
def compute_horne_slope(raster_array, 
                        longitude_array, 
                        latitude_array, 
                        planet_radius,
                        is_entire_planet_mapped = True):
    assert len(raster_array.shape) == 2, "Input raster is not two-dimensional"

    if is_entire_planet_mapped == True:
        raster_array = add_planet_border(raster_array)
        longitude_array = add_planet_border(longitude_array)
        latitude_array = add_planet_border(latitude_array)
    else:
        raster_array = add_symmetric_border(raster_array)
        longitude_array = add_symmetric_border(longitude_array)
        latitude_array = add_symmetric_border(latitude_array)
    
    slope_array = np.full(raster_array.shape, np.nan)

    for i in range(1, slope_array.shape[1] - 1):
        for j in range(1, slope_array.shape[0] - 1):
            dx = compute_gradient(j, i, raster_array,
                                  longitude_array, latitude_array, planet_radius)
            dy = compute_gradient(j, i, raster_array,
                                  longitude_array, latitude_array, planet_radius, 0)
            slope_array[j, i] = math.atan(math.sqrt(dx*dx + dy*dy))*180./math.pi
                
    return slope_array[1:-1, 1:-1]

@jit(nopython = True)
def compute_absolute_standard_deviation_filter(raster_array,
                                               window_size = 3,
                                               is_entire_planet_mapped = True):
    assert len(raster_array.shape) == 2, "Input raster is not two-dimensional"
    
    border_size = int(window_size/2)
    if is_entire_planet_mapped == True:
        raster_array = add_planet_border(raster_array, border_size)
    else:
        raster_array = add_symmetric_border(raster_array, border_size)

    std_array = np.full(raster_array.shape, np.nan)
    
    for i in range(border_size, std_array.shape[1] - border_size + 1):
        for j in range(border_size, std_array.shape[0] - border_size + 1):
            window = raster_array[j - border_size:j + border_size + 1, i - border_size:i + border_size + 1]
            std_array[j, i] = abs(np.nanstd(window))
            
    return std_array[border_size:-border_size, border_size:-border_size]