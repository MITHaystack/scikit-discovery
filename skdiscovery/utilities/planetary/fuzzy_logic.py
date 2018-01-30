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

import numpy as np

################################################################################
# Function definitions
################################################################################

# Membership functions

def trapezoidal_function(raster_array, 
                        x_start_rise, x_start_plateau, x_end_plateau, x_end_slope, 
                        bottom_value = 0.2, plateau_value = 1, nan_value = 0.1):
    
    assert x_start_rise <= x_start_plateau, 'x_start_plateau should be higher than x_start_rise'
    assert x_start_plateau <= x_end_plateau, 'x_end_plateau should be higher than x_start_plateau'
    assert x_end_plateau <= x_end_slope, 'x_end_slope should be higher than x_end_plateau'
    assert bottom_value <= plateau_value, 'bottom_value should be higher than plateau_value'
    
    filtered_raster_array = np.full(raster_array.shape, bottom_value, np.double)
    # NaN
    indices = np.where(np.isnan(raster_array))
    filtered_raster_array[indices] = nan_value
    # Left side
    if x_start_rise != x_start_plateau:
        indices = np.where(np.logical_and(x_start_rise < raster_array,
                                          raster_array < x_start_plateau))
        filtered_raster_array[indices] = \
            bottom_value + \
            (raster_array[indices] - x_start_rise)*(plateau_value - bottom_value)/ \
            float(x_start_plateau - x_start_rise)
    # Right side
    if x_end_plateau != x_end_slope:
        indices = np.where(np.logical_and(x_end_plateau < raster_array,
                                          raster_array < x_end_slope))
        filtered_raster_array[indices] = \
            plateau_value + \
            (raster_array[indices] - x_end_plateau)*(bottom_value - plateau_value)/ \
            float(x_end_slope - x_end_plateau)
    # Plateau
    indices = np.where(np.logical_and(x_start_plateau <= raster_array,
                                      raster_array <= x_end_plateau))
    filtered_raster_array[indices] = plateau_value
    
    return filtered_raster_array

# Fuzzy operators

def union(*args):
    return np.maximum.reduce(args)

def intersection(*args):
    return np.minimum.reduce(args)

def complement(raster_array_a):
    return 1 - raster_array_a

def algebraic_product(*args):
    return np.multiply.reduce(args)

def algebraic_sum(*args):
    return 1 - np.multiply.reduce([1 - arg for arg in args])

def gamma_operation(gamma, *args):
    assert gamma >= 0. and gamma <= 1., "Gamma should be between 0 and 1"
    return np.power(algebraic_sum(*args),
                    gamma)*np.power(algebraic_product(*args),
                                    1 - gamma)