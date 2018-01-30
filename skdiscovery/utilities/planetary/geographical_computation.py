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

################################################################################
# Function definitions
################################################################################

# Distance computation

@jit(nopython=True)
def haversine_distance_math(longitude_1, latitude_1, longitude_2, latitude_2, radius):
    rad_longitude_1 = longitude_1*math.pi/180.
    rad_latitude_1 = latitude_1*math.pi/180.
    rad_longitude_2 = longitude_2*math.pi/180.
    rad_latitude_2 = latitude_2*math.pi/180.
    
    diff_longitudes = rad_longitude_2 - rad_longitude_1
    diff_latitudes = rad_latitude_2 - rad_latitude_1
    
    root = math.sin(diff_latitudes/2.)**2 + math.cos(rad_latitude_1)*math.cos(rad_latitude_2)*math.sin(diff_longitudes/2.)**2
    return 2*radius*math.asin(math.sqrt(root))

# Computation of geographical measurements

def nvector_from_lonlat(longitude_1, latitude_1):
    rad_longitude_1 = longitude_1*np.pi/180.
    rad_latitude_1 = latitude_1*np.pi/180.
    return np.array([np.cos(rad_longitude_1)*np.cos(rad_latitude_1),
                     np.sin(rad_longitude_1)*np.cos(rad_latitude_1), 
                     np.sin(rad_latitude_1)])

def compute_great_circle_nvector(nvector_1,
                                 bearing, distance, 
                                 planet_radius):
    
    k_east = np.cross([0, 0, 1], nvector_1)
    magnitude_k_east = np.sqrt(k_east[0]**2 + k_east[1]**2 + k_east[2]**2)
    k_east /= magnitude_k_east
    k_north = np.cross(nvector_1, k_east)
    
    d_x = k_north[0]*np.cos(bearing) + k_east[0]*np.sin(bearing)
    d_y = k_north[1]*np.cos(bearing) + k_east[1]*np.sin(bearing)
    d_z = k_north[2]*np.cos(bearing) + k_east[2]*np.sin(bearing)
    
    nvector_2_x = nvector_1[0]*np.cos(distance/planet_radius) + d_x*np.sin(distance/planet_radius)
    nvector_2_y = nvector_1[1]*np.cos(distance/planet_radius) + d_y*np.sin(distance/planet_radius)
    nvector_2_z = nvector_1[2]*np.cos(distance/planet_radius) + d_z*np.sin(distance/planet_radius)
    
    return [nvector_2_x, nvector_2_y, nvector_2_z]

def lonlat_from_nvector(nvector_1):
    rad_longitude_1 = np.arctan2(nvector_1[1], nvector_1[0])
    rad_latitude_1 = np.arctan2(nvector_1[2], np.sqrt(nvector_1[0]**2 + nvector_1[1]**2))
    return rad_longitude_1*180/np.pi, rad_latitude_1*180/np.pi

def mod(y, x):
    return y - x*np.floor(y/x)

def compute_great_circle_distance_and_bearing(rad_longitude_1, rad_latitude_1, 
                                              rad_longitude_2, rad_latitude_2, 
                                              planet_radius):
    
    delta_longitude = rad_longitude_1 - rad_longitude_2
    delta_latitude = rad_latitude_1 - rad_latitude_2
    
    # Haversine distance
    root = np.sin(delta_latitude/2.)**2 + np.cos(rad_latitude_1)*np.cos(rad_latitude_2)*np.sin(delta_longitude/2.)**2
    distance = 2*planet_radius*np.arcsin(np.sqrt(root))
    
    bearing = np.arctan2(np.sin(delta_longitude)*np.cos(rad_latitude_2), 
                         np.cos(rad_latitude_1)*np.sin(rad_latitude_2) - np.sin(rad_latitude_1)*np.cos(rad_latitude_2)*np.cos(delta_longitude))
    bearing = mod(-bearing, 2*np.pi)
    
    return distance, bearing

@jit(nopython = True)
def nvector_from_lonlat_math(rad_longitude_1, rad_latitude_1):
    nvector = [0., 0., 0.]
    nvector[0] = math.cos(rad_longitude_1)*math.cos(rad_latitude_1)
    nvector[1] = math.sin(rad_longitude_1)*math.cos(rad_latitude_1)
    nvector[2] = math.sin(rad_latitude_1)
    return nvector

@jit(nopython = True)
def cross(vector_a, vector_b):
    vector = [0., 0., 0.]
    vector[0] = vector_a[1]*vector_b[2] - vector_a[2]*vector_b[1]
    vector[1] = vector_a[2]*vector_b[0] - vector_a[0]*vector_b[2]
    vector[2] = vector_a[0]*vector_b[1] - vector_a[1]*vector_b[0]
    return vector

@jit(nopython = True)
def scalar_division(vector_a, scalar):
    vector = [0., 0., 0.]
    vector[0] = vector_a[0]/scalar
    vector[1] = vector_a[1]/scalar
    vector[2] = vector_a[2]/scalar
    return vector

@jit(nopython = True)
def compute_great_circle_nvector_math(nvector_1,
                                      bearing, distance, 
                                      planet_radius):
    
    k_east = cross([0., 0., 1.], nvector_1)
    magnitude_k_east = math.sqrt(k_east[0]**2 + k_east[1]**2 + k_east[2]**2)
    k_east = scalar_division(k_east, magnitude_k_east)
    k_north = cross(nvector_1, k_east)
    
    d_x = k_north[0]*math.cos(bearing) + k_east[0]*math.sin(bearing)
    d_y = k_north[1]*math.cos(bearing) + k_east[1]*math.sin(bearing)
    d_z = k_north[2]*math.cos(bearing) + k_east[2]*math.sin(bearing)
    
    nvector_2 = [0., 0., 0.]
    nvector_2[0] = nvector_1[0]*math.cos(distance/planet_radius) + d_x*math.sin(distance/planet_radius)
    nvector_2[1] = nvector_1[1]*math.cos(distance/planet_radius) + d_y*math.sin(distance/planet_radius)
    nvector_2[2] = nvector_1[2]*math.cos(distance/planet_radius) + d_z*math.sin(distance/planet_radius)
    
    return nvector_2

@jit(nopython = True)
def lonlat_from_nvector_math(nvector_1):
    rad_longitude_1 = math.atan2(nvector_1[1], nvector_1[0])
    rad_latitude_1 = math.atan2(nvector_1[2], math.sqrt(nvector_1[0]**2 + nvector_1[1]**2))
    return rad_longitude_1*180/math.pi, rad_latitude_1*180/math.pi

@jit(nopython = True)
def mod_math(y, x):
    return y - x*math.floor(y/x)

@jit(nopython = True)
def compute_great_circle_distance_and_bearing_math(rad_longitude_1, rad_latitude_1, 
                                                   rad_longitude_2, rad_latitude_2, 
                                                   planet_radius):
    
    delta_longitude = rad_longitude_1 - rad_longitude_2
    delta_latitude = rad_latitude_1 - rad_latitude_2

    root = math.sin(delta_latitude/2.)**2 + math.cos(rad_latitude_1)*math.cos(rad_latitude_2)*math.sin(delta_longitude/2.)**2
    haversine_distance = 2*planet_radius*math.asin(math.sqrt(root))
    
    bearing = math.atan2(math.sin(delta_longitude)*math.cos(rad_latitude_2), 
                         math.cos(rad_latitude_1)*math.sin(rad_latitude_2) - math.sin(rad_latitude_1)*math.cos(rad_latitude_2)*math.cos(delta_longitude))
    bearing = mod_math(-bearing, 2*math.pi)
    
    return haversine_distance, bearing

# Computation of geographical maps

def compute_longitude_and_latitude_maps(lon_min, lon_max, 
                                        lat_min, lat_max, 
                                        raster_width, raster_height):
    pixel_lon_size = (lon_max - lon_min)/raster_width
    pixel_lat_size = (lat_max - lat_min)/raster_height
    longitudes = np.linspace(lon_min + 0.5*pixel_lon_size, 
                             lon_max - 0.5*pixel_lon_size, 
                             raster_width)
    # GDAL array starts from the max value of latitude, keep that here for coherency
    latitudes = np.linspace(lat_max - 0.5*pixel_lat_size, 
                            lat_min + 0.5*pixel_lat_size, 
                            raster_height)
    return np.meshgrid(longitudes, latitudes)

def compute_surface_area(raster_longitude_array, 
                         raster_latitude_array, 
                         lon_min, lon_max, 
                         lat_min, lat_max,
                         planet_radius):
    pixel_lon_size = (lon_max - lon_min)/raster_longitude_array.shape[1]
    pixel_lat_size = (lat_max - lat_min)/raster_longitude_array.shape[0]

    north_raster_latitude_array = (raster_latitude_array + 0.5*pixel_lat_size)*np.pi/180.
    south_raster_latitude_array = (raster_latitude_array - 0.5*pixel_lat_size)*np.pi/180.

    east_raster_longitude_array = (raster_longitude_array + 0.5*pixel_lon_size)*np.pi/180.
    west_raster_longitude_array = (raster_longitude_array - 0.5*pixel_lon_size)*np.pi/180.

    return planet_radius*planet_radius*(east_raster_longitude_array - 
                                        west_raster_longitude_array)*(np.sin(north_raster_latitude_array) - 
                                                                      np.sin(south_raster_latitude_array))