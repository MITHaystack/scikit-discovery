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

from matplotlib.path import Path
from matplotlib.path import get_paths_extents
from mpl_toolkits.basemap import Basemap

from numba import jit

from skdiscovery.utilities.planetary.geographical_computation import *

################################################################################
# Function definitions
################################################################################

# Build ellipse polygon

def coordinates_coding(ob):
    n = len(ob)
    codes = np.ones(n, dtype=Path.code_type)*Path.LINETO
    codes[0] = Path.MOVETO
    return codes

def create_path_from_coordinates(xy_outer_ring, xy_inner_rings = []):
    vertices = np.concatenate(
        [xy_outer_ring]
        + [r for r in xy_inner_rings])
    codes = np.concatenate(
        [coordinates_coding(xy_outer_ring)]
        + [coordinates_coding(r) for r in xy_inner_rings])
    return Path(vertices, codes)

def compute_ellipse_path(center_longitude, center_latitude, a, b, azimut, 
                         planet_radius, 
                         number_of_nodes = 100, basemap = None):
    bearings = np.linspace(0., 2*np.pi, number_of_nodes)
    radii = np.sqrt(2)*a*b/(np.sqrt((b**2 - a**2)*np.cos(2*bearings - 2*azimut) + a**2 + b**2))
    center_nvector = nvector_from_lonlat(center_longitude, center_latitude)
    nvector_1 = compute_great_circle_nvector(center_nvector,
                                             bearings, radii, 
                                             planet_radius)
    longitude_1, latitude_1 = lonlat_from_nvector(nvector_1)
    longitude_2 = []
    latitude_2 = []
    
    longitude_sign = np.sign(longitude_1)
    sz = longitude_sign == 0
    while sz.any():
        longitude_sign[sz] = np.roll(longitude_sign, 1)[sz]
        sz = longitude_sign == 0
        
    longitude_sign_change = ((np.roll(longitude_sign, 1) - longitude_sign) != 0).astype(int)
    indices = np.where(longitude_sign_change == 1)[0]

    if len(indices) == 2:
        lon_extremity_1 = np.roll(longitude_1, 1)[indices[0]]
        lat_extremity_1 = np.roll(latitude_1, 1)[indices[0]]
        lon_extremity_2 = longitude_1[indices[0]]
        lat_extremity_2 = latitude_1[indices[0]]
        
        lon_extremity_3 = np.roll(longitude_1, 1)[indices[1]]
        lat_extremity_3 = np.roll(latitude_1, 1)[indices[1]]
        lon_extremity_4 = longitude_1[indices[1]]
        lat_extremity_4 = latitude_1[indices[1]]

        if abs(lon_extremity_1 - lon_extremity_2) > 359 and \
           abs(lon_extremity_3 - lon_extremity_4) > 359 and \
           abs(lat_extremity_1 - lat_extremity_2) < 1 and \
           abs(lat_extremity_3 - lat_extremity_4) < 1:
            latitude_2 = latitude_1[longitude_sign == -1]
            longitude_2 = longitude_1[longitude_sign == -1]
            latitude_1 = latitude_1[longitude_sign == 1]
            longitude_1 = longitude_1[longitude_sign == 1]
        elif abs(lon_extremity_1 - lon_extremity_2) > 359 and \
             abs(lon_extremity_3 - lon_extremity_4) < 1 and \
             abs(lat_extremity_1 - lat_extremity_2) < 1 and \
             abs(lat_extremity_3 - lat_extremity_4) < 1:
            longitude_1 = np.insert(longitude_1, indices[0], 
                                    [np.sign(lon_extremity_1)*180, np.sign(lon_extremity_2)*180])
            latitude_1 = np.insert(latitude_1, indices[0], 
                                   [np.sign(lat_extremity_1)*90, np.sign(lat_extremity_2)*90])
        elif abs(lon_extremity_1 - lon_extremity_2) < 1 and \
             abs(lon_extremity_3 - lon_extremity_4) > 359 and \
             abs(lat_extremity_1 - lat_extremity_2) < 1 and \
             abs(lat_extremity_3 - lat_extremity_4) < 1:
            longitude_1 = np.insert(longitude_1, indices[1], 
                                    [np.sign(lon_extremity_3)*180, np.sign(lon_extremity_4)*180])
            latitude_1 = np.insert(latitude_1, indices[1], 
                                   [np.sign(lat_extremity_3)*90, np.sign(lat_extremity_4)*90])

    if basemap != None:
        longitude_1, latitude_1 = basemap(longitude_1, latitude_1)
    nodes_1 = list(zip(longitude_1, latitude_1))
    
    if len(longitude_2) != 0:
        if basemap != None:
            longitude_2, latitude_2 = basemap(longitude_2, latitude_2)
        nodes_2 = list(zip(longitude_2, latitude_2))
        return create_path_from_coordinates(nodes_1), create_path_from_coordinates(nodes_2)
    else:
        return create_path_from_coordinates(nodes_1), None

# TODO: Handle the case when x or y on a cell border, especially the max
def transform_to_pixel_coordinates(x, y, xmin, xmax, ymin, ymax, width, height):
    u = np.floor((x - xmin)*(width/(xmax - xmin)))
    v = np.floor((ymax - y)*(height/(ymax - ymin)))
    return np.int_(u), np.int_(v)

def compute_ellipse_path_bounding_box(ellipse_path,
                                      lon_min, lon_max, lat_min, lat_max, 
                                      raster_width, raster_height):
    if ellipse_path != None:
        point_0, point_1 = get_paths_extents([ellipse_path]).get_points()

        i_0, j_0 = transform_to_pixel_coordinates(point_0[0], point_0[1], 
                                                  lon_min, lon_max, lat_min, lat_max, 
                                                  raster_width, raster_height)
        if j_0 > raster_height - 1:
            j_0 = raster_height - 1
        i_1, j_1 = transform_to_pixel_coordinates(point_1[0], point_1[1], 
                                                  lon_min, lon_max, lat_min, lat_max, 
                                                  raster_width, raster_height)
        if j_1 < 0:
            j_1 = 0
        return (int(i_0), int(i_1) + 1), (int(j_1), int(j_0) + 1)
    else:
        return (None, 0), (None, 0)

def compute_ellipse_path_and_bounding_box(center_longitude, center_latitude,
                                          a, b, azimut, 
                                          lon_min, lon_max, lat_min, lat_max, 
                                          raster_width, raster_height, 
                                          planet_radius, number_of_nodes = 100):
    ellipse_path_1, ellipse_path_2 = compute_ellipse_path(center_longitude, center_latitude, a, b, azimut, 
                                                          planet_radius, 
                                                          number_of_nodes = number_of_nodes)

    slice_1_i, slice_1_j = compute_ellipse_path_bounding_box(ellipse_path_1,
                                                             lon_min, lon_max, lat_min, lat_max, 
                                                             raster_width, raster_height)
    slice_2_i, slice_2_j = compute_ellipse_path_bounding_box(ellipse_path_2,
                                                             lon_min, lon_max, lat_min, lat_max, 
                                                             raster_width, raster_height)
    
    return slice_1_i, slice_1_j, slice_2_i, slice_2_j

def compute_raster_ellipse(favorability_map_array, 
                           rad_center_longitude, rad_center_latitude, 
                           rad_longitudes, rad_latitudes,
                           planet_radius,
                           a, b, azimut,
                           ellipse_slice):
    center_distances, center_bearings = compute_great_circle_distance_and_bearing(rad_center_longitude, 
                                                                                  rad_center_latitude, 
                                                                                  rad_longitudes[ellipse_slice], 
                                                                                  rad_latitudes[ellipse_slice], 
                                                                                  planet_radius)
    ellipse_radii = np.sqrt(2)*a*b/(np.sqrt((b**2 - a**2)*np.cos(2*center_bearings - 2*azimut) + a**2 + b**2))
    center_distances[center_distances <= ellipse_radii] = 1.
    center_distances[center_distances > ellipse_radii] = np.nan
    center_distances *= favorability_map_array[ellipse_slice]
    
    return center_distances

# Computation of the landing ellipse uncertainty

@jit(nopython = True)
def compute_ellipse_coordinates(rad_center_longitude, rad_center_latitude, a, b, azimut, 
                                planet_radius, 
                                number_of_nodes = 100):
    longitudes = [-99999.]
    latitudes = [-99999.]
    bearing = 0.
    for node in range(number_of_nodes):
        bearing = 2*math.pi*node/(number_of_nodes - 1)
        radius = math.sqrt(2)*a*b/(math.sqrt((b**2 - a**2)*math.cos(2*bearing - 2*azimut) + a**2 + b**2))
        center_nvector = nvector_from_lonlat_math(rad_center_longitude, rad_center_latitude)
        nvector_1 = compute_great_circle_nvector_math(center_nvector,
                                                      bearing, radius, 
                                                      planet_radius)
        longitude, latitude = lonlat_from_nvector_math(nvector_1)
        longitudes.append(longitude)
        latitudes.append(latitude)
        
    return longitudes[1:], latitudes[1:]

@jit(nopython = True)
def min_list(list_a):
    min_value = list_a[0]
    for value in list_a:
        if value < min_value:
            min_value = value
    return min_value

@jit(nopython = True)
def max_list(list_a):
    max_value = list_a[0]
    for value in list_a:
        if value > max_value:
            max_value = value
    return max_value

@jit(nopython = True)
def compute_ellipse_extremities(ellipse_path_longitudes, ellipse_path_latitudes):
    longitudes_1 = [-99999.]
    latitudes_1 = [-99999.]
    longitudes_2 = [-99999.]
    latitudes_2 = [-99999.]
    longitudes_3 = [-99999.]
    latitudes_3 = [-99999.]
    previous_lon_sign = 1
    if abs(ellipse_path_longitudes[0]) != 0.:
        previous_lon_sign = ellipse_path_longitudes[0]/abs(ellipse_path_longitudes[0])
    list_status = 0
    for i in range(len(ellipse_path_longitudes)):
        lon_sign = previous_lon_sign
        if abs(ellipse_path_longitudes[i]) != 0.:
            lon_sign = ellipse_path_longitudes[i]/abs(ellipse_path_longitudes[i])
        if lon_sign != previous_lon_sign:
            list_status += 1
            previous_lon_sign = lon_sign
        if list_status == 0:
            longitudes_1.append(ellipse_path_longitudes[i])
            latitudes_1.append(ellipse_path_latitudes[i])
        elif list_status == 1:
            longitudes_2.append(ellipse_path_longitudes[i])
            latitudes_2.append(ellipse_path_latitudes[i])
        elif list_status == 2:
            longitudes_3.append(ellipse_path_longitudes[i])
            latitudes_3.append(ellipse_path_latitudes[i])
    longitudes_1 = longitudes_1[1:][::-1] + longitudes_3[1:][::-1]
    latitudes_1 = latitudes_1[1:][::-1] + latitudes_3[1:][::-1]
    longitudes_2 = longitudes_2[1:]
    latitudes_2 = latitudes_2[1:]

    extremities_1 = [(-99999., -99999.), (-99999., -99999.)]
    extremities_2 = [(-99999., -99999.), (-99999., -99999.)]
    if (abs(longitudes_1[0] - longitudes_2[0]) > 359 and
        abs(longitudes_1[-1] - longitudes_2[-1]) > 359 and
        abs(latitudes_1[0] - latitudes_2[0]) < 1 and
        abs(latitudes_1[-1] - latitudes_2[-1]) < 1):
        extremities_1[0] = (min_list(longitudes_1), min_list(latitudes_1))
        extremities_1[1] = (max_list(longitudes_1), max_list(latitudes_1))
        extremities_2[0] = (min_list(longitudes_2), min_list(latitudes_2))
        extremities_2[1] = (max_list(longitudes_2), max_list(latitudes_2))
    elif (((abs(longitudes_1[0] - longitudes_2[0]) > 359 and
            abs(longitudes_1[-1] - longitudes_2[-1]) < 1) or
           (abs(longitudes_1[0] - longitudes_2[0]) < 1 and
            abs(longitudes_1[-1] - longitudes_2[-1]) > 359)) and
          abs(latitudes_1[0] - latitudes_2[0]) < 1 and
          abs(latitudes_1[-1] - latitudes_2[-1]) < 1):
        min_latitude = min_list(ellipse_path_latitudes)
        max_latitude = 90
        if latitudes_1[0] < 0:
            min_latitude = -90
            max_latitude = max_list(ellipse_path_latitudes)
        extremities_1[0] = (-180, min_latitude)
        extremities_1[1] = (180, max_latitude)
    else:
        extremities_1[0] = (min_list(ellipse_path_longitudes),
                            min_list(ellipse_path_latitudes))
        extremities_1[1] = (max_list(ellipse_path_longitudes),
                            max_list(ellipse_path_latitudes))

    return extremities_1, extremities_2

@jit(nopython = True)
def compute_ellipse_bounding_box(ellipse_extremities,
                                 lon_min, lon_max, lat_min, lat_max, 
                                 raster_width, raster_height):
    if len(ellipse_extremities) > 0:
        i_0, j_0 = transform_to_pixel_coordinates_math(ellipse_extremities[0][0], ellipse_extremities[0][1], 
                                                       lon_min, lon_max, lat_min, lat_max, 
                                                       raster_width, raster_height)
        if j_0 > raster_height - 1:
            j_0 = raster_height - 1
        i_1, j_1 = transform_to_pixel_coordinates_math(ellipse_extremities[1][0], ellipse_extremities[1][1], 
                                                       lon_min, lon_max, lat_min, lat_max, 
                                                       raster_width, raster_height)
        if j_1 < 0:
            j_1 = 0
        return (i_0, i_1 + 1), (j_1, j_0 + 1)
    else:
        return (-99999, -99999), (-99999, -99999)

@jit(nopython = True)
def transform_to_pixel_coordinates_math(x, y,
                                        xmin, xmax,
                                        ymin, ymax,
                                        width, height):
    u = math.floor((x - xmin)*(width/(xmax - xmin)))
    v = math.floor((ymax - y)*(height/(ymax - ymin)))
    return int(u), int(v)

@jit(nopython = True)
def compute_ellipse_and_bounding_box(center_longitude, center_latitude,
                                     a, b, azimut, 
                                     lon_min, lon_max, lat_min, lat_max, 
                                     raster_width, raster_height, 
                                     planet_radius, number_of_nodes = 100):
    longitudes, latitudes = compute_ellipse_coordinates(center_longitude, center_latitude, a, b, azimut, 
                                                        planet_radius, 
                                                        number_of_nodes = number_of_nodes)
    extremities_1, extremities_2 = compute_ellipse_extremities(longitudes, latitudes)

    slice_1_i, slice_1_j = compute_ellipse_bounding_box(extremities_1,
                                                        lon_min, lon_max, lat_min, lat_max, 
                                                        raster_width, raster_height)
    slice_2_i, slice_2_j = compute_ellipse_bounding_box(extremities_2,
                                                        lon_min, lon_max, lat_min, lat_max, 
                                                        raster_width, raster_height)
    
    return slice_1_i, slice_1_j, slice_2_i, slice_2_j

@jit(nopython = True)
def get_favorability_inside_ellipse(favorability_map_array, 
                                    rad_center_longitude, rad_center_latitude, 
                                    rad_longitude_array, rad_latitude_array,
                                    planet_radius,
                                    a, b, azimut,
                                    slice_i, slice_j):
    ellipse_favorability = [-99999.]
    if slice_i[0] != -99999:
        for j in range(slice_j[0], slice_j[1]):
            for i in range(slice_i[0], slice_i[1]):

                center_distance, center_bearing = compute_great_circle_distance_and_bearing_math(rad_center_longitude, 
                                                                                                 rad_center_latitude, 
                                                                                                 rad_longitude_array[j, i], 
                                                                                                 rad_latitude_array[j, i], 
                                                                                                 planet_radius)
                ellipse_radius = math.sqrt(2)*a*b/(math.sqrt((b**2 - a**2)*math.cos(2*center_bearing - 2*azimut) + a**2 + b**2))
                if center_distance <= ellipse_radius:
                    ellipse_favorability.append(favorability_map_array[j, i])
    return ellipse_favorability[1:]

@jit(nopython = True)
def compute_number_of_ellipse_nodes(latitude, 
                                    min_number_of_nodes = 100, 
                                    max_number_of_nodes = 500, 
                                    sigmoid_midlatitude = 85,
                                    steepness = 0.75):
    number_of_nodes = min_number_of_nodes + \
        (max_number_of_nodes - min_number_of_nodes)/(1 + np.exp(-steepness*(abs(latitude) - abs(sigmoid_midlatitude))))
    return int(number_of_nodes)

@jit(nopython = True)
def compute_landing_ellipse_uncertainty(raster_rawfavorability_array,
                                        i,
                                        j,
                                        rad_longitude_array, 
                                        rad_latitude_array, 
                                        a, 
                                        b,
                                        azimuth,
                                        min_number_of_nodes = 100, 
                                        max_number_of_nodes = 500, 
                                        sigmoid_midlatitude = 85,
                                        steepness = 0.75,
                                        raster_lon_min = -180, 
                                        raster_lon_max = 180, 
                                        raster_lat_min = -90, 
                                        raster_lat_max = 90,    
                                        planet_radius = 3389.50):
    
    rad_center_longitude = rad_longitude_array[j, i]
    rad_center_latitude = rad_latitude_array[j, i]

    number_of_nodes = compute_number_of_ellipse_nodes(rad_center_latitude*180./math.pi, 
                                                      min_number_of_nodes,
                                                      max_number_of_nodes,
                                                      sigmoid_midlatitude,
                                                      steepness)

    slice_1_i, slice_1_j, \
    slice_2_i, slice_2_j = compute_ellipse_and_bounding_box(rad_center_longitude, rad_center_latitude, 
                                                            a, 
                                                            b, 
                                                            azimuth, 
                                                            raster_lon_min, 
                                                            raster_lon_max, 
                                                            raster_lat_min, 
                                                            raster_lat_max, 
                                                            rad_longitude_array.shape[1], 
                                                            rad_longitude_array.shape[0], 
                                                            planet_radius, 
                                                            number_of_nodes = number_of_nodes)

    center_distances_1 = get_favorability_inside_ellipse(raster_rawfavorability_array,
                                                         rad_center_longitude, 
                                                         rad_center_latitude, 
                                                         rad_longitude_array, 
                                                         rad_latitude_array,
                                                         planet_radius,
                                                         a, b, azimuth,
                                                         slice_1_i, slice_1_j)
    center_distances_2 = get_favorability_inside_ellipse(raster_rawfavorability_array,
                                                         rad_center_longitude, 
                                                         rad_center_latitude, 
                                                         rad_longitude_array, 
                                                         rad_latitude_array,
                                                         planet_radius,
                                                         a, b, azimuth,
                                                         slice_2_i, slice_2_j)
    center_distances = center_distances_1 + center_distances_2
            
    return np.mean(np.array(center_distances)), np.std(np.array(center_distances))

@jit(nopython = True)
def compute_landing_ellipse_uncertainties(raster_rawfavorability_array,
                                          ii,
                                          jj,
                                          rad_longitude_array, 
                                          rad_latitude_array, 
                                          a, 
                                          b,
                                          azimuth,
                                          min_number_of_nodes = 100, 
                                          max_number_of_nodes = 500, 
                                          sigmoid_midlatitude = 85,
                                          steepness = 0.75,
                                          raster_lon_min = -180, 
                                          raster_lon_max = 180, 
                                          raster_lat_min = -90, 
                                          raster_lat_max = 90,    
                                          planet_radius = 3389.50):
    
    ellipse_array = np.full((2, rad_longitude_array.shape[0], rad_longitude_array.shape[1]), np.nan)
    
    for j in jj:
        for i in ii:
            ellipse_array[0, j, i], ellipse_array[1, j, i] = compute_landing_ellipse_uncertainty(raster_rawfavorability_array,
                                                                                                 i,
                                                                                                 j,
                                                                                                 rad_longitude_array, 
                                                                                                 rad_latitude_array, 
                                                                                                 a, 
                                                                                                 b,
                                                                                                 azimuth,
                                                                                                 min_number_of_nodes, 
                                                                                                 max_number_of_nodes, 
                                                                                                 sigmoid_midlatitude,
                                                                                                 steepness,
                                                                                                 raster_lon_min, 
                                                                                                 raster_lon_max, 
                                                                                                 raster_lat_min, 
                                                                                                 raster_lat_max,    
                                                                                                 planet_radius)
            
    return ellipse_array