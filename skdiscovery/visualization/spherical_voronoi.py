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

import numpy as np
import pandas as pd
import matplotlib
from matplotlib.patches import Polygon
from scipy.spatial import SphericalVoronoi
import pyproj

# utility functions for generating the spherical voronoi tesselation.

def sphericalToXYZ(lat,lon,radius=1):
    ''' 
    Convert spherical coordinates to x,y,z

    @param lat: Latitude, scalar or array
    @param lon: Longitude, scalar or array
    @param radius: Sphere's radius

    @return Numpy array of x,y,z coordinates
    '''
    phi = np.deg2rad(90.0 - lat)
    theta = np.deg2rad(lon % 360)
    x = radius * np.cos(theta)*np.sin(phi)
    y = radius * np.sin(theta)*np.sin(phi)
    z = radius * np.cos(phi)
    
    if np.isscalar(x) == False:
        return np.vstack([x,y,z]).T
    else:
        return np.array([x,y,z])


def xyzToSpherical(x,y,z):
    ''' 
    Convert x,y,z to spherical coordinates

    @param x: Cartesian coordinate x
    @param y: Cartesian coordinate y
    @param z: Cartesian coordinate z

    @return numpy array of latitude,longitude, and radius
    '''
    radius = np.sqrt(x**2 + y**2 + z**2)
    theta = np.rad2deg(np.arctan2(y,x))
    phi = np.rad2deg(np.arccos(z/radius))
    # lon = (theta + 180) % 360 - 180
    # lon = (theta + 360) % 360
    lon = theta
    lat = 90 - phi
    
    return np.array([lat,lon,radius]).T


def find_match(region_index, region_list):
    ''' 
    Find neighboring regions

    @param region_index: Numeric index of region to find matches for (number between 0 and len(vertices))
    @param region_list: list of lists of vertices that define regions

    @return Numeric indices of regions that border the region specified by region_index
    '''
    regions = region_list[region_index]
    matches = []
    num_matched_list=[]
    
    for i in range(len(region_list)):
        test_regions = region_list[i]

        num_matches = 0
        found = False
        for region in regions:
            if region in test_regions:
                num_matches += 1
                found = True
                
        if found is True:
            matches.append(i)
            
        num_matched_list.append(num_matches)

    return matches


def getVoronoiCollection(data, lat_name, lon_name, bmap = None, v_name = None, full_sphere = False, 
                         max_v=.3, min_v=-0.3, cmap = matplotlib.cm.get_cmap('jet'), test_point = None,
                         proj1=None, proj2=None, **kwargs):
    '''
    Perform a Spherical Voronoi Tessellation on the input data.

    In the case where the data is restricted to one part of the globe, a polygon will not be returned
    for all objects, as matplotlib polygons won't be able to stretch over half the globe. 

    @param data: Input pandas data frame
    @param lat_name: Name of latitude column
    @param lon_name: Name of longitude column
    @param bmap: Basemap instance used to convert from lat, lon coordinates to projection coordinates
    @param v_name: Name of value column. Use this to color each cell according to a value.
    @param full_sphere: Set to true if the data spans the entire globe. 
                        If false, a fictional point is created during tessellation and 
                        removed later to work around issues when polygons are suppose to 
                        span the over half the globe.
    @param max_v: Specify a maximum value to use when assigning values to the tessellation
    @param min_v: Specify a minimum value to use when assigning values to the tessellation
    @param cmap: Matplotlib color map to use
    @param test_point: Tuple containing the latitude and longitude of the ficitonal point to used to remove polygons that
                       wrap around the earth. If none, a point is automatically chosen
    @param proj1: PyProj projection of input coordinates
    @param proj2: PyProj projection of sphere
    @param kwargs: Extra keyword arguments are passed to SphericalVoronoi class in scipy

    @return Matplotlib patch collection of tessellation, scipy.spatial.SphericalVoronoi object, integer index of objects in patch collection.
    '''

    data = data.copy()

    if full_sphere == False:
        if test_point == None:
            test_lat = -1*np.mean(data[lat_name])
            test_lon = np.mean(data[lon_name]) + 180

        else:
            test_lat = test_point[0]
            test_lon = test_point[1]
    
        full_data = data
        full_data = pd.concat([full_data, pd.DataFrame({lat_name: test_lat, lon_name: test_lon}, 
                                                        index=[full_data.index[0]])])
        
        full_data.set_index(np.arange(len(full_data)), inplace=True)

        
    else:
        full_data = data
        
    # print(full_data.tail())

    if proj1 != None and proj2 != None:
        results = pyproj.transform(proj1, proj2, full_data[lon_name].as_matrix(), full_data[lat_name].as_matrix())
        full_data[lon_name] = results[0]
        full_data[lat_name] = results[1]

    xyz = pd.DataFrame(sphericalToXYZ(full_data[lat_name], full_data[lon_name]),columns=['x','y','z'],index=full_data.index)
    
    if v_name != None:
        full_data = pd.concat([full_data.loc[:,[lat_name,lon_name, v_name]],xyz],axis=1)
    else:
        full_data = pd.concat([full_data.loc[:,[lat_name,lon_name, v_name]],xyz],axis=1)
        
    unique_index = np.unique(full_data.loc[:,lat_name] + 1j*full_data.loc[:,lon_name],return_index=True)[1]
    
    full_data = full_data.iloc[np.sort(unique_index)]
    
    voronoi = SphericalVoronoi(full_data.loc[:,['x','y','z']].as_matrix(), **kwargs)
    
    voronoi.sort_vertices_of_regions()
    
    latlon_verts = xyzToSpherical(voronoi.vertices[:,0],voronoi.vertices[:,1], voronoi.vertices[:,2])

    if proj1 != None and proj2 != None:
        results = pyproj.transform(proj2, proj1, latlon_verts[:,1], latlon_verts[:,0])
        latlon_verts[:, 1] = results[0]
        latlon_verts[:, 0] = results[1]
    
    matches = list(map(lambda x: find_match(x, voronoi.regions), range(len(voronoi.regions))))
    
    patch_list = []
    patch_index = []

    for i, (region,match,(station,row)) in enumerate(zip(voronoi.regions,matches,
                                                         full_data.iterrows())):

        if full_sphere or (len(matches)-1) not in match:
            # Keep an index of regions in patchcollection
            patch_index.append(i)

            if bmap != None:
                xy = np.array(bmap(latlon_verts[region,1],latlon_verts[region,0])).T
            else:
                xy = np.array([latlon_verts[region,1],latlon_verts[region,0]]).T
                
            
            if v_name != None:
                value = row[v_name]
                scaled_value = (value - min_v) / (max_v - min_v)
                if scaled_value > 1:
                    scaled_value = 1.0
                elif scaled_value < 0:
                    scaled_value = 0.0

                poly = Polygon(xy, fill=True,facecolor = cmap(scaled_value),edgecolor=cmap(scaled_value))
                
            else:
                poly = Polygon(xy, fill=False)
                
            patch_list.append(poly)

    return matplotlib.collections.PatchCollection(patch_list,match_original=True), voronoi, patch_index
