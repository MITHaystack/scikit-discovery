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

import random
import numpy as np

from osgeo import gdal, ogr, osr
import csv

from shapely.geometry import Polygon
from shapely.geometry import LineString
from shapely.wkt import loads

from matplotlib.path import Path
from matplotlib.patches import PathPatch
from matplotlib import colors
from mpl_toolkits.basemap import Basemap

################################################################################
# Function definitions
################################################################################

# Get the data with OGR and csv

def open_shapefile(shapefile_path, writeable = False):
    driver = ogr.GetDriverByName('ESRI Shapefile')
    return driver.Open(shapefile_path, writeable)

def get_latitude_longitude_from_csv_file(csv_file_location, 
                                         longitude_column_index = 0, 
                                         latitude_column_index = 1,
                                         other_data_column_indexes = []):
    longitudes = []
    latitudes = []
    other_data = [[] for i in range(len(other_data_column_indexes))]
    with open(csv_file_location) as csvfile:
        has_header = csv.Sniffer().has_header(csvfile.read(1024))
        csvfile.seek(0)
        reader = csv.reader(csvfile)
        if has_header:
            next(reader)
        for row in reader:
            longitude = float(row[longitude_column_index])
            latitude = float(row[latitude_column_index])
            longitudes.append(longitude)
            latitudes.append(latitude)
            for array_index in range(len(other_data_column_indexes)):
                data_value = row[other_data_column_indexes[array_index]]
                other_data[array_index].append(data_value)
    return [longitudes, latitudes] + other_data

# Get some information about the data
  
def print_shapefile_field_names(shapefile):
    layer = shapefile.GetLayer()
    layer_definition = layer.GetLayerDefn()
    print(" _______________________________________________________________________________")
    print("/ Shapefile field names:")
    print("| ")
    for i in range(layer_definition.GetFieldCount()):
        print("|", layer_definition.GetFieldDefn(i).GetName(), end = " ")
    print("\n\_______________________________________________________________________________")

def get_field_values(shapefile, field_name):
    field_values = []
    layer = shapefile.GetLayer()
    for i in range(layer.GetFeatureCount()):
        feature = layer.GetFeature(i)
        field_values.append(feature.GetField(field_name))
    return field_values

def print_shapefile_unique_field_values(shapefile, field_name):
    unique_field_values = np.unique(get_field_values(shapefile, field_name))
    print(" _______________________________________________________________________________")
    print("/ ", field_name, ":", sep = "")
    print("| ")
    for value in unique_field_values:
        print("|", value, end = " ")
    print("\n\_______________________________________________________________________________")
    
# Build polygons from shapefile geometries

def shape_coding(ob):
    n = len(ob.coords)
    codes = np.ones(n, dtype=Path.code_type)*Path.LINETO
    codes[0] = Path.MOVETO
    return codes

def create_path_from_shape(shape):
    path = None
    if isinstance(shape, LineString):
        vertices = np.asarray(shape)
        codes = shape_coding(shape)
        path = Path(vertices, codes)
    elif isinstance(shape, Polygon):
        vertices = np.concatenate(
            [np.asarray(shape.exterior)]
            + [np.asarray(r) for r in shape.interiors])
        codes = np.concatenate(
            [shape_coding(shape.exterior)]
            + [shape_coding(r) for r in shape.interiors])
        path = Path(vertices, codes)
    return path

def get_geometry_coordinates(geometry, 
                             xy_outer_path, 
                             xy_inner_paths, 
                             basemap = None):
    # Just one geometry
    outer_path = geometry
    # If nested geometries
    if geometry.GetGeometryCount() != 0:
        outer_path = geometry.GetGeometryRef(0)
    longitudes = [outer_path.GetX(i) for i in range(outer_path.GetPointCount())]
    latitudes = [outer_path.GetY(i) for i in range(outer_path.GetPointCount())]
    x = longitudes
    y = latitudes
    if basemap is not None:
        x, y = basemap(longitudes, latitudes)
    xy_outer_path.extend(list(zip(x, y)))
    # Iterate over the nested geometries
    for count in range(1, geometry.GetGeometryCount()):
        inner_path = geometry.GetGeometryRef(count)
        longitudes = [inner_path.GetX(i) for i in range(inner_path.GetPointCount())]
        latitudes = [inner_path.GetY(i) for i in range(inner_path.GetPointCount())]
        x = longitudes
        y = latitudes
        if basemap is not None:
            x, y = basemap(longitudes, latitudes)
        xy_inner_paths.append(list(zip(x, y)))
        
def build_shape_from_geometry(geometry, 
                              basemap = None):
    xy_outer_path = []
    xy_inner_paths = []
    get_geometry_coordinates(geometry, 
                             xy_outer_path, 
                             xy_inner_paths,  
                             basemap)
    if len(xy_outer_path) > 0 and "LINESTRING" in geometry.GetGeometryName():
        shape = LineString(xy_outer_path)
    elif len(xy_outer_path) > 0 and "POLYGON" in geometry.GetGeometryName():
        shape = Polygon(xy_outer_path, xy_inner_paths)
    else:
        shape = None
    return shape

# Add the data to a basemap

def add_shape_to_map(axes, 
                     shape, 
                     legend_label, 
                     facecolor = '#cccccc',
                     alpha = 1.,
                     hatch = None,
                     edgecolor = '#999999', 
                     linewidth = 0.25,
                     linestyle = '-'):
    if shape != None:
        path = create_path_from_shape(shape)
        fill = True
        if isinstance(shape, LineString) or facecolor is None:
            fill = False
        patch = PathPatch(path,
                          fill = fill,
                          facecolor = facecolor,
                          alpha = alpha,
                          hatch = hatch,
                          edgecolor = edgecolor, 
                          linewidth = linewidth,
                          linestyle = linestyle,
                          label = legend_label)
        axes.add_patch(patch)
    
def add_geometry_to_map(axes, 
                        basemap, 
                        geometry, 
                        legend_label, 
                        facecolor = '#cccccc',
                        alpha = 1.,
                        hatch = None,
                        edgecolor = '#999999', 
                        linewidth = 0.25,
                        linestyle = '-'):
    if "MULTI" in geometry.GetGeometryName():
        for count in range(geometry.GetGeometryCount()):
            sub_geometry = geometry.GetGeometryRef(count)
            shape = build_shape_from_geometry(sub_geometry, basemap)
            add_shape_to_map(axes, 
                             shape, 
                             legend_label, 
                             facecolor = facecolor,
                             alpha = alpha,
                             hatch = hatch,
                             edgecolor = edgecolor,
                             linewidth = linewidth,
                             linestyle = linestyle)
    else:
        shape = build_shape_from_geometry(geometry, basemap)
        add_shape_to_map(axes, 
                         shape, 
                         legend_label, 
                         facecolor = facecolor,
                         alpha = alpha,
                         hatch = hatch,
                         edgecolor = edgecolor, 
                         linewidth = linewidth,
                         linestyle = linestyle)
    
def add_vector_to_map(axes,
                      basemap,
                      shapefile,
                      field_name,
                      random_colors = False,
                      facecolor = '#08519c',
                      alpha = 1.,
                      hatch = None,
                      edgecolor = '#252525', 
                      linewidth = 0.25,
                      linestyle = '-'):
    layer = shapefile.GetLayer()
    for i in range(layer.GetFeatureCount()):
        feature = layer.GetFeature(i)
        feature_geometry = feature.GetGeometryRef()
        if feature_geometry != None:
            final_facecolor = facecolor
            if random_colors == True:
                r = lambda: random.randint(0,255)
                final_facecolor = '#%02X%02X%02X'%(r(), r(), r())
            legend_label = None
            if field_name != None:
                legend_label = feature.GetField(field_name)
            add_geometry_to_map(axes, basemap, feature_geometry, legend_label,
                                facecolor = final_facecolor, alpha = alpha, hatch = hatch,
                                edgecolor = edgecolor, linewidth = linewidth,
                                linestyle = linestyle)
        
def add_path_to_map(axes, 
                    path, 
                    legend_label, 
                    facecolor = '#cccccc',
                    alpha = 1.,
                    edgecolor = '#999999',
                    linestyle = '-',
                    linewidth = 0.25,
                    zorder = 1):
    patch = PathPatch(path, 
                      facecolor = facecolor,
                      alpha = alpha,
                      edgecolor = edgecolor, 
                      linewidth = linewidth,
                      linestyle = linestyle,
                      label = legend_label,
                      zorder = zorder)
    axes.add_patch(patch)

# Modify Shapefiles

def filter_shapefile(shapefile, 
                     field_name, 
                     field_filter_values, 
                     file_type = 'Memory',
                     file_path = '',
                     geom_type = None):

    input_layer = shapefile.GetLayer()
    input_spatial_reference = input_layer.GetSpatialRef()
    if geom_type is None:
        geom_type = input_layer.GetGeomType()

    output_driver = ogr.GetDriverByName(file_type)
    output_shapefile = output_driver.CreateDataSource(file_path)
    output_layer = output_shapefile.CreateLayer('',
                                                srs = input_spatial_reference,
                                                geom_type = geom_type)
    output_layer_defn = output_layer.GetLayerDefn()
    
    for i in range(input_layer.GetFeatureCount()):
        input_feature = input_layer.GetFeature(i)
        if (input_feature is not None
            and input_feature.GetField(field_name) in field_filter_values):
            input_geometry = input_feature.GetGeometryRef()
            output_feature = ogr.Feature(output_layer_defn)
            output_feature.SetGeometry(input_geometry)
            output_layer.CreateFeature(output_feature)
        
    return output_shapefile

def get_shapefile_borders(shapefile, 
                          file_type = 'Memory',
                          file_path = '',
                          geom_type = ogr.wkbLineString):

    input_layer = shapefile.GetLayer()
    input_spatial_reference = input_layer.GetSpatialRef()

    output_driver = ogr.GetDriverByName(file_type)
    output_shapefile = output_driver.CreateDataSource(file_path)
    output_layer = output_shapefile.CreateLayer('',
                                                srs = input_spatial_reference,
                                                geom_type = geom_type)
    output_layer_defn = output_layer.GetLayerDefn()
    
    for i in range(input_layer.GetFeatureCount()):
        input_feature = input_layer.GetFeature(i)
        if input_feature is not None:
            input_geometry = input_feature.GetGeometryRef()
            input_boundary = input_geometry.Boundary()
            output_feature = ogr.Feature(output_layer_defn)
            output_feature.SetGeometry(input_boundary)
            output_layer.CreateFeature(output_feature)
        
    return output_shapefile

def buffer_shapefile(shapefile,
                     buffer_distance, 
                     file_type = 'Memory',
                     file_path = '',
                     geom_type = ogr.wkbPolygon):

    input_layer = shapefile.GetLayer()
    input_spatial_reference = input_layer.GetSpatialRef()

    output_driver = ogr.GetDriverByName(file_type)
    output_shapefile = output_driver.CreateDataSource(file_path)
    output_layer = output_shapefile.CreateLayer('',
                                                srs = input_spatial_reference,
                                                geom_type = geom_type)
    output_layer_defn = output_layer.GetLayerDefn()

    for i in range(input_layer.GetFeatureCount()):
        feature = input_layer.GetFeature(i)
        input_geometry = feature.GetGeometryRef()
        buffer_geometry = input_geometry.Buffer(buffer_distance)

        output_feature = ogr.Feature(output_layer_defn)
        output_feature.SetGeometry(buffer_geometry)
        output_layer.CreateFeature(output_feature)
        output_feature = None
        
    return output_shapefile

def clip_shapefile(shapefile, 
                   polygon_clip, 
                   file_type = 'Memory',
                   file_path = '',
                   geom_type = ogr.wkbPolygon):

    input_layer = shapefile.GetLayer()
    input_spatial_reference = input_layer.GetSpatialRef()

    output_driver = ogr.GetDriverByName(file_type)
    output_shapefile = output_driver.CreateDataSource(file_path)
    output_layer = output_shapefile.CreateLayer('',
                                                srs = input_spatial_reference,
                                                geom_type = geom_type)

    # Add the fields to the new shapefile
    input_layer_defn = input_layer.GetLayerDefn()
    for i in range(input_layer_defn.GetFieldCount()):
        field_defn = input_layer_defn.GetFieldDefn(i)
        output_layer.CreateField(field_defn)

    # Add the features to the new shapefile
    output_layer_defn = output_layer.GetLayerDefn()
    for i in range(input_layer.GetFeatureCount()):
        input_feature = input_layer.GetFeature(i)
        geometry = input_feature.GetGeometryRef()
        geometry = loads(geometry.ExportToWkt())
        geometry = geometry.intersection(polygon_clip)

        if geometry.is_empty == False:
            output_feature = ogr.Feature(output_layer_defn)
            output_geometry = ogr.CreateGeometryFromWkt(geometry.wkt)
            output_feature.SetGeometry(output_geometry)
            for i in range(0, output_layer_defn.GetFieldCount()):
                output_feature.SetField(output_layer_defn.GetFieldDefn(i).GetNameRef(), 
                                        input_feature.GetField(i))
            output_layer.CreateFeature(output_feature)
        
    return output_shapefile

def union_shapefiles(shapefile_1, 
                     shapefile_2, 
                     file_type = 'Memory',
                     file_path = '',
                     geom_type = ogr.wkbPolygon):

    input_layer_1 = shapefile_1.GetLayer()
    input_layer_2 = shapefile_2.GetLayer()
    input_spatial_reference = input_layer_1.GetSpatialRef()

    output_driver = ogr.GetDriverByName(file_type)
    output_shapefile = output_driver.CreateDataSource(file_path)
    output_layer = output_shapefile.CreateLayer('',
                                                srs = input_spatial_reference,
                                                geom_type = geom_type)
    
    input_layer_1.Union(input_layer_2, output_layer)
        
    return output_shapefile

def intersect_shapefiles(shapefile_1, 
                         shapefile_2, 
                         file_type = 'Memory',
                         file_path = '',
                         geom_type = ogr.wkbPolygon):

    input_layer_1 = shapefile_1.GetLayer()
    input_layer_2 = shapefile_2.GetLayer()
    input_spatial_reference = input_layer_1.GetSpatialRef()

    output_driver = ogr.GetDriverByName(file_type)
    output_shapefile = output_driver.CreateDataSource(file_path)
    output_layer = output_shapefile.CreateLayer('',
                                                srs = input_spatial_reference,
                                                geom_type = geom_type)
    
    input_layer_1.Intersection(input_layer_2, output_layer)
        
    return output_shapefile

def get_intersected_features_from_shapefile(input_shapefile,
                                            method_shapefile,
                                            look_for_intersection = True,
                                            file_type = 'Memory',
                                            file_path = ''):

    input_layer = input_shapefile.GetLayer()
    method_layer = method_shapefile.GetLayer()
    input_spatial_reference = input_layer.GetSpatialRef()
    geom_type = input_layer.GetGeomType()

    output_driver = ogr.GetDriverByName(file_type)
    output_shapefile = output_driver.CreateDataSource(file_path)
    output_layer = output_shapefile.CreateLayer('',
                                                srs = input_spatial_reference,
                                                geom_type = geom_type)
    output_layer_defn = output_layer.GetLayerDefn()

    for i in range(input_layer.GetFeatureCount()):
        input_feature = input_layer.GetFeature(i)
        if input_feature is not None:
            input_geometry = input_feature.GetGeometryRef()
            count_not_intersected_features = 0
            for i in range(method_layer.GetFeatureCount()):
                method_feature = method_layer.GetFeature(i)
                if method_feature is not None:
                    method_geometry = method_feature.GetGeometryRef()
                    if input_geometry.Intersects(method_geometry) == False:
                        count_not_intersected_features += 1

            if ((look_for_intersection == False
                 and count_not_intersected_features == method_layer.GetFeatureCount())
                or (look_for_intersection == True
                    and count_not_intersected_features < method_layer.GetFeatureCount())):
                output_feature = ogr.Feature(output_layer_defn)
                output_feature.SetGeometry(input_geometry)
                output_layer.CreateFeature(output_feature)
                    
    return output_shapefile

def modify_shapefile_extent(shapefile, 
                            x_min, x_max, y_min, y_max, 
                            new_x_min, new_x_max, new_y_min, new_y_max, 
                            file_type = 'Memory',
                            file_path = '',
                            geom_type = ogr.wkbPolygon):

    input_layer = shapefile.GetLayer()
    input_spatial_reference = input_layer.GetSpatialRef()

    output_driver = ogr.GetDriverByName(file_type)
    output_shapefile = output_driver.CreateDataSource(file_path)
    output_layer = output_shapefile.CreateLayer('',
                                                srs = input_spatial_reference,
                                                geom_type = geom_type)

    # Add the fields to the new shapefile
    input_layer_defn = input_layer.GetLayerDefn()
    for i in range(input_layer_defn.GetFieldCount()):
        field_defn = input_layer_defn.GetFieldDefn(i)
        output_layer.CreateField(field_defn)

    # Add the features to the new shapefile
    output_layer_defn = output_layer.GetLayerDefn()
    for i in range(input_layer.GetFeatureCount()):
        input_feature = input_layer.GetFeature(i)
        geometry = input_feature.GetGeometryRef()

        # Just one geometry
        outer_ring = geometry
        # If nested geometries
        if geometry.GetGeometryCount() != 0:
            outer_ring = geometry.GetGeometryRef(0)
        for i in range(outer_ring.GetPointCount()):
            x = outer_ring.GetX(i)
            y = outer_ring.GetY(i)
            new_x = new_x_min + (new_x_max - new_x_min)*(x - x_min)/(x_max - x_min)
            new_y = new_y_max - (new_y_max - new_y_min)*(y_max - y)/(y_max - y_min)
            outer_ring.SetPoint(i, new_x, new_y)
        # Iterate over the nested geometries
        for count in range(1, geometry.GetGeometryCount()):
            inner_ring = geometry.GetGeometryRef(count)
            for i in range(inner_ring.GetPointCount()):
                x = inner_ring.GetX(i)
                y = inner_ring.GetY(i)
                new_x = new_x_min + (new_x_max - new_x_min)*(x - x_min)/(x_max - x_min)
                new_y = new_y_max - (new_y_max - new_y_min)*(y_max - y)/(y_max - y_min)
                inner_ring.SetPoint(i, new_x, new_y)

        output_feature = ogr.Feature(output_layer_defn)
        output_feature.SetGeometry(geometry)
        for i in range(0, output_layer_defn.GetFieldCount()):
            output_feature.SetField(output_layer_defn.GetFieldDefn(i).GetNameRef(), 
                                    input_feature.GetField(i))
        output_layer.CreateFeature(output_feature)
        
#     output_shapefile.Destroy()
        
    return output_shapefile

# Rasterization

def rasterize_geometries(shapes,
                         data_type,
                         raster_x_size, 
                         raster_y_size,
                         geotransform,
                         spatial_reference,
                         fill_value = 0,
                         background_value = 1,
                         no_data_value = -99999, 
                         scale = 1, 
                         offset = 0,
                         all_touched = False,
                         file_type = 'MEM',
                         file_path = '',
                         number_of_bands = 1):
    
    shape_driver = ogr.GetDriverByName('Memory')
    shape_dataset = shape_driver.CreateDataSource('')
    shape_layer = shape_dataset.CreateLayer('', srs = spatial_reference) #, geom_type = ogr.wkbPolygon
    shape_layer_defn = shape_layer.GetLayerDefn()
        
    if isinstance(shapes, list):
        for shape in shapes:
            feature = ogr.Feature(shape_layer_defn)
            geometry = ogr.CreateGeometryFromWkt(shape.wkt)
            feature.SetGeometry(geometry)
            shape_layer.CreateFeature(feature)
    elif isinstance(shapes, Polygon):
        feature = ogr.Feature(shape_layer_defn)
        geometry = ogr.CreateGeometryFromWkt(shapes.wkt)
        feature.SetGeometry(geometry)
        shape_layer.CreateFeature(feature)
        
    # MEM: temporary memory file
    raster_driver = gdal.GetDriverByName(file_type)
    raster_dataset = raster_driver.Create(file_path,
                                          raster_x_size,
                                          raster_y_size,
                                          number_of_bands,
                                          data_type)
    raster_dataset.SetProjection(spatial_reference.ExportToWkt())
    raster_dataset.SetGeoTransform(geotransform)

    raster_band = raster_dataset.GetRasterBand(1)
    raster_band.SetNoDataValue(no_data_value)
    raster_band.SetScale(scale)
    raster_band.SetOffset(offset)
    raster_band.Fill(background_value)
    
    rasterize_options = []
    if all_touched == True:
        rasterize_options.append("ALL_TOUCHED=TRUE")

    error = gdal.RasterizeLayer(raster_dataset, 
                                [1], 
                                shape_layer, 
                                burn_values = [fill_value], 
                                options = rasterize_options)

    if error != 0:
        print('Non-zero gdal.RasterizeLayer error:', error)
        
    shape_layer = None
    
    return raster_dataset

def rasterize_shapefile(shapefile, 
                        field_name,
                        data_type,
                        raster_x_size, 
                        raster_y_size,            
                        geotransform, 
                        projection,
                        fill_value = 0,
                        background_value = 1,
                        no_data_value = -99999, 
                        scale = 1, 
                        offset = 0,
                        file_type = 'MEM',
                        file_path = '',
                        number_of_bands = 1):
    layer = shapefile.GetLayer()
    
    # MEM: temporary memory file
    driver = gdal.GetDriverByName(file_type)
    raster_dataset = driver.Create(file_path, raster_x_size, raster_y_size, number_of_bands, data_type)
    raster_dataset.SetProjection(projection)
    raster_dataset.SetGeoTransform(geotransform)

    raster_band = raster_dataset.GetRasterBand(1)
    raster_band.SetNoDataValue(no_data_value)
    raster_band.SetScale(scale)
    raster_band.SetOffset(offset)
    raster_band.Fill(background_value)
    if field_name is not None:
        #raster_band.Fill(raster_band.GetNoDataValue())
        error = gdal.RasterizeLayer(raster_dataset, [1], layer, options = ["ATTRIBUTE=" + field_name])
    else:
        #raster_band.Fill(background_value)
        error = gdal.RasterizeLayer(raster_dataset, [1], layer, burn_values = [fill_value])
    
    if error != 0:
        print('Non-zero gdal.RasterizeLayer error:', error)
    
    return raster_dataset