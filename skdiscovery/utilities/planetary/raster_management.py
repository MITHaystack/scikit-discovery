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

from osgeo import gdal, osr
from gdalconst import GA_ReadOnly, GA_Update

from mpl_toolkits.basemap import Basemap

################################################################################
# Function definitions
################################################################################

# Get the data with GDAL

def open_raster(gdal_raster_path, read_only = True):
    read_status = GA_ReadOnly
    if read_only == False:
        read_status = GA_Update
    gdal_raster = gdal.Open(gdal_raster_path, read_status)
    if gdal_raster is None:
        print('Could not open file')
    return gdal_raster

def get_raster_array(gdal_raster, remove_ndv = True):
    assert gdal_raster is not None, 'No raster available'

    number_of_bands = gdal_raster.RasterCount
    raster_array = gdal_raster.ReadAsArray().astype(np.float)
    for i_band in range(number_of_bands):
        raster_band = gdal_raster.GetRasterBand(i_band + 1)
        no_data_value = raster_band.GetNoDataValue()
        if no_data_value is not None and remove_ndv == True:
            if number_of_bands > 1:
                raster_array[i_band, :, :][raster_array[i_band, :, :] == no_data_value] = np.nan
            else:
                raster_array[raster_array == no_data_value] = np.nan
        scale = raster_band.GetScale()
        if scale is None:
            scale = 1.
        offset = raster_band.GetOffset()
        if offset is None:
            offset = 0.
        if number_of_bands > 1:
            raster_array[i_band, :, :] = raster_array[i_band, :, :]*scale + offset
        else:
            raster_array = raster_array*scale + offset
            
    return raster_array

# Get some information about the data

def get_raster_extent(gdal_raster):
    assert gdal_raster is not None, 'No raster available'

    raster_x_size = gdal_raster.RasterXSize
    raster_y_size = gdal_raster.RasterYSize
    geotransform = gdal_raster.GetGeoTransform()
    xmin = geotransform[0]
    ymax = geotransform[3]
    xmax = xmin + geotransform[1]*raster_x_size
    ymin = ymax + geotransform[5]*raster_y_size
    return xmin, xmax, ymin, ymax

def print_raster_info(gdal_raster):
    assert gdal_raster is not None, 'No raster available'

    print('Driver: ', gdal_raster.GetDriver().ShortName, '/', gdal_raster.GetDriver().LongName)
    print('Size of the cube is ', gdal_raster.RasterXSize, 'x', gdal_raster.RasterYSize, 'x', gdal_raster.RasterCount)
    print('Projection is ', gdal_raster.GetProjection())
    geotransform = gdal_raster.GetGeoTransform()
    if not geotransform is None:
        print('Origin = (', geotransform[0], ',', geotransform[3], ')')
        print('Pixel Size = (', geotransform[1], ',', geotransform[5], ')')

def define_geotransform(xmin, xmax, ymin, ymax, raster_x_size, raster_y_size):
    pixel_x_size = (xmax - xmin)/raster_x_size
    pixel_y_size = (ymax - ymin)/raster_y_size
    return (xmin, pixel_x_size, 0, ymax, 0, -pixel_y_size)

# Add the data to a basemap

class DiscreteColormap:
    def __init__(self, cmap, norm, boundaries, ticks):
        self.cmap = cmap
        self.norm = norm
        self.boundaries = boundaries
        self.ticks = ticks

def add_raster_to_map(basemap,
                      raster_array,
                      raster_name,
                      min_longitude = -180, max_longitude = 180, 
                      min_latitude = -90, max_latitude = 90,
                      colormap = 'viridis',
                      add_colorbar = True,
                      zorder = 1,
                      use_latlon = True,
                      use_pcolormesh = True):
    # TODO: the coordinates arrays are one cell larger than the data array. 
    # This is the only solution that appears to place the whole GDAL raster correctly 
    # in the basemap. But it does not work with latlon set to True and min_longitude < -180.
    longitudes = np.linspace(min_longitude, max_longitude, raster_array.shape[1] + 1)
    latitudes = np.linspace(max_latitude, min_latitude, raster_array.shape[0] + 1)
    longitudes, latitudes = np.meshgrid(longitudes, latitudes)
    if use_latlon == False:
        longitudes, latitudes = basemap(longitudes, latitudes)
    masked_raster_array = np.ma.masked_invalid(raster_array)
    
    cmap = colormap
    norm = None
    boundaries = None
    ticks = None
    if isinstance(colormap, DiscreteColormap ) == True:
        cmap = colormap.cmap
        norm = colormap.norm
        boundaries = colormap.boundaries
        ticks = colormap.ticks
    
    if use_pcolormesh == True:
        raster_map = basemap.pcolormesh(longitudes, 
                                        latitudes, 
                                        masked_raster_array,
                                        latlon = use_latlon, 
                                        cmap = cmap,
                                        norm = norm,
                                        zorder = zorder)
    else:
        raster_map = basemap.pcolor(longitudes, 
                                    latitudes, 
                                    masked_raster_array,
                                    latlon = use_latlon, 
                                    cmap = cmap,
                                    norm = norm,
                                    zorder = zorder)
    if add_colorbar == True:
        raster_map_colorbar = basemap.colorbar(raster_map, 
                                               cmap = cmap, 
                                               norm = norm, 
                                               boundaries = boundaries, 
                                               ticks = ticks)
        raster_map_colorbar.set_label(raster_name)

# Modify Rasters

def create_raster_from_array(raster_array,           
                             geotransform, 
                             projection,
                             file_type = 'MEM',
                             file_path = '',
                             data_type = gdal.GDT_Float64,
                             no_data_value = -99999.,
                             scale = 1.,
                             offset = 0.,
                             options = []):
    raster_x_size = raster_array.shape[1]
    raster_y_size = raster_array.shape[0]
    number_of_bands = 1
    if len(raster_array.shape) >= 3:
        raster_x_size = raster_array.shape[2]
        raster_y_size = raster_array.shape[1]
        number_of_bands = raster_array.shape[0]

    driver = gdal.GetDriverByName(file_type)
    new_raster = driver.Create(file_path,
                               raster_x_size,
                               raster_y_size,
                               number_of_bands,
                               data_type,
                               options = options)
    new_raster.SetGeoTransform(geotransform)
    new_raster.SetProjection(projection)
    
    for band_number in range(1, number_of_bands + 1):
        new_raster_band = new_raster.GetRasterBand(band_number)
        new_raster_band.SetNoDataValue(no_data_value)
        new_raster_band.SetScale(scale)
        new_raster_band.SetOffset(offset)
        # Fill the raster band, otherwise no data values are not set in the new raster
        new_raster_band.Fill(no_data_value)

        if len(raster_array.shape) >= 3:
            new_raster_band.WriteArray(raster_array[band_number - 1, :, :])
        else:
            new_raster_band.WriteArray(raster_array)

    return new_raster

def transform_to_i_coordinate(x, xmin, xmax, width):
    return int((x - xmin)*(width/(xmax - xmin)))

def recenter_raster_array(raster_array, 
                          old_central_meridian, 
                          new_central_meridian, 
                          old_lon_min, 
                          old_lon_max):

    raster_width = raster_array.shape[1]
    old_central_i = transform_to_i_coordinate(old_central_meridian, 
                                              old_lon_min, 
                                              old_lon_max, 
                                              raster_width)
    new_central_i = transform_to_i_coordinate(new_central_meridian, 
                                              old_lon_min, 
                                              old_lon_max, 
                                              raster_width)
    inter_meridian_width = abs(old_central_i - new_central_i)

    i_first = 0
    i_middle = raster_array.shape[1] - inter_meridian_width
    i_last = raster_array.shape[1]
    temp1 = np.copy(raster_array[:, i_first:i_middle])
    temp2 = np.copy(raster_array[:, i_middle:i_last])

    return np.concatenate((temp2, temp1), axis = 1)

def recenter_raster(raster,
                    old_central_meridian, 
                    new_central_meridian, 
                    old_lon_min, 
                    old_lon_max,
                    file_type = 'MEM',
                    file_path = ''):
    
    raster_band = raster.GetRasterBand(1)
    
    raster_width = raster.RasterXSize
    raster_height = raster.RasterYSize
    raster_geotransform = raster.GetGeoTransform()
    spatial_reference_system = osr.SpatialReference(wkt = raster.GetProjection())
    spatial_reference_system.SetProjParm("central_meridian", new_central_meridian)
    raster_projection = spatial_reference_system.ExportToWkt()

    driver = gdal.GetDriverByName(file_type)
    new_raster = driver.Create(file_path,
                               raster_width,
                               raster_height,
                               raster.RasterCount,
                               raster_band.DataType)
    new_raster.SetGeoTransform(raster_geotransform)
    new_raster.SetProjection(raster_projection)
    
    for band_number in range(1, raster.RasterCount + 1):
        raster_band = raster.GetRasterBand(band_number)
        new_raster_band = new_raster.GetRasterBand(band_number)
        if raster_band.GetNoDataValue() is not None:
            new_raster_band.SetNoDataValue(raster_band.GetNoDataValue())
        else:
            new_raster_band.SetNoDataValue(-99999.)
        new_raster_band.SetScale(raster_band.GetScale())
        new_raster_band.SetOffset(raster_band.GetOffset())
        # Fill the raster band, otherwise no data values are not set in the new raster
        new_raster_band.Fill(new_raster_band.GetNoDataValue())
        
        raster_array = raster_band.ReadAsArray()
        new_raster_array = recenter_raster_array(raster_array, 
                                                 old_central_meridian, 
                                                 new_central_meridian, 
                                                 old_lon_min, 
                                                 old_lon_max)

        new_raster_band.WriteArray(new_raster_array)
        
    if file_type != 'MEM':
        new_raster = None
        new_raster = open_raster(file_path)
        print_raster_info(new_raster)
    return new_raster