# The MIT License (MIT)
# Copyright (c) 2017 Massachusetts Institute of Technology
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

# Standard library imports
from collections import OrderedDict

# scikit discovery imports
from skdiscovery.data_structure.framework.base import PipelineItem

# scikit data access imports
from skdaccess.utilities.support import progress_bar
from skdaccess.utilities.image_util import SplineLatLon

# Pyinsar imports
from pyinsar.processing.instruments.sentinel import retrieve_azimuth_time, read_geolocation, update_geolocation_lines

# 3rd party imports
from scipy.interpolate import SmoothBivariateSpline


class Deburst(PipelineItem):
    ''' Debursts Sentinel-1 TOPSAR data '''


    def __init__(self, str_description, cut_on_master=True):
        self._cut_on_master = True

        super(Deburst, self).__init__(str_description)

    def process(self, obj_data):
        '''
        Preprocesses sentinel 1 data

        @param obj_data: Data wrapper
        '''

        for index, (label, image) in enumerate(obj_data.getIterator()):

            tree = obj_data.info(label)['Tree']

            azimuth_time, line_index, split_indicies = retrieve_azimuth_time(tree)

            if self._cut_on_master and index==0:
                master_azimuth_time = azimuth_time
                master_line_index = line_index
                master_split_indicies = split_indicies

            elif self._cut_on_master:
                line_index = master_line_index

            obj_data.info(label)['Azimuth Time'] = azimuth_time[line_index].reset_index(drop=True)
            obj_data.info(label)['Split Indicies'] = split_indicies
            obj_data.info(label)['Line Index'] = line_index


            geo_info = read_geolocation(tree)

            updated_lines = update_geolocation_lines(tree, azimuth_time[line_index], geo_info)

            lat_spline = SmoothBivariateSpline(updated_lines,
                                               geo_info['Samples'],
                                               geo_info['Latitudes'], kx=1, ky=1)

            lon_spline = SmoothBivariateSpline(updated_lines,
                                               geo_info['Samples'],
                                               geo_info['Longitudes'], kx=1, ky=1)


            obj_data.info(label)['Geolocation'] = SplineLatLon(lat_spline, lon_spline)

            obj_data.updateData(label, image[line_index,:])
