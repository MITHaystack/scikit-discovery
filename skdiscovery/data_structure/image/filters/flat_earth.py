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

# Standard library imports
from collections import OrderedDict

# scikit discovery imports
from skdiscovery.data_structure.framework.base import PipelineItem

# pyinsar imports
from pyinsar.processing.corrections.topography import ellipsoidal_earth_slant_ranges
from pyinsar.processing.utilities.generic import OrbitInterpolation
from pyinsar.processing.utilities.insar_simulator_utils import wrap

# 3rd party imports
from geodesy import wgs84
import numpy as np

class FlatEarth(PipelineItem):
    ''' Remove flat Earth contribution from interferogram '''

    def __init__(self, str_description, pairing='neighbor', x_range = None, y_range = None, k = 2,
                 remove_topopgraphy = False, save_correction=False):

        self._pairing = pairing
        self._x_range = x_range
        self._y_range = y_range
        self._save_correction = save_correction
        self.k = k
        
        super(FlatEarth, self).__init__(str_description)

    def process(self, obj_data):

        if self._save_correction:
            flat_earth_dict = OrderedDict()
        
        for label, data in obj_data.getIterator():

            latlon = obj_data.info(label)['Geolocation']
            wavelength = obj_data.info(label)['Wavelength']            

            image_names = ['image1', 'image2']

            if self._x_range == None:
                x_start = 0
                x_end = data.shape[1]
                
            else:
                x_start = self._x_range[0]
                x_end = self._x_range[1]
                
                
                
            if self._y_range == None:
                y_start = 0
                y_end = data.shape[0]
                
            else:
                y_start = self._y_range[0]
                y_end = self._y_range[1]

                

            slant_range_dict = OrderedDict()

            for image_name in image_names:
                orbit_interp = OrbitInterpolation(obj_data.info(label)[image_name]['Orbit'])
                az_time = obj_data.info(label)[image_name]['Azimuth Time']

                slant_range_dict[image_name] = ellipsoidal_earth_slant_ranges(az_time, latlon, orbit_interp,
                                                                              x_start, x_end, y_start, y_end)[0]

            flat_earth_inteferogram = wrap(-2 * np.pi * self.k *(slant_range_dict[image_names[0]] - slant_range_dict[image_names[1]])/wavelength)
            if self._save_correction:
                flat_earth_dict[label] = flat_earth_inteferogram

            data[y_start:y_end,x_start:x_end] = wrap(data[y_start:y_end,x_start:x_end] - flat_earth_inteferogram)

        if self._save_correction:
            obj_data.addResult(self.str_description, flat_earth_dict)
