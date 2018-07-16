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

from collections import OrderedDict

from skdiscovery.data_structure.framework.base import PipelineItem
from skdiscovery.utilities.patterns.image_tools import divideIntoSquares

import numpy as np


class TileImage(PipelineItem):
    '''
    Create several smaller images from a larger image
    '''
    def __init__(self, str_description, ap_paramList, size, min_deviation=None, min_fraction=None, deviation_as_percent=False):
        '''
        Initialize TileImage item

        @param str_description: String description of item
        @param ap_paramList[stride]: Distance between neighboring tiles
        @param size: Size of tile (length of one side of a square)
        @param min_deviation = Minimum deviation to use when determining to keep tile
        @param min_fraction: Minimum fraction of pixels above min_deviation needed to keep tile
        @param deviation_as_percent: Treat min_deviation as a percentage of the max value of the original image
        '''

        if deviation_as_percent and min_deviation is None:
            raise RuntimeError('Must supply min_deviation when deviation_as_percent is True')


        self.size = size
        self._min_deviation = min_deviation
        self._min_fraction = min_fraction
        self._deviation_as_percent = deviation_as_percent
        

        super(TileImage, self).__init__(str_description, ap_paramList)

    def process(self, obj_data):
        '''
        Genrate new images by tiling input images

        @param obj_data: Input image wrapper
        '''

        stride = self.ap_paramList[0]()

        if len(self.ap_paramList) > 1:
            threshold_function = self.ap_paramList[1]()
        else:
            threshold_function = None

        if threshold_function is not None and self._min_fraction is None:
            raise RuntimeError('Must supply min_fraction with threshold function')

        if threshold_function is not None and self._min_deviation is not None:
            raise RuntimeError('Cannot supply both min_deviation and threshold function')

        results = OrderedDict()
        metadata = OrderedDict()
        
        for label, data in obj_data.getIterator():

            extents, patches = divideIntoSquares(data, self.size, stride)

            if self._deviation_as_percent:
                min_deviation = self._min_deviation * np.max(data)

            else:
                min_deviation = self._min_deviation

            if self._min_fraction is not None:

                if min_deviation is not None:
                    valid_index = np.count_nonzero(np.abs(patches) < min_deviation, axis=(1,2)) / np.prod(patches.shape[1:]) > self._min_fraction

                else:
                    threshold = threshold_function(np.abs(data))
                    threshold_data = np.abs(patches.copy())

                    threshold_data[threshold_data < threshold] = np.nan

                    valid_index = np.count_nonzero(~np.isnan(threshold_data), axis=(1,2)) / np.prod(patches.shape[1:]) > self._min_fraction


                patches = patches[valid_index]
                extents = extents[valid_index]


            try:
                metadata[label] = obj_data.info(label)
            except TypeError:
                pass
                
            for index in range(0,patches.shape[0]):
                new_label = label + ', Square: ' + str(index)
                results[new_label] = patches[index, ...]
                metadata[new_label] = OrderedDict()
                metadata[new_label]['extent'] = extents[index,...]
            

        obj_data.update(results)
        obj_data.updateMetadata(metadata)
