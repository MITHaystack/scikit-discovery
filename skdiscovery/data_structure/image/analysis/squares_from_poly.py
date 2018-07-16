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
from skdiscovery.utilities.patterns.image_tools import generateSquaresAroundPoly

import numpy as np


class SquaresFromPoly(PipelineItem):
    '''
    Generate shapely squares that intersect with a shapely polygon
    '''

    def __init__(self, str_description, polygon_name, size=100, stride=20, required_fraction = 0.5):
        '''
        Create a pipeline item to generate a shapely squares from a polygon

        @param str_description: String description of pipeline item
        @param polygon_name: Name of polygon pipeline item
        @param size: Length of a side of the shapely squares that will be generated
        @param stride: Distance between squares
        @param required_fraction: Fraction of overlap between polygon and square
        '''

        self._polygon_name = polygon_name
        self._size = size
        self._stride = stride
        self._required_fraction = required_fraction

        super(SquaresFromPoly, self).__init__(str_description)

    def process(self, obj_data):
        '''
        Process data in an image data wrapper

        @param obj_data: Image data wrapper
        '''

        result_dict = OrderedDict()

        poly_dict = obj_data.getResults()[self._polygon_name]

        for label, data in obj_data.getIterator():
            poly = poly_dict[label]

            square_list = generateSquaresAroundPoly(poly, self._size, self._stride)

            good_square_list = []
            for square in square_list:
                if poly.intersection(square).area / square.area > self._required_fraction:
                    good_square_list.append(square)

            result_dict[label] = good_square_list



        obj_data.addResult(self.str_description, result_dict)
