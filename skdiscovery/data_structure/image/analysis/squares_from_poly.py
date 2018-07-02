from collections import OrderedDict

from skdiscovery.data_structure.framework.base import PipelineItem
from skdiscovery.utilities.patterns.image_tools import generateSquaresAroundPoly

import numpy as np


class SquaresFromPoly(PipelineItem):

    def __init__(self, str_description, polygon_name, size=100, stride=20, required_fraction = 0.5):

        self._polygon_name = polygon_name
        self._size = size
        self._stride = stride
        self._required_fraction = required_fraction

        super(SquaresFromPoly, self).__init__(str_description)

    def process(self, obj_data):

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
