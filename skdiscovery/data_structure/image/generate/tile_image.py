from collections import OrderedDict

from skdiscovery.data_structure.framework.base import PipelineItem
from skdiscovery.utilities.patterns.image_tools import divideIntoSquares

import numpy as np


class TileImage(PipelineItem):

    def __init__(self, str_description, ap_paramList, size, min_deviation=None, min_fraction=None):

        assert min_deviation is None and min_fraction is None or \
               min_deviation is not None and min_fraction is not None, \
               'Both min_deviation and min_fraction have to be None or not None'

        self.size = size
        self._min_deviation = min_deviation
        self._min_fraction = min_fraction
        

        super(TileImage, self).__init__(str_description, ap_paramList)

    def process(self, obj_data):

        stride = self.ap_paramList[0]()

        if len(self.ap_paramList) > 1:
            threshold_function = self.ap_paramList[1]()
        else:
            threshold_function = None

        results = OrderedDict()
        metadata = OrderedDict()
        
        for label, data in obj_data.getIterator():

            extents, patches = divideIntoSquares(data, self.size, stride)

            if self._min_fraction is not None:
                valid_index = np.count_nonzero(np.abs(patches) < self._min_deviation, axis=(1,2)) / np.prod(patches.shape[1:]) > self._min_fraction
                print(np.count_nonzero(valid_index))
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
