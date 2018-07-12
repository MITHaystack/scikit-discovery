from collections import OrderedDict

from skdiscovery.data_structure.framework.base import PipelineItem

import numpy as np


class RotateImage(PipelineItem):

    def process(self, obj_data):

        new_data = OrderedDict()
        new_meta = OrderedDict()
        for label, data in obj_data.getIterator():

            for i in range(0,4):
                new_label = label + ', Rotated: ' + str(i*90)
                new_data[new_label] = np.rot90(data, i)

                new_meta[new_label] = OrderedDict()
                new_meta[new_label]['Rotated'] = i*90

        obj_data.update(new_data)
        obj_data.updateMetadata(new_meta)
