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

import numpy as np

class RotateImage(PipelineItem):
    '''
    Create new images by rotating 90, 180, and 270 degrees
    '''

    def process(self, obj_data):
        '''
        Generate new images by rotate input images

        @param obj_data: Image data wrapper
        '''

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
