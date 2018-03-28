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
from pyinsar.processing.instruments.sentinel import SentinelRamp, selectValidLines, transformSLC


# 3rd party imports
import cv2
import imreg_dft as ird
import numpy as np

class Coregister(PipelineItem):

    def __init__(self, str_description, image_limits = None, num_iterations=3):
        self._image_limits = image_limits
        self._num_iterations = num_iterations
        super(Coregister, self).__init__(str_description)

    def process(self, obj_data):

        master_burst_list = None
        for label, data in obj_data.getIterator():
            if master_burst_list == None:

                master_burst_list = selectValidLines(data, obj_data.info(label)['Tree'], cut=False)
            else:

                burst_list = selectValidLines(data, obj_data.info(label)['Tree'], cut=False)
                lines_per_burst = int(obj_data.info(label)['Tree'].find('swathTiming/linesPerBurst').text)
                samples_per_burst = int(obj_data.info(label)['Tree'].find('swathTiming/samplesPerBurst').text)
                lines, samples = np.meshgrid(np.arange(lines_per_burst), np.arange(samples_per_burst), indexing = 'ij')
                ramp = SentinelRamp(obj_data.info(label))

                for index, burst in enumerate(burst_list):

                    if self._image_limits == None:
                        line_slice = slice(0, lines_per_burst)
                        sample_slice = slice(0, samples_per_burst)

                    elif self._image_limits[index] != None:
                        line_slice = self._image_limits[index][0]
                        sample_slice = self._image_limits[index][1]

                    else:
                        continue


                    shift = ird.translation(np.abs(master_burst_list[index][line_slice, sample_slice]),
                                                   np.abs(burst[line_slice, sample_slice]))

                    transform_matrix = np.array([[1, 0, shift['tvec'][1]],
                                                 [0, 1, shift['tvec'][0]]])

                    shifted_data = transformSLC(burst[line_slice, sample_slice],
                                                 -ramp(lines[line_slice, sample_slice], samples[line_slice,sample_slice], index),
                                                 transform_matrix)

                    if line_slice.start == None:
                        line_start = 0
                    elif line_slice.start < 0:
                        line_start = lines_per_burst + line_slice.start

                    if line_slice.stop == None:
                        line_end = lines_per_burst
                    elif line_slice.stop < 0:
                        line_end = lines_per_burst + line_slice.stop

                    full_data_slice = slice(lines_per_burst*index + line_start, lines_per_burst*(index) + line_end)

                    data[full_data_slice,sample_slice] = shifted_data
