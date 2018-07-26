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
import os

# scikit discovery imports
from skdiscovery.data_structure.framework.base import PipelineItem

# 3rd party imports
import numpy as np
import h5py


class Saver(PipelineItem):
    ''' Write images out to a hdf5 file '''

    def __init__(self, str_description, folder_name, data_type = None):
        '''
        Initialize coherence pipeline item

        @param str_description: String identifier for item
        @param folder_name: Name to save hdf fils
        @param data_type: Data type to save data as (None defaults to input data type)
        '''

        self.folder_name = folder_name
        self.data_type = data_type

        super(Saver, self).__init__(str_description,[])


    def process(self, obj_data):
        '''
        Save images to hdf files

        @param obj_data: Data wrapper
        '''

        run_id = obj_data.getRunID()

        name = 'run_' + str(run_id).zfill(6)

        filepath = os.path.join(self.folder_name, name + '.h5')

        if len(obj_data.get()) > 0:

            new_data = np.stack([data for label, data in obj_data.get().items()])

            h5_file = h5py.File(filepath, mode = 'w-')

            if self.data_type is None:
                h5_file.create_dataset(name, data=new_data)
            else:
                h5_file.create_dataset(name, data=new_data, dtype=self.data_type)

            h5_file.close()
