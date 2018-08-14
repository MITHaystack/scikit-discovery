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

# Standad library imports
import itertools

# Scikit Discovery imports
from skdiscovery.data_structure.framework.base import PipelineItem

# Third party imports
import numpy as np

class AddRandomTile(PipelineItem):
    """
    Add a randomly selected tile to daa
    """

    def __init__(self, str_description, tiles):
        """
        Initialize Item

        @param str_description: String description of item
        @param tiles: Numpy array here the first dimension is the one randomly selected from
        """

        self.tiles = tiles
        super(AddRandomTile, self).__init__(str_description)

    def process(self, obj_data):
        """
        Add a random tile to the data

        @param obj_data: Image data wrapper
        """

        tile_indices = np.arange(self.tiles.shape[0])
        np.random.shuffle(tile_indices)

        for tile_index, (label, data) in zip(itertools.cycle(tile_indices), obj_data.getIterator()):

            if data.ndim == 2:
                combine_function = np.stack
            elif data.ndim == 3:
                combine_funciton = np.concatenate
            else:
                raise RuntimeError('Only 2 and 3d images supported')


            new_data = combine_function([data, self.tiles[tile_index,:,:]])

            obj_data.updateData(label, new_data)
