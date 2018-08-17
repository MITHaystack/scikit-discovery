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

# skdiscovery imports
from skdiscovery.data_structure.framework.base import PipelineItem

# 3rd party imports
import numpy as np

class Apply(PipelineItem):
    """
    Apply a function to data on
    """

    def __init__(self, str_description, function):
        """
        Initialize Apply object

        @param str_description: String describing item
        @param function: Function to apply to data
        """

        self.function = function
        super(Apply, self).__init__(str_description)

    def process(self, obj_data):
        """
        Apply a function to data in an Image Wrapper

        @param obj_data: Image Wrapper
        """
        for label, data in obj_data.getIterator():
            obj_data.updateData(label, self.function(data))

