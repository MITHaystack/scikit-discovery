# The MIT License (MIT)
# Copyright (c) 2017 Massachusetts Institute of Technology
#
# Authors: Victor Pankratius, Justin Li, Cody Rude
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

from skdiscovery.data_structure.framework.base import PipelineItem

from skdiscovery.utilities.planetary.map_util import wgs84_distance


class GeoLocationFilter(PipelineItem):
    '''
    Removes objects not located in a specified region
    '''
    
    def __init__(self, str_description, ap_paramList):
        '''
        Initialize GeolocationFilter

        @param str_description: String describing filter
        @param ap_paramList[ap_lat]: Latitude coordinate
        @param ap_paramList[ap_lon]: Longitude coordinate
        @param ap_paramList[ap_radius]: cut objects whose distance from lat/lon 
                                        is greater than ap_radius
        '''
        super(GeoLocationFilter, self).__init__(str_description, ap_paramList)
        
    def process(self, obj_data): 
        ''' 
        Apply geolocation filter to data set
        
        @param obj_data: Table data wrapper
        '''
       
        lat = self.ap_paramList[0]()
        lon = self.ap_paramList[1]()
        radius = self.ap_paramList[2]()
        
        remove_list = []
        
        for label, data in obj_data.getIterator():
            data_lat = obj_data.info(label)['Lat']
            data_lon = obj_data.info(label)['Lon']
            
            if wgs84_distance((lat,lon),(data_lat,data_lon)) > radius:
                remove_list.append(label)
                
        obj_data.removeFrames(remove_list) 
