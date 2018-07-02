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

from six.moves import configparser
import os

def getConfig():
    '''
    Retrieve skdiscovery configuaration

    @return skdiscovery configparser
    '''

    config_location = os.path.expanduser('~') + '/.skdiscovery.conf'

    conf = configparser.ConfigParser()
    conf.read(config_location)

    return conf

def writeConfigValue(section,key,value):
    '''
    Write config to disk

    @param section: Name of section
    @param key: Name of key
    @param value: Value to write
    '''

    conf = getConfig()

    if not conf.has_section(section):
        conf.add_section(section)

    conf.set(section, key, value)    
    
    config_location = os.path.expanduser('~') + '/.skdiscovery.conf'
    config_handle = open(config_location, "w")
    conf.write(config_handle)
    config_handle.close()


def getConfigValue(section, key):
    '''
    Retrieve a value from the config file

    @param section: Section of the configuration file that contains the value
    @param key: Key of the value

    @return value in the specified section associated with given key
    '''

    conf = getConfig()

    return conf.get(section, key, fallback=None)

def getDispyPassword():
    '''
    Get dispy password

    @return dispy password
    '''
    conf = getConfig()

    return conf.get('Dispy', 'password', fallback=None)

def getHostName():
    '''
    Get Host name for displaying link to dispy status

    @return Hostname
    '''
    conf = getConfig()

    return conf.get('Dispy','hostname')
