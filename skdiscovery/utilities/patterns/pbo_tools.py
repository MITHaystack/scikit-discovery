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

# """@package pbo_tools
# Tools for working with PBO GPS data
# """

import time

import numpy as np
import pandas as pd

import skdiscovery.utilities.planetary.map_util as mo

def mogi(xdata, y, x, source_depth, amplitude, latlon=True):
    '''
    Compute the surface deformation due to changes in a mogi source

    @param xdata: List of the position data with each array element containing [ direction (x, y, or z), lat, lon ]
    @param y: Source y Position of (default: latitude)
    @param x: Source x Position (default longitude)
    @param source_depth: Depth of source
    @param amplitude: Amplitude of mogi source
    @param latlon: Source y is latitude and source x is longitude

    @return list of resulting deformation for each point in xdata
    '''
    source_coords = (y, x)

    results = []

    for data in xdata:

        dim = data[0]
        station_coords = (float(data[1]),float(data[2]))
        # print(station_coords)

        if latlon==True:
            y_distance = mo.wgs84_distance( source_coords, (station_coords[0], source_coords[1]) )
            x_distance = mo.wgs84_distance( source_coords, (source_coords[0], station_coords[1]) )
            x_distance = x_distance * np.sign(station_coords[1] - source_coords[1])
            y_distance = y_distance * np.sign(station_coords[0] - source_coords[0])
        else:
            y_distance = station_coords[0] - source_coords[0]
            x_distance = station_coords[1] - source_coords[1]

        R3 = (x_distance**2 + y_distance**2 + source_depth**2)**(3/2)

        result = None

        if dim == 'x':
            result = amplitude * x_distance / R3
        elif dim == 'y':
            result = amplitude * y_distance / R3
        elif dim == 'z':
            result = amplitude * source_depth / R3
        else:
            print("Did not understand dimension")

        results.append(result)
    return results



def finite_sphere(xdata, lat, lon, source_depth, amplitude, alpha_rad):
    '''
    Compute the surface deformation due to changes in a finite sphere source

    For reference, see "Volcano Deformation", Dzurisin 2006, pg 290
    (http://link.springer.com/book/10.1007/978-3-540-49302-0)

    @param xdata: List of the position data with each array element containing [ direction (x, y, or z), lat, lon ]
    @param lat: Latitude of source
    @param lon: Longitude of source
    @param source_depth: Depth of source
    @param amplitude: Ampltiude of source
    @param alpha_rad: Alpha radius of the source

    @return list of resulting deformation for each point in xdata

    '''
    nu_v = .25
    C1 = (1+nu_v)/(2*(-7+5*nu_v))
    C2 = 15*(-2+nu_v)/(4*(-7+5*nu_v))
    source_coords = (lat, lon)

    results = []

    for data in xdata:

        dim = data[0]
        station_coords = (float(data[1]),float(data[2]))
        # print(station_coords)

        y_distance = mo.wgs84_distance( source_coords, (station_coords[0], source_coords[1]) )
        x_distance = mo.wgs84_distance( source_coords, (source_coords[0], station_coords[1]) )
        x_distance = x_distance * np.sign(station_coords[1] - source_coords[1])
        y_distance = y_distance * np.sign(station_coords[0] - source_coords[0])

        R3 = (x_distance**2 + y_distance**2 + source_depth**2)**(3/2)
        result = None
        if dim == 'x':
            result = amplitude *alpha_rad**3*(1+(alpha_rad/source_depth)**3*(C1+C2*source_depth**2/R3**(2/3))) * x_distance / R3
        elif dim == 'y':
            result = amplitude *alpha_rad**3*(1+(alpha_rad/source_depth)**3*(C1+C2*source_depth**2/R3**(2/3))) * y_distance / R3
        elif dim == 'z':
            result = amplitude *alpha_rad**3*(1+(alpha_rad/source_depth)**3*(C1+C2*source_depth**2/R3**(2/3))) * source_depth / R3
        else:
            print("Did not understand dimension")

        results.append(result)
    return results


def closed_pipe(xdata, lat, lon, source_depth, amplitude, pipe_delta):
    '''
    Compute the surface deformation due to changes in a closed pipe source

    For reference, see "Volcano Deformation", Dzurisin 2006, pg 292
    (http://link.springer.com/book/10.1007/978-3-540-49302-0)

    @param xdata: List of the position data with each array element containing [ direction (x, y, or z), lat, lon ]
    @param lat: Latitude of source
    @param lon: Longitude of source
    @param source_depth: Depth of source
    @param amplitude: Ampltiude of source
    @param pipe_delta: Pipe delta from source depth to top/bottom

    @return list of resulting deformation for each point in xdata
    '''
    nu_v = .25
    source_coords = (lat, lon)

    results = []

    for data in xdata:

        dim = data[0]
        station_coords = (float(data[1]),float(data[2]))
        # print(station_coords)

        y_distance = mo.wgs84_distance( source_coords, (station_coords[0], source_coords[1]) )
        x_distance = mo.wgs84_distance( source_coords, (source_coords[0], station_coords[1]) )
        x_distance = x_distance * np.sign(station_coords[1] - source_coords[1])
        y_distance = y_distance * np.sign(station_coords[0] - source_coords[0])

        result = None
        c1 = source_depth - pipe_delta
        c2 = source_depth + pipe_delta
        R_1 = (x_distance**2 + y_distance**2 + c1**2)**(1/2)
        R_2 = (x_distance**2 + y_distance**2 + c2**2)**(1/2)
        r2  = (x_distance**2 + y_distance**2)
        if dim == 'x':
            result = amplitude *((c1/R_1)**3+2*c1*(-3+5*nu_v)/R_1+(5*c2**3*(1-2*nu_v)-2*c2*r2*(-3+5*nu_v))/R_2**3) * x_distance / r2
        elif dim == 'y':
            result = amplitude *((c1/R_1)**3+2*c1*(-3+5*nu_v)/R_1+(5*c2**3*(1-2*nu_v)-2*c2*r2*(-3+5*nu_v))/R_2**3) * y_distance / r2
        elif dim == 'z':
            result = - amplitude *(c1**2/R_1**3+2*(-2+5*nu_v)/R_1+(c2**2*(3-10*nu_v)-2*r2*(-2+5*nu_v))/R_2**3)
        else:
            print("Did not understand dimension")

        results.append(result)
    return results


def constant_open_pipe(xdata, lat, lon, source_depth, amplitude, pipe_delta):
    '''
    Compute the surface deformation due to changes in a constant width open pipe source

    For reference, see "Volcano Deformation", Dzurisin 2006, pg 295
    (http://link.springer.com/book/10.1007/978-3-540-49302-0)

    @param xdata: List of the position data with each array element containing [ direction (x, y, or z), lat, lon ]
    @param lat: Latitude of source
    @param lon: Longitude of source
    @param source_depth: Depth of source
    @param amplitude: Ampltiude of source
    @param pipe_delta: Pipe delta from source depth to top/bottom

    @return list of resulting deformation for each point in xdata
    '''
    nu_v = .25
    source_coords = (lat, lon)

    results = []

    for data in xdata:

        dim = data[0]
        station_coords = (float(data[1]),float(data[2]))
        # print(station_coords)

        y_distance = mo.wgs84_distance( source_coords, (station_coords[0], source_coords[1]) )
        x_distance = mo.wgs84_distance( source_coords, (source_coords[0], station_coords[1]) )
        x_distance = x_distance * np.sign(station_coords[1] - source_coords[1])
        y_distance = y_distance * np.sign(station_coords[0] - source_coords[0])

        result = None
        c1 = source_depth + pipe_delta
        c2 = source_depth - pipe_delta
        R_1 = (x_distance**2 + y_distance**2 + c1**2)**(1/2)
        R_2 = (x_distance**2 + y_distance**2 + c2**2)**(1/2)
        r2  = (x_distance**2 + y_distance**2)
        if dim == 'x':
            result = amplitude *((c1/R_1)**3-2*c1*(1+nu_v)/R_1+(c2**3*(1+2*nu_v)+2*c2*r2*(1+nu_v))/R_2**3)* x_distance / r2
        elif dim == 'y':
            result = amplitude *((c1/R_1)**3-2*c1*(1+nu_v)/R_1+(c2**3*(1+2*nu_v)+2*c2*r2*(1+nu_v))/R_2**3)* y_distance / r2
        elif dim == 'z':
            result = - amplitude *(c1**2/R_1**3-2*nu_v/R_1+(-c2**2+2*R_2**2*nu_v)/R_2**3)
        else:
            print("Did not understand dimension")

        results.append(result)
    return results


def rising_open_pipe(xdata, lat, lon, source_depth, amplitude, pipe_delta,open_pipe_top):
    '''
    Compute the surface deformation due to changes in a rising width amplitude open pipe source

    For reference, see "Volcano Deformation", Dzurisin 2006, pg 295
    (http://link.springer.com/book/10.1007/978-3-540-49302-0)

    @param xdata: List of the position data with each array element containing [ direction (x, y, or z), lat, lon ]
    @param lat: Latitude of source
    @param lon: Longitude of source
    @param source_depth: Depth of source
    @param amplitude: Ampltiude of source
    @param pipe_delta: Pipe delta from source depth to top/bottom
    @param open_pipe_top: Depth of the top of the open pipe

    @return list of resulting deformation for each point in xdata
    '''
    nu_v = .25
    source_coords = (lat, lon)

    results = []

    for data in xdata:

        dim = data[0]
        station_coords = (float(data[1]),float(data[2]))
        # print(station_coords)

        y_distance = mo.wgs84_distance( source_coords, (station_coords[0], source_coords[1]) )
        x_distance = mo.wgs84_distance( source_coords, (source_coords[0], station_coords[1]) )
        x_distance = x_distance * np.sign(station_coords[1] - source_coords[1])
        y_distance = y_distance * np.sign(station_coords[0] - source_coords[0])

        result = None
        c0 = open_pipe_top
        c1 = source_depth + pipe_delta
        R_0 = (x_distance**2 + y_distance**2 + c0**2)**(1/2)
        R_1 = (x_distance**2 + y_distance**2 + c1**2)**(1/2)
        r2  = (x_distance**2 + y_distance**2)
        if dim == 'x':
            result = amplitude *(-(c0**2/R_0**3)+2*nu_v/R_0+(c1**2-2*(c1**2+r2)*nu_v)/R_1**3)* x_distance / c1
        elif dim == 'y':
            result = amplitude *(-(c0**2/R_0**3)+2*nu_v/R_0+(c1**2-2*(c1**2+r2)*nu_v)/R_1**3)* y_distance / c1
        elif dim == 'z':
            result = -amplitude *((c0**3/R_0**3)-c1**3/R_1**3+c1*(-1+2*nu_v)/R_1+c0*(1-2*nu_v)/R_0+(-1+2*nu_v)*np.log(c0+R_0)-(-1+2*nu_v)*np.log(c1+R_1))/ c1
        else:
            print("Did not understand dimension")

        results.append(result)
    return results


def sill(xdata, lat, lon, source_depth, amplitude):
    '''
    Compute the surface deformation due to changes in a sill-like source

    For reference, see "Volcano Deformation", Dzurisin 2006, pg 297
    (http://link.springer.com/book/10.1007/978-3-540-49302-0)

    @param xdata: List of the position data with each array element containing [ direction (x, y, or z), lat, lon ]
    @param lat: Latitude of source
    @param lon: Longitude of source
    @param source_depth: Depth of source
    @param amplitude: Ampltiude of source

    @return list of resulting deformation for each point in xdata
    '''
    source_coords = (lat, lon)
    results = []

    for data in xdata:

        dim = data[0]
        station_coords = (float(data[1]),float(data[2]))
        # print(station_coords)

        y_distance = mo.wgs84_distance( source_coords, (station_coords[0], source_coords[1]) )
        x_distance = mo.wgs84_distance( source_coords, (source_coords[0], station_coords[1]) )
        x_distance = x_distance * np.sign(station_coords[1] - source_coords[1])
        y_distance = y_distance * np.sign(station_coords[0] - source_coords[0])

        R5 = (x_distance**2 + y_distance**2 + source_depth**2)**(5/2)
        result = None
        if dim == 'x':
            result = amplitude * x_distance * source_depth**2 / R5
        elif dim == 'y':
            result = amplitude * y_distance * source_depth**2 / R5
        elif dim == 'z':
            result = amplitude * source_depth**3 / R5
        else:
            print("Did not understand dimension")

        results.append(result)
    return results


def dirEigenvectors(coord_list, pca_comps,pdir='H'):
    '''
    Takes eigenvectors (north and east) and forces them to point "outward"

    Flips the sign of the projection if needed so that eigenvectors point
    outward. Needed because the "positive" direction for PCA is arbitrary

    @param coord_list: Location of stations for projecting the eigenvectors
    @param pca_comps: PCA components
    @param pdir: PCA direction, vertical or horizontal

    @return station_lat_list: the station latitude coordinates
    @return station_lon_list: the station longitude coordinates
    @return ev_lat_list: the properly origented corresponding eigenvector latitude component
    @return ev_lon_list: the properly origented corresponding eigenvector longitude component
    @return direction scale factor (1 for no flip, or -1 for flip)
    '''
    station_lat_list = []
    station_lon_list = []
    ev_lat_list = []
    ev_lon_list = []
    if pdir=='H':
        for i in range(int(len(pca_comps)/2)):
            station_lat_list.append(coord_list[i][0])
            station_lon_list.append(coord_list[i][1])
            ev_lat_list.append(pca_comps[2*i])
            ev_lon_list.append(pca_comps[2*i+1])

    elif pdir=='V':
        for i in range(int(len(pca_comps))):
            station_lat_list.append(coord_list[i][0])
            station_lon_list.append(coord_list[i][1])
            ev_lat_list.append(pca_comps[i])
            ev_lon_list.append(0)

    station_lat_list = np.array(station_lat_list)
    station_lon_list = np.array(station_lon_list)
    ev_lat_list = np.array(ev_lat_list)
    ev_lon_list = np.array(ev_lon_list)

    if pdir=='H':
        # mean coordinate
        mean_lon = np.mean(station_lon_list)
        mean_lat = np.mean(station_lat_list)

        # project the eigen vector along the direction relative to mean
        mean_angle = np.arctan2(station_lat_list-mean_lat,station_lon_list-mean_lon)
        eigv_angle = np.arctan2(ev_lat_list,ev_lon_list)
        proj_angle = np.abs(mean_angle - eigv_angle)

        if np.sum(proj_angle<np.pi/2)<np.sum(proj_angle>=np.pi/2):
            return station_lat_list, station_lon_list, ev_lat_list*-1, ev_lon_list*-1, -1
        else:
            return station_lat_list, station_lon_list, ev_lat_list, ev_lon_list, 1
    else:
        return station_lat_list, station_lon_list, ev_lat_list, ev_lon_list, 1


def datetimeToNumber(in_time):
    '''
    Converts input pandas Timestamp or pandas DatetimeIndex to unix time

    @param in_time: Input pandas timestamp or pandas DatetimeIndex

    @return unix time
    '''
    if isinstance (in_time, pd.Timestamp):
        return time.mktime(in_time.timetuple())
    elif isinstance(in_time, pd.DatetimeIndex):
        return in_time.map(lambda t: time.mktime(t.timetuple()))
    else:
        return in_time
