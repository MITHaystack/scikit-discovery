# The MIT License (MIT)
# Copyright (c) 2017, 2018 Massachusetts Institute of Technology
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

class SourceWrapper(object):
    """
    Wrapper for using old interface with updated source interfaces
    """
    def __init__(self, source_method):
        """
        Initialize source wrapper

        @param source_method: Source function that will be wrapped
        """
        self.source_method = source_method

    def __call__(self, *args):
        """
        Call the source function using the old interface

        @param args: Arguments for the wrapped source function
        @return return list of resulting deformation for each point requested point
        """
        xdata = args[0]
        args = args[1:]

        results = []

        for data in xdata:
            direction = data[0]
            lat = float(data[1])
            lon = float(data[2])

            function_res = self.source_method(lat, lon, *args)

            if direction == 'x':
                results.append(function_res[0,0])
            elif direction == 'y':
                results.append(function_res[0,1])
            elif direction == 'z':
                results.append(function_res[0,2])

        return results

def getLength(position_y, position_x):
    """
    Get the length of the input position y and position x data

    @param position_y: y positions
    @param position_x: x positions

    @return The maximum length between the x and y positions
    """
    length_y = 1
    length_x = 1

    if not np.isscalar(position_y):
        length_y = len(position_y)

    if not np.isscalar(position_x):
        length_x = len(position_x)

    return max(length_y, length_x)


def compute_distances(position_y, position_x, source_y, source_x, latlon=True):
    """
    Compute the y and x distance between the observation location and the source location

    @param position_y: Obsevation y position
    @param position_x: Observation x position
    @param source_y: Source y position
    @param source_x: Source x position
    @param latlon: Interpret positions as latitudes and longitudes

    @return The y and x distance between observation location and source locaiton
    """
    if latlon==True:
        y_distance = mo.wgs84_distance( (source_y, source_x), (position_y, source_x) )
        x_distance = mo.wgs84_distance( (source_y, source_x), (source_y, position_x) )
        x_distance = x_distance * np.sign(position_x - source_x)
        y_distance = y_distance * np.sign(position_y - source_y)

    else:
        x_distance = position_x - source_x
        y_distance = position_y - source_y


    return y_distance, x_distance



def mogi(position_y, position_x, source_y, source_x, source_depth, amplitude, latlon=True):
    '''
    Compute the surface deformation due to changes in a mogi source

    @param position_y: Observation positions in the y coordinate
    @param position_x: Observation positions in the x coordinate
    @param source_y: Position of the source in the y coordinate
    @param source_x: Position of the source in the x coordinate
    @param source_depth: Depth of source
    @param amplitude: Amplitude of mogi source
    @param latlon: If true, then position_y, position_x, source_y, and source_x
                   are given in latitude and longitude coordinates

    @return Array containing the x, y, and z deformations
    '''

    y_distance, x_distance = compute_distances(position_y, position_x, source_y, source_x, latlon)

    R3 = (x_distance**2 + y_distance**2 + source_depth**2)**(3/2)

    result_x = amplitude * x_distance / R3
    result_y = amplitude * y_distance / R3
    result_z = amplitude * source_depth / R3

    return np.column_stack([result_x, result_y, result_z])


def finite_sphere(position_y, position_x, source_y, source_x, source_depth, amplitude, alpha_rad, latlon=True):
    '''
    Compute the surface deformation due to changes in a finite sphere source

    For reference, see "Volcano Deformation", Dzurisin 2006, pg 290
    (http://link.springer.com/book/10.1007/978-3-540-49302-0)

    @param position_y: Observation positions in the y coordinate
    @param position_x: Observation positions in the x coordinate
    @param source_y: Position of the source in the y coordinate
    @param source_x: Position of the source in the x coordinate
    @param source_depth: Depth of source
    @param amplitude: Ampltiude of source
    @param alpha_rad: Alpha radius of the source
    @param latlon: If true, then position_y, position_x, source_y, and source_x
                   are given in latitude and longitude coordinates

    @return Array containing the x, y, and z deformations
    '''
    nu_v = .25
    C1 = (1+nu_v)/(2*(-7+5*nu_v))
    C2 = 15*(-2+nu_v)/(4*(-7+5*nu_v))

    y_distance, x_distance = compute_distances(position_y, position_x, source_y, source_x, latlon)


    R3 = (x_distance**2 + y_distance**2 + source_depth**2)**(3/2)

    result_x = amplitude *alpha_rad**3*(1+(alpha_rad/source_depth)**3*(C1+C2*source_depth**2/R3**(2/3))) * x_distance / R3
    result_y = amplitude *alpha_rad**3*(1+(alpha_rad/source_depth)**3*(C1+C2*source_depth**2/R3**(2/3))) * y_distance / R3
    result_z = amplitude *alpha_rad**3*(1+(alpha_rad/source_depth)**3*(C1+C2*source_depth**2/R3**(2/3))) * source_depth / R3

    return np.column_stack([result_x, result_y, result_z])


def closed_pipe(position_y, position_x, source_y, source_x, source_depth, amplitude, pipe_delta, latlon=True):
    '''
    Compute the surface deformation due to changes in a closed pipe source

    For reference, see "Volcano Deformation", Dzurisin 2006, pg 292
    (http://link.springer.com/book/10.1007/978-3-540-49302-0)

    @param position_y: Observation positions in the y coordinate
    @param position_x: Observation positions in the x coordinate
    @param source_y: Position of the source in the y coordinate
    @param source_x: Position of the source in the x coordinate
    @param source_depth: Depth of source
    @param amplitude: Ampltiude of source
    @param pipe_delta: Pipe delta from source depth to top/bottom
    @param latlon: If true, then position_y, position_x, source_y, and source_x
                   are given in latitude and longitude coordinates

    @return Array containing the x, y, and z deformations
    '''
    nu_v = .25

    y_distance, x_distance = compute_distances(position_y, position_x, source_y, source_x, latlon=latlon)

    c1 = source_depth - pipe_delta
    c2 = source_depth + pipe_delta
    R_1 = (x_distance**2 + y_distance**2 + c1**2)**(1/2)
    R_2 = (x_distance**2 + y_distance**2 + c2**2)**(1/2)
    r2  = (x_distance**2 + y_distance**2)


    result_x = amplitude *((c1/R_1)**3+2*c1*(-3+5*nu_v)/R_1+(5*c2**3*(1-2*nu_v)-2*c2*r2*(-3+5*nu_v))/R_2**3) * x_distance / r2
    result_y = amplitude *((c1/R_1)**3+2*c1*(-3+5*nu_v)/R_1+(5*c2**3*(1-2*nu_v)-2*c2*r2*(-3+5*nu_v))/R_2**3) * y_distance / r2
    result_z = -1*amplitude *(c1**2/R_1**3+2*(-2+5*nu_v)/R_1+(c2**2*(3-10*nu_v)-2*r2*(-2+5*nu_v))/R_2**3)

    return np.column_stack([result_x, result_y, result_z])


def constant_open_pipe(position_y, position_x, source_y, source_x, source_depth, amplitude, pipe_delta, latlon=True):
    '''
    Compute the surface deformation due to changes in a constant width open pipe source

    For reference, see "Volcano Deformation", Dzurisin 2006, pg 295
    (http://link.springer.com/book/10.1007/978-3-540-49302-0)

    @param position_y: Observation positions in the y coordinate
    @param position_x: Observation positions in the x coordinate
    @param source_y: Position of the source in the y coordinate
    @param source_x: Position of the source in the x coordinate
    @param source_depth: Depth of source
    @param amplitude: Ampltiude of source
    @param pipe_delta: Pipe delta from source depth to top/bottom
    @param latlon: If true, then position_y, position_x, source_y, and source_x
                   are given in latitude and longitude coordinates

    @return Array containing the x, y, and z deformations
    '''
    nu_v = .25
    C1 = (1+nu_v)/(2*(-7+5*nu_v))
    C2 = 15*(-2+nu_v)/(4*(-7+5*nu_v))


    result_length = getLength(position_y, position_x)
    result = np.zeros([result_length, 3])

    y_distance, x_distance = compute_distances(position_y, position_x, source_y, source_x, latlon=latlon)

    c1 = source_depth + pipe_delta
    c2 = source_depth - pipe_delta
    R_1 = (x_distance**2 + y_distance**2 + c1**2)**(1/2)
    R_2 = (x_distance**2 + y_distance**2 + c2**2)**(1/2)
    r2  = (x_distance**2 + y_distance**2)

    result[:,0] = amplitude *((c1/R_1)**3-2*c1*(1+nu_v)/R_1+(c2**3*(1+2*nu_v)+2*c2*r2*(1+nu_v))/R_2**3)* x_distance / r2
    result[:,1] = amplitude *((c1/R_1)**3-2*c1*(1+nu_v)/R_1+(c2**3*(1+2*nu_v)+2*c2*r2*(1+nu_v))/R_2**3)* y_distance / r2
    result[:,2] = - amplitude *(c1**2/R_1**3-2*nu_v/R_1+(-c2**2+2*R_2**2*nu_v)/R_2**3)

    # if result_length == 1:
    #     result = result.ravel()

    return result



def rising_open_pipe(position_y, position_x, source_y, source_x, source_depth, amplitude, pipe_delta, latlon=True):
    '''
    Compute the surface deformation due to changes in a rising width amplitude open pipe source

    For reference, see "Volcano Deformation", Dzurisin 2006, pg 295
    (http://link.springer.com/book/10.1007/978-3-540-49302-0)

    @param position_y: Observation positions in the y coordinate
    @param position_x: Observation positions in the x coordinate
    @param source_y: Position of the source in the y coordinate
    @param source_x: Position of the source in the x coordinate
    @param source_depth: Depth of source
    @param amplitude: Ampltiude of source
    @param pipe_delta: Pipe delta from source depth to top/bottom
    @param open_pipe_top: Depth of the top of the open pipe
    @param latlon: If true, then position_y, position_x, source_y, and source_x
                   are given in latitude and longitude coordinates

    @return Array containing the x, y, and z deformations
    '''
    nu_v = .25

    result_length = getLength(position_y, position_x)
    result = np.zeros([result_length, 3])


    y_distance, x_distance = compute_distances(position_y, position_x, source_y, source_x, latlon=latlon)

    c0 = source_depth - pipe_delta
    c1 = source_depth + pipe_delta
    R_0 = (x_distance**2 + y_distance**2 + c0**2)**(1/2)
    R_1 = (x_distance**2 + y_distance**2 + c1**2)**(1/2)
    r2  = (x_distance**2 + y_distance**2)

    result[:,0] = amplitude *(-(c0**2/R_0**3)+2*nu_v/R_0+(c1**2-2*(c1**2+r2)*nu_v)/R_1**3)* x_distance / c1
    result[:,1] = amplitude *(-(c0**2/R_0**3)+2*nu_v/R_0+(c1**2-2*(c1**2+r2)*nu_v)/R_1**3)* y_distance / c1
    result[:,2] = -amplitude *((c0**3/R_0**3)-c1**3/R_1**3+c1*(-1+2*nu_v)/R_1+c0*(1-2*nu_v)/R_0+(-1+2*nu_v)*np.log(c0+R_0)-(-1+2*nu_v)*np.log(c1+R_1))/ c1

    # if result_length == 1:
    #     result = result.ravel()

    return result


def sill(position_y, position_x, source_y, source_x, source_depth, amplitude, latlon=True):
    '''
    Compute the surface deformation due to changes in a sill-like source

    For reference, see "Volcano Deformation", Dzurisin 2006, pg 297
    (http://link.springer.com/book/10.1007/978-3-540-49302-0)

    @param position_y: Station y location
    @param position_x: Station x location
    @param source_y: y position of source
    @param source_x: x position of source
    @param source_depth: Depth of source
    @param amplitude: Ampltiude of source
    @param latlon: If true, then position_y, position_x, source_y, and source_x
                   are given in latitude and longitude coordinates

    @return Array containing the x, y, and z deformations
    '''
    y_distance, x_distance = compute_distances(position_y, position_x, source_y, source_x, latlon)

    R5 = (x_distance**2 + y_distance**2 + source_depth**2)**(5/2)

    result_x = amplitude * x_distance * source_depth**2 / R5

    result_y = amplitude * y_distance * source_depth**2 / R5

    result_z = amplitude * source_depth**3 / R5

    return np.column_stack([result_x, result_y, result_z])


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



def MogiVectors(mogi_res,station_lat_list,station_lon_list,flag3D=False):
    '''
    Creates a set of Mogi vectors for plotting

    @param mogi_res: Magma source inversion results
    @param station_lat_list: List of station latitudes
    @param station_lon_list: List of station longitudes
    @param flag3D: Flag for generating 3 dimensional vectors instead of only horizontal

    @return x and y Mogi vectors scaled by pca amplitude change
    '''
    # grab the correct magma forward function
    mag_source = SourceWrapper(globals()[mogi_res['source_type']])

    mogi_x_disp = []
    mogi_y_disp = []
    for lat, lon in zip(station_lat_list, station_lon_list):
        mogi_data = []

        mogi_data.append( ('x', lat, lon) )
        mogi_data.append( ('y', lat, lon) )
        if flag3D:
            mogi_data.append( ('z', lat, lon) )

        if np.isnan(mogi_res['ex_params']).any()==True:
            res = mag_source(mogi_data, mogi_res['lat'], mogi_res['lon'],
                             mogi_res['depth'], mogi_res['amplitude'])
        elif len(mogi_res['ex_params'])==1:
            res = mag_source(mogi_data, mogi_res['lat'], mogi_res['lon'],
                             mogi_res['depth'], mogi_res['amplitude'],mogi_res['ex_params'][0])
        elif len(mogi_res['ex_params'])==2:
            res = mag_source(mogi_data, mogi_res['lat'], mogi_res['lon'],
                             mogi_res['depth'], mogi_res['amplitude'],mogi_res['ex_params'][0],mogi_res['ex_params'][1])

        mogi_x_disp.append(res[0])
        mogi_y_disp.append(res[1])

    # the scaling here is because mogi distance vectors are in km
    # the factor converts to mm and scales by the pca change
    mogi_x_disp = np.array(mogi_x_disp) * 1e6 / (mogi_res['pca_amplitude'])
    mogi_y_disp = np.array(mogi_y_disp) * 1e6 / (mogi_res['pca_amplitude'])

    return mogi_x_disp, mogi_y_disp
