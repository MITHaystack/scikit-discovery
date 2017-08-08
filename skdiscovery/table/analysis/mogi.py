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

from skdiscovery.framework import PipelineItem
import skdiscovery.utilities.pbo_tools as pbo_tools
import skdaccess.utilities.pbo_util as pbo_utils
import numpy as np
import scipy.optimize as optimize
import collections


class Mogi_Inversion(PipelineItem):
    '''
    Perform a mogi source inversion on a set of gps table data

    The source is assumed to be a mogi source (point source), but other source models can be selected.
    Assumes directions are named ('dN', 'dE', 'dU').
    '''
    def __init__(self, str_description, ap_paramList, pca_name, column_names=['dN', 'dE', 'dU']):
        '''
        Initialize Mogi analysis item

        @param str_description: Description of item
        @param ap_paramList[source_type]: Type of magma chamber source model to use (default-mogi,finite_sphere,closed_pipe,constant_open_pipe,rising_open_pipe,sill)
        @param pca_name: Name of pca result
        @param column_names: The data direction column names
        ''' 

        self.pca_name = pca_name
        self.column_names = column_names
        super(Mogi_Inversion, self).__init__(str_description, ap_paramList)
        self.ap_paramNames = ['source_type']
        
        
    def FitPCA(self, hPCA_Proj):
        '''
        Determine the timing of the inflation event from the first component of the pca projection

        fits A * arctan( (t - t0) / c ) + B to the first pca projection, in order to estimate
        source amplitude parameters

        @param hPCA_Proj: The sklearn PCA 

        @return ct: the t0, c, and B parameters from the fit
        @return pA[0]: the fitted amplitude parameter
        '''

        fitfunc = lambda p,t: p[0]*np.arctan((t-p[1])/p[2])+p[3]
        errfunc = lambda p,x,y: fitfunc(p,x) - y
    
        dLen = len(hPCA_Proj[:,0])
        pA, success = optimize.leastsq(errfunc,[1.,dLen/2.,1.,0.],args=(np.arange(dLen),hPCA_Proj[:,0]))
        ct = pA[1:3]

        return ct, pA[0]
        

    def FitTimeSeries(self, pd_series, ct):
        '''
        Fits the amplitude and offset of an inflation event given the time and length of the event.
        Fits A and B in A * arctan( (t - t0) / c) + B
        
        @param pd_series: Time series to be fit
        @param ct: the time constants for the arctan

        @return res: Amplitude of the fit
        @return perr_leastsq: Error of the fit
        '''

        fitfunc2 = lambda p,c,t: p[0]*np.arctan((t-c[0])/c[1])+p[1]
        errfunc2 = lambda p,c,x,y: fitfunc2(p,c,x) - y

        dLen = len(pd_series)
        
        pA, pcov = optimize.leastsq(errfunc2,[1.,0.],args=(ct,np.arange(dLen),pd_series))
        # res = fitfunc2(pA,ct,np.arange(dLen))[-1]-fitfunc2(pA,ct,np.arange(dLen))[0]
        res = pA[0]*np.pi
        
        s_sq = (errfunc2(pA,ct,np.arange(dLen),pd_series)**2).sum()/(len(pd_series)-2)
        pcov = pcov * s_sq
        error = [] 
        for i in range(len(pA)):
            try:
              error.append(np.absolute(pcov[i][i])**0.5)
            except:
              error.append( 0.00 )
        perr_leastsq = np.array(error) 

        return res, perr_leastsq


    def process(self, obj_data):
        '''
        Finds the magma source (default-mogi) from PBO GPS data.

        Assumes time series columns are named ('dN', 'dE', 'dU'). Predicts the
        location of the magma source using scipy.optimize.curve_fit

        The result is added to the data wrapper as a list, with the four
        elements describing the location of the magma source:
        res[0] = latitude
        res[1] = longitude
        res[2] = source depth (km)
        res[3] = volume change (meters^3)

        @param obj_data: 
        '''
        h_pca_name = self.pca_name
        
        exN = {'mogi':0,'finite_sphere':1,'closed_pipe':1,'constant_open_pipe':1,'rising_open_pipe':2,'sill':0}
        try:
            mag_source = getattr(pbo_tools,self.ap_paramList[0]().lower())
            ExScParams = tuple(np.ones((exN[self.ap_paramList[0]().lower()],)))
        except:
            mag_source = pbo_tools.mogi
            ExScParams = ()
            print('No source type called '+self.ap_paramList[0]()+', defaulting to a Mogi source.')
                
        projection = obj_data.getResults()[h_pca_name]['Projection']
        start_date = obj_data.getResults()[h_pca_name]['start_date']
        end_date = obj_data.getResults()[h_pca_name]['end_date']        

        ct, pca_amp = self.FitPCA(projection)
        pca_amp *= np.pi

        xvs = []
        yvs = []
        zvs = []

        label_list = []

        for label, data in obj_data.getIterator():
            label_list.append(label)
            for column in self.column_names:
                distance,f_error = self.FitTimeSeries(data.loc[start_date:end_date,column], ct)
                if column == self.column_names[1]:
                    xvs.append(distance)
                elif column == self.column_names[0]:
                    yvs.append(distance)
                elif column == self.column_names[2]:
                    zvs.append(distance)
                else:
                    print('Ignoring column: ', label)

        xvs = np.array(xvs)*1e-6
        yvs = np.array(yvs)*1e-6
        zvs = np.array(zvs)*1e-6

        ydata = np.hstack((xvs, yvs,zvs)).T
        station_list = obj_data.get().keys()
        meta_data = obj_data.info()
        station_coords = pbo_utils.getStationCoords(meta_data, station_list)
        
        dimensions = ('x','y','z')
        xdata = []
        for dim in dimensions:
            for coord in station_coords:
                xdata.append((dim, coord[0], coord[1]))

        coord_range = np.array(pbo_utils.getLatLonRange(meta_data, station_list))

        lat_guess = np.mean(coord_range[0,:])
        lon_guess = np.mean(coord_range[1,:])

        fit = optimize.curve_fit(mag_source, xdata, ydata, (lat_guess, lon_guess, 5, 1e-4)+ExScParams)

        res = collections.OrderedDict()

        res['lat'] = fit[0][0]
        res['lon'] = fit[0][1]
        res['depth'] = fit[0][2]
        res['amplitude'] = fit[0][3]
        res['labels'] = label_list
        if len(fit[0])>4:
            res['ex_params'] = fit[0][4:]
        else:
            res['ex_params'] = np.nan
        res['pca_amplitude'] = pca_amp
        res['source_type'] = self.ap_paramList[0]().lower()

        obj_data.addResult(self.str_description, res)

        # lat_fit_range = (np.min(lat_list)-0.15, np.max(lat_list)+0.15)
        # lon_fit_range = (np.min(lon_list)-0.15, np.max(lon_list)+0.15)

        

        # res = optimize.brute(self.mogi, (lat_fit_range, lon_fit_range,
        #                                  (1,10), (1e-5, 1e-3)),
        #                      args = (xvs*1e-6, yvs*1e-6, zvs*1e-6,
        #                              station_list, meta_data))


def MogiVectors(mogi_res,station_lat_list,station_lon_list,flag3D=False):
    '''
    Creates a set of mogi vectors for plotting

    @param mogi_res: Magma source
    @param station_lat_list: List of station latitudes
    @param station_lon_list: List of station longitudes
    @param flag3D: Make vectors 3 dimensional, defaults to False (2D)
    '''
    
    # grab the correct magma forward function
    mag_source = getattr(pbo_tools,mogi_res['source_type'])

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
