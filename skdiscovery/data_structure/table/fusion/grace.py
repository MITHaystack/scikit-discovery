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


from skdiscovery.data_structure.framework.base import PipelineItem
from skdiscovery.data_structure.framework import DiscoveryPipeline
from skdiscovery.data_structure.generic.accumulators import DataAccumulator
from skdiscovery.data_structure.table.filters import CalibrateGRACE, Resample, CalibrateGRACEMascon
from skdiscovery.data_structure.framework.stagecontainers import *
from skdaccess.framework.param_class import *
from skdaccess.geo.grace import DataFetcher as GDF
from skdaccess.geo.grace.mascon.cache import DataFetcher as MasconDF
from skdaccess.geo.gldas import DataFetcher as GLDASDF

import numpy as np


class GraceFusion(PipelineItem):
    '''
    Fuses GRACE equivelent water depth time series

    Works on table data (original data from http://grace.jpl.nasa.gov/data/get-data/monthly-mass-grids-land/)
    '''

    def __init__(self, str_description, ap_paramList, metadata, column_data_name = 'Grace', column_error_name = 'Grace_Uncertainty'):
        '''
        Initialize Grace Fusion item

        @param str_description: String describing item
        @param ap_paramList[gldas]: How to use of the global land data assimilation water model
        @param ap_paramList[mascons]: Boolean indicating if the mascon solution should be used
        @param ap_paramList[apply_scale_factor]: Boolean indicating if the scaling factors shoud be applied
        @param metadata: Metadata that contains lat,lon coordinates based on data labels
        @param column_data_name: Name of column for GRACE data
        @param column_error_name: Grace Uncertainty column name
        '''

        super(GraceFusion, self).__init__(str_description, ap_paramList)
        self.metadata = metadata.copy()
        self.column_data_name = column_data_name
        self.column_error_name = column_error_name
        # remove_sm_and_snow
        self._tileCache = None


    def process(self, obj_data):
        '''
        Adds columns for GRACE data and uncertainties

        @param obj_data: Input DataWrapper, will be modified in place
        '''

        # Only perform fusion if data exists

        gldas = self.ap_paramList[0]()
        use_mascons = self.ap_paramList[1]()
        apply_scale_factor = self.ap_paramList[2]()

        if obj_data.getLength() > 0:
            start_date = None
            end_date = None

            for label, data in obj_data.getIterator():
                try:
                    lat = self.metadata[label]['Lat']
                    lon = self.metadata[label]['Lon']
                except:
                    lat = self.metadata.loc[label,'Lat']
                    lon = self.metadata.loc[label,'Lon']

                locations = [(lat,lon)]

                if start_date == None:
                    start_date = data.index[0]
                    end_date = data.index[-1]
                else:
                    if start_date != data.index[0] \
                    or end_date != data.index[-1]:
                        raise RuntimeError("Varying starting and ending dates not supported")


                al_locations = AutoList(locations)
                al_locations_gldas = AutoList(locations)
                if use_mascons == False:
                    graceDF = GDF([al_locations], start_date, end_date)

                else:
                    graceDF = MasconDF([al_locations], start_date, end_date)

                gldasDF = GLDASDF([al_locations_gldas], start_date, end_date)


                def getData(datafetcher, pipe_type):

                    ac_data = DataAccumulator('Data',[])
                    sc_data = StageContainer(ac_data)

                    fl_grace = CalibrateGRACE('Calibrate', apply_scale_factor = apply_scale_factor)
                    sc_grace = StageContainer(fl_grace)

                    fl_mascon = CalibrateGRACEMascon('CalibrateMascon', apply_scale_factor = apply_scale_factor)
                    sc_mascon = StageContainer(fl_mascon)

                    fl_resample = Resample('Resample',start_date, end_date)
                    sc_resample = StageContainer(fl_resample)

                    if pipe_type == 'grace':
                        pipeline = DiscoveryPipeline(datafetcher, [sc_grace, sc_resample, sc_data])

                    elif pipe_type == 'mascon':
                        pipeline = DiscoveryPipeline(datafetcher, [sc_mascon, sc_resample, sc_data])

                    elif pipe_type == 'gldas':
                        pipeline = DiscoveryPipeline(datafetcher, [sc_resample, sc_data])

                    else:
                        raise RuntimeError('pipe_type: ' + str(pipe_type) + ' not understood')

                    pipeline.run(num_cores=1)

                    key = list(pipeline.getResults(0)['Data'].keys())[0]

                    return pipeline.getResults(0)['Data'][key]



                # Load GRACE data
                if use_mascons == False:
                    grace_data = getData(graceDF, 'grace')
                else:
                    grace_data = getData(graceDF, 'mascon')

                if gldas.lower() == 'off':
                # We are not removing sm and snow
                    obj_data.addColumn(label, self.column_data_name, grace_data['EWD'])
                    obj_data.addColumn(label, self.column_error_name, grace_data['EWD_Error'])

                elif gldas.lower() == 'remove':
                    # If we are removing sm and snow
                    gldas_data = getData(gldasDF, 'gldas')


                    grace = grace_data['Data']
                    gldas = gldas_data['Data']

                    grace_index = grace.index
                    grace.dropna(inplace=True)

                    # If no grace data available, no need to remove gldas
                    if len(grace) == 0:
                        continue

                    # Get matching start and end
                    start_grace = grace.index[0]
                    end_grace = grace.index[-1]
                    start_gldas = gldas.index[0]
                    end_gldas = gldas.index[-1]

                    start_month = np.max([start_gldas,start_grace]).strftime('%Y-%m')
                    end_month = np.min([end_gldas,end_grace]).strftime('%Y-%m')

                    # Convert gldas to a data frame
                    # and save index
                    # gldas = gldas.loc[:,:,'GLDAS'].copy()
                    gldas.loc[:,'Date'] = gldas.index

                    # Index GLDAS data by month
                    new_index = [date.strftime('%Y-%m') for date in gldas.index]
                    gldas.loc[:,'Month'] = new_index
                    gldas.set_index('Month',inplace=True)

                    # select only months that are also in GRACE
                    cut_gldas = gldas.loc[[date.strftime('%Y-%m') for date in grace.loc[start_month:end_month,:].index],:]

                    # index GLDAS data to GRACE dates
                    cut_gldas.loc[:,'Grace Index'] =  grace.loc[start_month:end_month,:].index
                    cut_gldas.set_index('Grace Index', inplace=True)


                    # Calculate distance between days
                    offset_days = cut_gldas.index - cut_gldas.loc[:,'Date']
                    offset_days = offset_days.apply(lambda t: t.days)
                    cut_gldas.loc[:,'Offset'] = offset_days


                    # Remove any data where the difference between gldas and grace are > 10 days
                    cut_gldas = cut_gldas[np.abs(cut_gldas.loc[:,'Offset']) < 10].copy()

                    # Select appropriate Grace Data
                    cut_grace = grace.loc[cut_gldas.index,:]

                    # Remove contribution of snow + sm to GRACE
                    cut_grace.loc[:,'Grace'] = cut_grace.loc[:,'Grace'] - cut_gldas['GLDAS']

                    # Now restore to original index, filling in with NaN's
                    grace = cut_grace.reindex(grace_index)

                    # index, place the result back into the
                    grace_data[key]['Data'] = grace

                    # All the snow and sm contribution has been removed,
                    # so the dictionary can now be returned
                    obj_data.addColumn(label, self.column_data_name, grace_data['EWD'])
                    obj_data.addColumn(label, self.column_error_name, grace_data['EWD_Error'])

                elif gldas.lower() == 'only':
                    obj_data.addColumn(label, self.column_data_name, ['EWD'])

                else:
                    raise ValueError('Did not understand gldas option: ' + gldas.lower())
