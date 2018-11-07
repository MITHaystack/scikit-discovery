import numpy as np
import pandas as pd
from filterpy.kalman import KalmanFilter


from skdiscovery.data_structure.framework.base import TablePipelineItem
from skdiscovery.utilities.patterns import kalman_smoother


class RTS_Smoother(TablePipelineItem):
    """
    ** In development ** Initialize RTS Smoother

    @param str_description: String describing filter
    @param ap_paramList[ap_tau]: the correlation time
    @param ap_paramList[ap_sigmaSq]: the data noise
    @param ap_paramList[ap_R]: the process noise
    """


    def process(self, obj_data):
        """
        Apply RTS smoother to table data

        @pama obj_data: Table wrapper
        """

        ap_tau = self.ap_paramList[0]()
        ap_sigma_sq = self.ap_paramList[1]()
        ap_R = self.ap_paramList[2]()


        column_list = self._getColumns(obj_data)

        error_column_list = self._getErrorColumns(obj_data)

        for label, data in obj_data.getIterator():
            for column in column_list:

                column_data = data[column]

                # Get first non zero value:
                good_index = np.argmax(~np.isnan(column_data.as_matrix()))
                x0 = column_data.iloc[good_index]

                # Create kalman filter
                kf = KalmanFilter(dim_x=1, dim_z=1)
                kf.x = np.array([x0])
                kf.F = np.array([[np.exp((-1.0)/ap_tau)]])
                kf.H = np.array([[1.0]])
                kf.P *= 1e3
                kf.R = ap_R
                kf.Q = np.array([[ap_sigma_sq*(1-np.exp(-2/ap_tau))]])


                input_data = column_data.tolist()
                for nan_index in np.nonzero(np.isnan(column_data.as_matrix()))[0]:
                    input_data[nan_index] = None

                kf_results = kf.batch_filter(input_data)

                rts_results = kf.rts_smoother(*kf_results[:2])

                obj_data.updateData(label, data.index, column, rts_results[0].reshape(-1))
