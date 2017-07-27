from collections import defaultdict

from skdiscovery.framework import PipelineItem
from scipy.stats import skew
import numpy as np

class Skew(PipelineItem):
    ''' Calculates the skew of table data '''
    def process(self, obj_data):
        '''
        Apply Skew analysis with results added to the data wrapper

        @param obj_data: Data wrapper
        '''

        column_names = obj_data.getDefaultColumns()
        
        results = defaultdict(dict)
        # for label, frame in tqdm(obj_data.getIterator()):
        for label, frame in obj_data.getIterator():
            for column in column_names:
                # dropping missing data in order to remove top and bottom 2%
                data = frame[column].dropna()
                # Remove top and bottom 2%
                rem_num = round(len(data)*0.02)
                res = skew(data.sort_values(ascending=True)[rem_num:-rem_num])
                if isinstance(res, np.ma.masked_array):
                    res = np.float(res.data)
                results[label][column] = res
                obj_data.addResult(self.str_description, results)
