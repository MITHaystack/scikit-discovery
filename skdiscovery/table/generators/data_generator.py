from skdaccess.framework.data_class import DataFetcherBase
from skdaccess.framework.data_class import TableWrapper
import pandas as pd
import numpy as np

class DataGenerator(DataFetcherBase):

    ''' In Development: Class for generating random data '''
    
    def __init__(self, length, *args, seed = None, final_function = None):
        '''
        Initialize Random data generator


        @param length: Number of rows to generate
        @param *args: Dictionaries containing entries: 'name',,'start', 'end', and optionally 'func'
        @param seed: Seed to use when generating random data
        @final_function: Final function to call on random data
        '''

        self.length = length
        self.seed = seed
        self.args = args
        self.final_function = final_function

    def output(self):
        if self.seed is not None:
            np.random.seed(self.seed)
            

        new_data = dict()
        name_list = []
        for arg in self.args:
            new_data[arg['name']] = np.random.rand(self.length) * (arg['end'] - arg['start']) + arg['start']

            name_list.append(arg['name'])

            if 'func' in arg:
                new_data[arg['name']] = arg['func'](new_data[arg['name']])


        if self.final_function is not None:
            new_data = pd.DataFrame.from_dict(new_data)
            new_data, updated_column_names = self.final_function(new_data)
            if updated_column_names is not None:
                default_columns = updated_column_names


        data = {'generated_data' : new_data}
        data_wrapper = TableWrapper(data, default_columns = name_list)
        
        return data_wrapper
