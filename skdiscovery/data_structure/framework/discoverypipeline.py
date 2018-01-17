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


import copy
import re
from contextlib import ExitStack
from multiprocessing import Pool
from multiprocessing import cpu_count

import graphviz as gv
from IPython.display import SVG
from IPython.display import display
from graphviz import Source
from tqdm import tqdm

from . import config


def _cluster_run(data_fetcher, stage_containers, run_id=-1, verbose = False):
    ''' 
    Run the pipeline

    @param data_fetcher: Data fetcher to use as a data source (from skdaccess)
    @param stage_containers: List of stage containers
    @param run_id: Run_id to use, for offloading
    @param verbose: Display the tqdm progress bar for this run

    @return results from pipeline run
    '''
    data_container = data_fetcher.output()
    data_container.run_id = run_id

    if verbose == False:
        iterator = stage_containers
    else:
        iterator = tqdm(stage_containers)
    
    for s in iterator:
        s.run(data_container)

    return data_container.getResults()

def _wrap_cluster(args):
    ''' Wrap cluster run 

    @param args: Arguments to pass to _cluster_run

    @return results from _cluster_run
    '''
    return _cluster_run(*args)


class DiscoveryPipeline:
    ''' Pipeline for running the analysis '''
    
    def __init__(self, data_fetcher, list_StageContainers):
        '''
        Initialize a new pipeline.
        
        @param data_fetcher: Data fetcher to use as a data source (from skdaccess)
        @param list_StageContainers: List of stage containers
        '''
        
        self.stage_containers = list_StageContainers
        self.data_fetcher   = data_fetcher
        self.stageConfigurationHistory = []
        self.RA_results = []
        self.__cluster = None
        self._dispy_http = None        

    def run(self, num_runs=1, perturb_data = False, num_cores = 1, amazon=False, verbose=False):
        '''
        Run the pipeline

        @param num_runs: Number of times to run the pipeline
        @param perturb_data: Boolean flag. If running the pipeline multiple times then
                            perturb the data instead of the pipeline
        @param num_cores: Number of cores on the local machine to use. Defaults
                          to 1 core. Use 0 to select the minimum between the
                          number of runs and cpu cores.
        @param amazon: Offload the pipeline on amazon
        @param verbose: Display the pipeline for each run
        '''


        # Check to make sure data fetcher supports amazon or multicore if requested
        if (amazon == True or num_cores != 1) and self.data_fetcher.multirun_enabled() == False:
            raise RuntimeError('Data Fetcher type does not support amazon or using multiple cores')

        # Run the job on Amazon
        if amazon == True:

            if self.__cluster == None:
                self._startCluster(config.getDispyPassword())

            self._createDispyLink()

            # Run Jobs on Amazon
            jobs = []

            for i in range(0,num_runs):
                if num_runs > 1:
                    run_index = i
                else:
                    run_index = -1

                if verbose:
                    self.plotPipelineInstance()
                job = self.__cluster.submit(copy.deepcopy(self.data_fetcher), copy.deepcopy(self.stage_containers), run_id = run_index)
                job.id = i
                jobs.append(job)
                # save metadata configuration for history
                self.stageConfigurationHistory.append(self.getMetadata())

                if num_runs > 1:
                    if perturb_data == False:
                        self.perturb()
                    else:
                        self.perturbData()


            for job in jobs:
                results = job()
                # save results in the result accumulator
                if(job.exception is not None):
                    print(job.exception)
                self.RA_results.append(results)

            self.__cluster.stats()
            self.__cluster.close()
            self.__cluster = None

        # Run the job on the local machine
        else:

            # Function to generate inputs for running
            def generatePipelineInputs():
                self.stageConfigurationHistory.append(self.getMetadata())
                if verbose:
                    self.plotPipelineInstance()
                yield copy.deepcopy(self.data_fetcher), copy.deepcopy(self.stage_containers), 0

                for i in range(1, num_runs):
                    if perturb_data == False:
                        self.perturb()
                    else:
                        self.perturbData()

                    self.stageConfigurationHistory.append(self.getMetadata()) 
                    if verbose:
                        self.plotPipelineInstance()
                    yield copy.deepcopy(self.data_fetcher), copy.deepcopy(self.stage_containers), i

                # Perturb data as many times as runing more than once
                if num_runs > 1:
                    if perturb_data == False:
                        self.perturb()
                    else:
                        self.perturbData()
                    
            # Run the jobs
            with ExitStack() as stack:

                # Create workers for processing multiple pipelines locally
                if num_runs != 1:
                    if num_cores == 0:
                        # Create a pool with either the number of cores on
                        # the machine, or the number of runs, whichever is
                        # less
                        pool = Pool(min(cpu_count(), num_runs))
                        stack.callback(pool.close)
                    elif num_cores > 1:
                        # If number of cores is specified generate that many
                        # workers
                        pool = Pool(num_cores)
                        stack.callback(pool.close)

                if (num_cores == 0 or num_cores > 1) and num_runs != 1:
                    results = pool.map(_wrap_cluster, generatePipelineInputs())
                else:
                    results = list(map(_wrap_cluster, generatePipelineInputs()))

                self.RA_results += results                
        
    def perturb(self):
        ''' Perturb the paramters in the stage containers '''
        for s in self.stage_containers:
            s.perturb()

    def reset(self):
        ''' 
        Reset the stage containers to their default values and clear previous runs.
        '''
        for s in self.stage_containers:
            s.reset()

        self.stageConfigurationHistory = []
        self.RA_results = []

    def getMetadata(self):
        '''
        Retrieve Metadata from stage containers

        @return list of metadata for the current run
        '''

        data_metadata = self.data_fetcher.getMetadata()
        
        metadataList = [s.getMetadata() for s in self.stage_containers]
        returnList = [data_metadata]
        for info  in metadataList:
            if isinstance(info, list):
                returnList += info
            else:
                returnList.append(info)
                
        return returnList

    def getMetadataHistory(self):
        ''' 
        Get the metadata for each run in the pipeline 

        @return list of metadata configurations for all runs
        '''
        return self.stageConfigurationHistory

    def perturbData(self):
        ''' Perturb the input data '''
        self.data_fetcher.perturb()

    def getResults(self,index=None):
        '''
        Return results from previous runs.

        @param index: Index of run. If None, return all previous results

        @return results from a run at index. If index=None, returns list of all results
        '''
        if index==None:
            return self.RA_results
        else:
            return self.RA_results[index]

    def resultIter(self):
        ''' 
        Retrieves and iterator to the results and history of the pipeline.

        @return A 2 component iterator to the results and history of previous runs
        '''
        return zip(self.getResults(), self.getMetadataHistory())

    def plotPipelineInstance(self):
        '''
        Plot current instance of pipeline stages with metadata

        @return iPython display object
        '''
            
        #Graph setup
        g1 = gv.Digraph(format='svg')
        g1.graph_attr['rankdir'] = 'LR'
        g1.node_attr['shape'] = 'rounded'
        g1.node_attr['fontname'] = 'Arial'
        g1.node_attr['fontsize'] = '9'
        g1.node_attr['style'] = 'filled'
        g1.node_attr['margin'] = '0.1'
        g1.node_attr['height'] = '0.1'
        g1.node_attr['fillcolor'] = '#d8e9fd'
        #g1.node_attr['shape'] = 'plaintext'  #use this to remove boxes around nodes

        # Some stagecontainers directly return the metadata as string,
        # others return a list of strings for each item in the stage
        # container.
        metadata_list = []
        for s in self.stage_containers:
            first_metadata = s.getMetadata()
            if isinstance(first_metadata, str):
                metadata_list.append(first_metadata)

            else:
                for second_metadata in first_metadata:
                    metadata_list.append(second_metadata)
        
        nodelist = [re.sub('\:',';', metadata) for metadata in metadata_list]

        for s in nodelist:
            g1.node(s)
        
        g1.edges(zip(nodelist, nodelist[1:]))
        
        g1.render('img/plotPipelineInstance')
        
        #print(g1.source)
        #return display(Image('img/pipelinePlot.svg'))
        return display(SVG('img/plotPipelineInstance.svg'))

    def plotPipelineStructure(self):
        '''
        Plot pipeline structure

        @return iPython display object
        '''
            
        #Graph setup
        g1 = gv.Digraph(format='svg')
        g1.graph_attr['rankdir'] = 'LR'
        g1.node_attr['shape'] = 'rounded'
        g1.node_attr['fontname'] = 'Arial'
        g1.node_attr['fontsize'] = '9'
        g1.node_attr['style'] = 'filled'
        g1.node_attr['margin'] = '0.1'
        g1.node_attr['height'] = '0.1'
        g1.node_attr['fillcolor'] = '#fff7da'
        #g1.node_attr['shape'] = 'plaintext'  #use this to remove boxes around nodes
        
        nodelist= self.getMetadataNestedGraph()
        
        print (nodelist)
        src = Source(nodelist)
        print (dir(src))
        src.format='svg'
        src.render('img/plotPipelineStructure')
        
#        for s in nodelist:
#            g1.node(s)
#        g1.edges(zip(nodelist, nodelist[1:]))
#        g1.render('img/plotPipelineStructure')       
#        print(nodelist)
#        print(g1.source)
#        return display(SVG('img/plotPipelineStructure.svg'))

        return display(SVG('img/plotPipelineStructure.svg'))

    #++++++++++ NEW Victor +++++++++++
    def getMetadataNestedTypes(self):
        '''
        Get the Metadata Nested Types

        @return String: Metadata Nested types
        '''
        returnString=""

        for e in self.stage_containers:
            returnString = returnString + e.getMetadataNestedTypes()
            
        return returnString
    
    #++++++++++ NEW Victor +++++++++++
    def getMetadataNestedGraph(self):
        '''
        Retrieve the metadata nested graph

        @return String: Metadata nested graph
        '''
        returnString="digraph { \n graph [rankdir=LR] \n node [fillcolor=\"#d8e9fd\" fontname=Arial fontsize=9 height=0.1 margin=0.1 shape=rounded style=filled]\n"

        for e in self.stage_containers:
            returnString = returnString + e.getMetadataNestedGraph()
        
        returnString = returnString + "}"
        
        return returnString


    def _startCluster(self, secret):
        ''' 
        Start dispy cluster

        @param secret: Password for dispy nodes
        '''

        import dispy
        
        # Amazon run function
        def amazon_run(data_fetcher, stage_containers, run_id=-1):
            data_container = data_fetcher.output()
            data_container.run_id = run_id
            for s in stage_containers:
                s.run(data_container)

            return data_container.getResults()

        self.__cluster = dispy.SharedJobCluster(amazon_run, secret=secret, port=0,
                                                scheduler_node='127.0.0.1',
                                                ip_addr='127.0.0.1')

        # for i in range(1000):
        #     try:
        #         self._dispy_http = DispyHTTPServer(self.__cluster,port = 8181 + i)
        #         self._dispy_http.port = 8181 + i
        #         break
        #     except OSError:
        #         pass

            
    def _createDispyLink(self):
        ''' Create a link to the Dispy Cluster Information website '''
        
        print("Access cluster monitor at http://hostname:8181 where hostname is the address of the primary host")
        
    def __str__(self):
        '''
        String representation of the pipeline.

        @return String of current metadata of pipeline containers.
        '''
        return str(self.getMetadata())
