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

import random

class StageContainer:
    ''' Container to hold a stage for the DiscoveryPipeline. '''
    
    def __init__(self, obj_content, obj_runmethod = None, obj_perturbmethod = None, obj_reset = None):
        '''
        Get the object and its run method into this conainer.
        
        @param obj_content: filter, analysis, or accumlator
        @param obj_runmethod: Run method of the obj_content (default process)
        @param obj_perturbmethod: Perturb method of the obj_content (default peturbParams)
        @param obj_reset: Reset method of the obj_content (default resetParams)
        '''
        self.obj_content    = obj_content

        if obj_runmethod != None:
            self.runmethod      = obj_runmethod
        else:
            self.runmethod = self.obj_content.process

        if obj_perturbmethod != None:
            self.perturbmethod  = obj_perturbmethod
        else:
            self.perturbmethod  = self.obj_content.perturbParams


        if obj_reset != None:
            self.resetmethod  = obj_reset
        else:
            self.resetmethod  = self.obj_content.resetParams

        
    def run(self, obj_data_container):
        '''
        Execute the obj_content run method

        @param obj_data_container: Data container to be passed to the held obj_content's run method
        '''
        self.runmethod(obj_data_container)

    def perturb(self):
        ''' Execute the obj_content peturb method '''
        self.perturbmethod()

    def reset(self):
        ''' Execute the obj_content reset method '''
        self.resetmethod()

    def getMetadata(self):
        '''
        Retrieves the obj_content metadata
        
        @return obj_content metadata
        '''
        return self.obj_content.getMetadata()

    def getObjects(self):
        '''
        Return the obj_content in a list 

        @return Contained object in a list
        '''
        return [self.obj_content]

    #+++++++ NEW Victor ++++++++
    def getMetadataType(self):
        '''
        Get metadata type

        @return String: container type
        '''
        return "SC"


    #+++++++ NEW Victor ++++++++
    def getMetadataNestedTypes(self):
        '''
        Get the metadata along with container type

        @return string of container and metadata
        '''
        returnString= "#SC("+ self.getMetadata() +")"        
        return returnString 

    #+++++++ NEW Victor ++++++++
    def getMetadataNestedGraph(self):
        '''
        Get the nested graph for the container

        @return String: Stage container subgraph
        '''
        returnString= "\n subgraph cluster { \n \""+ self.getMetadata() +"\"; \n label=\"SC\" }\n"        
        return returnString 

    
class StageContainerAlternative:
    ''' 
    Stage Container that holds a list of stage containers and randomly chooses one to use.
    '''
    
    currentContainer =[]
    
    def __init__(self, list_stagecontainers):
        '''
        Initialize the StageContainerAlternative

        @param list_stagecontainers: List of stage containers
        '''
        
        self.list_stagecontainers = list_stagecontainers
        self.currentContainer = self.list_stagecontainers[0]
#        self.currentContainer = random.choice(self.list_stagecontainers)
#        self.perturb()
        
    def run(self, obj_data_container):
        ''' 
        Run the currently selected stage container.

        @param obj_data_container: Data container to be passed to the current stagecontainer
        '''
        self.currentContainer.run(obj_data_container)

    def perturb(self):
        ''' choose one of the containers as an alternative and perturb its parameters'''
        self.currentContainer = random.choice(self.list_stagecontainers)
        self.currentContainer.perturb()

    def getMetadata(self):
        ''' 
        Return metadata from the current container 

        @return metadata from the currently selected container
        '''
        return self.currentContainer.getMetadata()

    def getObjects(self):
        '''
        retrieve the current container as a list

        @return Current container being used as a list
        '''

        return self.currentContainer.getObjects()

    def reset(self):
        '''
        Reset the current chosen StageContainer

        self.currentContainer.reset()
        '''

    #++++++++ NEW Victor +++++++++++++++++
    def getMetadataType(self):
        '''
        Get metadata type

        @return String: container type
        '''
        return "SC_Alternative"

    #++++++++ NEW Victor +++++++++++++++++
    def getMetadataNestedTypes(self):
        '''
        Get the metadata along with container type

        @return string of container and metadata
        '''
        
        returnString="#ALT("
        for e in self.list_stagecontainers:
                returnString = returnString + e.getMetadataNestedTypes()
        returnString=returnString + ")"
        
        return returnString 

    #++++++++ NEW Victor +++++++++++++++++
    def getMetadataNestedGraph(self):
        '''
        Get the nested graph for the container

        @return String: Container subgraph
        '''
        
        returnString="\n subgraph cluster { \n"
        for e in self.list_stagecontainers:
                returnString = returnString + "\"" + e.getMetadataNestedTypes() +"\"; \n"
        returnString=returnString + " label=\"ALT\" } \n"
        
        return returnString 
        
        


class StageContainerIncrementalAdd:
    ''' In each perturb call, it incrementally adds one of the filters specified in the constructor. '''
    

    def __init__(self, list_stagecontainers):
        '''
        Initialize the container

        @param list_stagecontainers: List of stage containers.
        '''
        self.length = len(list_stagecontainers)
        self.list_AllStagecontainers = list_stagecontainers

        self.currentindex = 0
        self.list_currentContainers      = []

        self.perturb()
        
    def reset(self):
        ''' Reset the container so that it will only run the first stage container again. '''
        self.list_currentContainers     = []
        self.currentindex =0
        self.perturb()
    
    def run(self, obj_data_container):
        ''' Run the current list of stage containers'''
        for c in self.list_currentContainers:
            c.run(obj_data_container)

    def perturb(self):
        ''' Add another stage container to the current list of stage containers '''
        if self.currentindex < self.length:
            self.list_currentContainers.append(self.list_AllStagecontainers[self.currentindex])
            self.currentindex+=1
        
        # for c in self.list_currentContainers:
        #     c.perturb()
            
    def getMetadata(self):
        '''
        Return the metadata from the currently used stage containers.
        
        @return List of metadata from current containers
        '''
        metadataList = [s.getMetadata() for s in self.list_currentContainers]
        return metadataList

    def getObjects(self):
        '''
        Retrieve objects in the current list of stage containers

        @return List of current obj_content from the current list of stage containers
        '''
        obj_list = []
        for container in self.list_currentContainers:
            for obj in container.getObjects():
                obj_list.append(obj)

        return obj_list

    #++++++++ NEW Victor +++++++++++++++++
    def getMetadataType(self):
        '''
        Get metadata type

        @return String: container type
        '''        
        
        return "SC_IncrementalAdd"
    
    #++++++++ NEW Victor +++++++++++++++++
    def getMetadataNestedTypes(self):
        '''
        Get the metadata along with container type

        @return string of container and metadata
        '''
        returnString="#INCA("
        for e in self.list_stagecontainers:
                returnString = returnString + e.getMetadataNestedTypes()
        returnString=returnString + ")"

    #++++++++ NEW Victor +++++++++++++++++
    def getMetadataNestedGraph(self):
        '''
        Get the nested graph for the container

        @return String: Container subgraph
        '''
        
        returnString="\n subgraph cluster { \n"
        for e in self.list_stagecontainers:
                returnString = returnString + "\""+e.getMetadataNestedTypes()+"\"; \n"
        returnString=returnString + "\n label=\"INCA\" } \n"
        
        return returnString 
