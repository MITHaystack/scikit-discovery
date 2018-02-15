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


import time
from collections import OrderedDict

import ipywidgets as widgets
from IPython.display import display

import skdiscovery.utilities.cloud.amazon_control as ac
from skdiscovery.data_structure.framework import config

widget_dict = OrderedDict()
disable_list = ['initialize_button', 'cache_button', 'restore_button']
initialized_disabled_list = ['new_num_instances_widget', 'execute_instances_button']
key_value_list = ['aws_id_widget', 'aws_secret_widget', 'aws_region_widget','aws_security_widget',
                  'aws_keyname_widget','aws_pem_widget','aws_image_id', 'instance_type_widget']

initialized=False


# Creates a GUI for setting up the Amazon instances when offloading pipelines.
# The GUI can just be called directly, taking no inputs or outputs.
def init():
    '''Initialize GUI for controlling Amazon instances'''
    
    global widget_dict


    # Setup Options
    widget_dict['aws_id_widget'] = widgets.Text(
        value=None,
        placeholder='AWS Access ID',
        description='AWS ID:',
        disabled=False,
        width='100%'
    )

    widget_dict['aws_secret_widget'] = widgets.Text(
        value = '',
        placeholder='AWS Secret Access Key',
        description='Secret:',
        width='100%'
    )
    widget_dict['aws_region_widget'] = widgets.Text(
        value = '',
        placeholder='us-west-2',
        description='Region:',
        width='100%'
    )
    widget_dict['aws_security_widget'] = widgets.Text(
        value = '',
        placeholder='Security Group Name',
        description='Security: ',
        width='100%'
    )

    widget_dict['aws_keyname_widget'] = widgets.Text(
        value = '',
        placeholder='Name of Key Pair',
        description='Key Name:',
        width='100%'
    )

    widget_dict['aws_pem_widget'] = widgets.Text(
        value = '',
        placeholder='/home/username/security_key.pem',
        description='PEM File:',
        width='100%'
    )
    
    widget_dict['aws_image_id'] = widgets.Text(
        value = '',
        placeholder='ami-xxxxxxxx',
        description='Image ID:',
        width='100%'
    )    
    widget_dict['instance_type_widget'] = widgets.Text(
        placeholder = 't2.medium',
        description='Type:',
        disabled=False,
        width='100%'
    )

    widget_dict['initialize_button'] = widgets.Button(
        description='Initialize',
        disabled = False,
        button_style='success',
        )

    widget_dict['cache_button'] = widgets.Button(
        description='Save credentials',
        disabled=False,
        button_style='success',        
    )

    widget_dict['restore_button'] = widgets.Button(
        description='Restore credentials',
        disabled=False,
        button_style='success',        
    )

    def cache_cred(b):
        
        changeButtonState(enabled=False)        

        def cache_value(in_widget, key):
            #if in_widget.value != '':
            config.writeConfigValue('AWS', key, in_widget.value)

        cache_value(widget_dict['aws_id_widget'], 'ID')
        cache_value(widget_dict['aws_secret_widget'], 'Secret')
        cache_value(widget_dict['aws_region_widget'], 'Region')
        cache_value(widget_dict['aws_security_widget'], 'Security')
        cache_value(widget_dict['aws_keyname_widget'], 'Key_Name')
        cache_value(widget_dict['aws_pem_widget'], 'PEM_File')
        cache_value(widget_dict['aws_image_id'], 'image_id')
        cache_value(widget_dict['instance_type_widget'],'instance_type')

        changeButtonState(enabled=True)        

    widget_dict['cache_button'].on_click(cache_cred)

    def restore_cred(b):
        changeButtonState(enabled=False)
        conf = config.getConfig()
        widget_dict['aws_id_widget'].value = conf.get('AWS','ID',fallback='')
        widget_dict['aws_secret_widget'].value = conf.get('AWS','Secret',fallback='')
        widget_dict['aws_region_widget'].value = conf.get('AWS','Region',fallback='')
        widget_dict['aws_security_widget'].value = conf.get('AWS','Security',fallback='')
        widget_dict['aws_keyname_widget'].value = conf.get('AWS','Key_Name',fallback='')
        widget_dict['aws_pem_widget'].value = conf.get('AWS', 'PEM_File', fallback='')
        widget_dict['aws_image_id'].value = conf.get('AWS', 'image_id', fallback='')
        widget_dict['instance_type_widget'].value = conf.get('AWS','instance_type',fallback='')
        changeButtonState(enabled=True)        

    widget_dict['restore_button'].on_click(restore_cred)
    # shared_cluster_button = widgets.Button(icon='ban')


    def initialize_cluster(b):
        '''
        Initialize connection to Amazon Web Services and any offloading instances

        @param b: Button that was clicked button
        '''

        global initialized

        changeButtonState(enabled=False)
        widget_dict['aws_status_widget'].value = 'Initializing'

        if checkValidValues():
            ac.init(widget_dict['aws_id_widget'].value,
                    widget_dict['aws_secret_widget'].value,
                    widget_dict['aws_region_widget'].value,
                    widget_dict['aws_security_widget'].value,
                    widget_dict['aws_keyname_widget'].value,
                    widget_dict['aws_pem_widget'].value)

            widget_dict['aws_status_widget'].value = str(len(ac.amazon_list)) +' Amazon node(s) running'
            widget_dict['execute_instances_button'].disabled = False
            widget_dict['initialize_button'].disabled = False
            widget_dict['new_num_instances_widget'].disabled = False
            widget_dict['initialize_button'].description = 'Re-Initialize'
            widget_dict['new_num_instances_widget'].value = len(ac.amazon_list)

            initialized=True

        changeButtonState(enabled=True)


    widget_dict['initialize_button'].on_click(initialize_cluster)
    

    # Controlling Instances
    # Label converted to text    
    widget_dict['label_num_instances'] = widgets.Label(value='Set Number of Instances: ')
    widget_dict['new_num_instances_widget'] = widgets.IntSlider(
        disabled = True,
        description = '',
        value = 0,
        min=0,
        max=10,
        step = 1
    )

    # Label converted to text
    # widget_dict['run_label'] = widgets.Label(value='Status:')
    widget_dict['aws_status_widget'] = widgets.Text(
        value="Not Initialized",
        description='Status:',
        disabled=True,
    )

    widget_dict['execute_instances_button'] = widgets.Button(
        description='Execute',
        disabled=True,
        width='25%'
    )


    def execute_instances(b):
        '''
        Changes number of running instances to number in slider widget

        @param b: Button that was clicked
        '''
        changeButtonState(enabled=False)
        if checkValidValues():
            widget_dict['aws_status_widget'].value = 'Initializing'
            ac.setNumInstances(widget_dict['new_num_instances_widget'].value,
                               widget_dict['instance_type_widget'].value, widget_dict['aws_image_id'].value)
            widget_dict['aws_status_widget'].value = str(len(ac.amazon_list)) +' Amazon node(s) running'
        changeButtonState(enabled=True)

    widget_dict['execute_instances_button'].on_click(execute_instances)

    restore_cred(None)
    drawGUI()
    
def drawGUI():
    '''
    Draw the GUI on the screen
    '''
    display(widget_dict['aws_id_widget'])
    display(widget_dict['aws_secret_widget'])
    display(widget_dict['aws_region_widget'])
    display(widget_dict['aws_security_widget'])
    display(widget_dict['aws_keyname_widget'])
    display(widget_dict['aws_pem_widget'])
    display(widget_dict['aws_image_id'])    
    display(widget_dict['instance_type_widget'])
    display(widgets.HBox([widget_dict['initialize_button'], widget_dict['cache_button'],widget_dict['restore_button']]))

    # display(widgets.HBox([widget_dict['label_num_instances'], widget_dict['new_num_instances_widget']]))
    # display(widgets.HBox([widget_dict['run_label'], widget_dict['aws_status_widget']]))

    display(widgets.HBox([widget_dict['label_num_instances'], widget_dict['new_num_instances_widget']]))
    display(widget_dict['aws_status_widget'])
    display(widget_dict['execute_instances_button'])
    # Need to disable after displaying in ipywidgets 7
    widget_dict['execute_instances_button'].disabled=True
    widget_dict['execute_instances_button'].button_style=''
    widget_dict['new_num_instances_widget'].disabled=True


def changeButtonState(enabled=True):
    '''
    Enable or disable the buttons and slider in the GUI

    @param enabled: State to change the buttons to.
    '''
    if enabled==True:
        new_disabled = False
        button_style = 'success'
    else:
        new_disabled = True
        button_style = ''
        
    global widget_dict
    global disable_list
    for label in disable_list:
        widget_dict[label].disabled = new_disabled
        widget_dict[label].button_style = button_style

    if initialized:
        for label in initialized_disabled_list:
            widget_dict[label].disabled = new_disabled
            widget_dict[label].button_style = button_style

        
def checkValidValues():
    '''
    Check if Amazon information is valid

    @return True if all AWS text fields have data in them, false otherwise
    '''
    is_valid = True
    for key in key_value_list:
        if widget_dict[key].value == '':
            is_valid = False
            break

    if is_valid == False:
        prev_message = widget_dict['aws_status_widget'].value
        widget_dict['aws_status_widget'].value = 'Please fill in all information'
        time.sleep(5)

    return is_valid
