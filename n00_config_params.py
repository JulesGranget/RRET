
import numpy as np


################################
######## GENERAL PARAMS ######## 
################################

srate = 500


sujet_list = ['NS086_02', 'NS131_02', 'NS196', 'NS202', 'NS203', 'NS205', 'NS206', 'NS208',
              'NS211',  'NS213',  'NS215',
              ]


cond_list = ['VS', 'CHARGE']





########################################
######## PATH DEFINITION ########
########################################

import socket
import os
import platform
 
PC_OS = platform.system()
PC_ID = socket.gethostname()
init_workdir = os.getcwd()

if PC_ID == 'jules-ubuntu1':

    path_main_workdir = '/home/jules/Documents/RRET_JULES/Scripts'
    path_general = '/home/jules/Documents/RRET_JULES'
    n_core = 25

    
path_data = os.path.join(path_general, 'Data')
path_prep = os.path.join(path_general, 'Analyses', 'preprocessing')
path_precompute = os.path.join(path_general, 'Analyses', 'precompute') 
path_results = os.path.join(path_general, 'Analyses', 'results') 

os.chdir(init_workdir)











