
import numpy as np


################################
######## GENERAL PARAMS ######## 
################################

srate = 500

# initial sujet_list
# sujet_list = ['NS086_02', 'NS131_02', 'NS196', 'NS202', 'NS203', 'NS205', 'NS206', 'NS208',
#               'NS211',  'NS213',  'NS215',
#               ]

#NS196, trials num in data do not correspond to trials num in df_resp_cycle, removed
#NS086_02, no flow removed

#sujet list corrected
sujet_list = ['NS131_02', 'NS202', 'NS203', 'NS205', 'NS206', 'NS208',
              'NS211',  'NS213',  'NS215']

sujet_list_interaction_analysis = ['NS131_02', 'NS206', 'NS211']


ROI_short_list = ['Amygdala', 'insula-ant', 'lateralorbitofrontal', 'Hippocampus', 'insula-pos', 'medialorbitofrontal', 'insula', 'postcentral',
                      'precentral', 'Brain-Stem']
sujet_thresh = 3


cond_list = ['VS', 'CHARGE']
cond_list_interaction = ['ctrl', 'CO2', 'MECA', 'BOTH_MC']


sujet_list_allcond = {'rsp_ctrl' : ['NS086_02', 'NS131_02', 'NS202', 'NS203', 'NS205', 'NS206', 'NS208', 'NS211', 'NS213', 'NS215'], 
                      'rsp_chl' : ['NS131_02', 'NS202', 'NS203', 'NS205', 'NS206', 'NS208', 'NS211', 'NS213', 'NS215'], 
                      'oc_ctrl' : ['NS086_02', 'NS131_02', 'NS202', 'NS203', 'NS205', 'NS206', 'NS208', 'NS211', 'NS213', 'NS215'], 
                      'oc_chl' : ['NS131_02', 'NS202', 'NS203', 'NS205', 'NS206', 'NS208', 'NS211', 'NS213', 'NS215']}

conditions = ['rsp_ctrl', 'rsp_chl', 'oc_ctrl', 'oc_chl']
phase_cycle_list = ['whole', 'inspi', 'expi']

freq_band_dict = {'theta' : [4,8], 'alpha' : [8,12], 'beta' : [12,50], 'gamma' : [60, 150]}



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
    n_core = 18

    
path_data = os.path.join(path_general, 'Data')
path_prep = os.path.join(path_general, 'Analyses', 'preprocessing')
path_precompute = os.path.join(path_general, 'Analyses', 'precompute') 
path_results = os.path.join(path_general, 'Analyses', 'results') 
path_memmap = os.path.join(path_general, 'memmap') 

path_paper_figure_export = os.path.join(path_results, 'figure_paper')

os.chdir(init_workdir)



########################################
######## ECG PARAMS ########
########################################

ecg_sujet_pin_choice = {'NS131_02' : [0], 'NS202' : [3], 'NS203' : [2], 'NS205' : ['HR'], 
                        'NS206' : ['HR'], 'NS208' : ['HR'], 'NS211' : [0],  'NS213' : ['HR'],  
                        'NS215' : ['HR']}

sujet_list_ecg = ['NS131_02', 'NS203'] #NS202 removed no baseline, NS211 too bad ECG

ecg_extract_time = 30 #seconds
ecg_lower_limit_nHR = 10

ecg_chunk_time = 5 #seconds

########################################
######## CONSTRUCT LINEAR SIG ########
########################################

respi_cycle_search_time = [-2, 6] #sec





########################################
######## RF analysis ########
########################################

rf_metric_list = ['duration', 'cycle_freq', 'amplitude', 'volume', 'oc_ratio', 'oc_val']





################################
######## PRECOMPUTE TF ########
################################

#### chunk data
chunk_time = 2 #sec
exclude_cycle_duration_jules = 2 #sec
exclude_inspi_duration_jules = 1 #sec

#### stretch
stretch_point_TF = 250
stretch_TF_auto = False
ratio_stretch_TF = 0.50

#### TF & ITPC
nfrex = 150
ncycle_list = [7, 41]
freq_list = [2, 150]
srate_dw = 10
wavetime = np.arange(-3,3,1/srate)
frex = np.logspace(np.log10(freq_list[0]), np.log10(freq_list[1]), nfrex) 
cycles = np.logspace(np.log10(ncycle_list[0]), np.log10(ncycle_list[1]), nfrex).astype('int')
Pxx_wavelet_norm = 1000

#### STATS
n_surrogates_tf = 1000
tf_percentile_sel_stats = 2 # for both side
tf_stats_percentile_cluster = 95
norm_method = 'rscore'# 'zscore', 'dB'
tf_stats_percentile_cluster_allplot = 99
df_extraction_Cxy = 0.02 #Hz around respi median


#### plot
tf_plot_percentile_scale = 1 #for one side





########################################
######## POWER EXTRACTION ########
########################################

pre_extraction_time = -4 #second
peak_time_jitter = 0.05 #second


