



from A_config.n01_O_params import *
from A_config.n02_O_analysis_functions import *
from A_config.n03_X_manip_data import *









############################################
######## OC EXTRACTION FUNCTION ########
############################################



def get_oc_values(sig):

    sig_diff = np.diff(sig)
    start_inspi_dec = np.argmin(sig_diff[:int(stretch_point_TF/4)])
    start_inspi_asc = np.argmax(sig_diff[int(stretch_point_TF/4):int(stretch_point_TF/2)]) + int(stretch_point_TF/4)

    oc_top_i = np.where(sig_diff[start_inspi_dec:] > 0)[0][0] + start_inspi_dec
    oc_top_val = sig[oc_top_i]
    median_inspi_trough = np.median(sig[start_inspi_dec:start_inspi_asc])

    _oc_ratio = (oc_top_val - median_inspi_trough) / median_inspi_trough

    return _oc_ratio, oc_top_val








####################################
######## EXTRACT OC SIZE ########
####################################

#sujet = sujet_list[0]
def extract_oc_size_general(sujet):

    print(sujet)

    #### load
    os.chdir(os.path.join(path_precompute, 'RESP', 'respfeatures')) 
    respfeatures = pd.read_excel(f'{sujet}_respfeatures_cleaned_label.xlsx')

    #### extract
    df_oc = []

    #cond = conditions[-1]
    for cond in conditions:

        os.chdir(os.path.join(path_precompute, 'RESP', 'session'))
        resp_stretch = np.load(f"{sujet}_{cond}_stretch_resp_post.npy")
        _respfeature = respfeatures.query(f"cond == '{cond}'")
        
        oc_ratio = []
        oc_val = []

        for cycle_i in range(resp_stretch.shape[0]): 

            sig = resp_stretch[cycle_i]

            _oc_ratio, oc_top_val = get_oc_values(sig)
            
            oc_ratio.append(_oc_ratio)
            oc_val.append(oc_top_val)

            # if _oc_ratio > 10 and cond == 'rsp_chl':
            #     raise

            if debug:

                time_vec = np.arange(sig.size)[:-1]
                plt.plot(time_vec, sig[:-1])
                plt.vlines([start_inspi_asc, start_inspi_dec], ymin=sig.min(), ymax=sig.max(), colors='g')
                plt.hlines([median_inspi_trough], xmin=0, xmax=stretch_point_TF, colors='b')
                plt.scatter([oc_top_i], sig[oc_top_i], color='r')
                plt.plot(time_vec,sig_diff)
                plt.show()

        if debug:

            for cycle_i in range(resp_stretch.shape[0]):

                plt.plot(resp_stretch[cycle_i])
            
            plt.show()

            for cycle_i in range(resp_stretch.shape[0]):

                plt.plot(np.diff(resp_stretch[cycle_i]))
            
            plt.show()

            cycle_i = 0

            for cycle_i in range(resp_stretch.shape[0]): 

                sig = resp_stretch[cycle_i][:-1]
                sig_diff = np.diff(resp_stretch[cycle_i])
                start_inspi_dec = np.argmin(sig_diff[:int(stretch_point_TF/4)])
                start_inspi_asc = np.argmax(sig_diff[int(stretch_point_TF/4):int(stretch_point_TF/2)]) + int(stretch_point_TF/4)

                oc_top_i = np.where(sig_diff[start_inspi_dec:] > 0)[0][0] + start_inspi_dec

                median_inspi_trough = np.median(sig[start_inspi_dec:start_inspi_asc])

                time_vec = np.arange(sig.size)
                plt.plot(time_vec, sig)
                plt.vlines([start_inspi_asc, start_inspi_dec], ymin=sig.min(), ymax=sig.max(), colors='g')
                plt.hlines([median_inspi_trough], xmin=0, xmax=stretch_point_TF, colors='b')
                plt.scatter([oc_top_i], sig[oc_top_i], color='r')
                plt.plot(time_vec,sig_diff)
                plt.show()

            for cycle_i in range(resp_stretch.shape[0]): 
            
                sig = resp_stretch[cycle_i][:-1]
                sig_diff = np.diff(resp_stretch[cycle_i])
                start_inspi_dec = np.argmin(sig_diff[:int(stretch_point_TF/4)])
                start_inspi_asc = np.argmax(sig_diff[int(stretch_point_TF/4):int(stretch_point_TF/2)]) + int(stretch_point_TF/4)

                oc_top_i = np.where(sig_diff[start_inspi_dec:] > 0)[0][0] + start_inspi_dec

                median_inspi_trough = np.median(sig[start_inspi_dec:start_inspi_asc])

                time_vec = np.arange(sig.size)
                plt.plot(time_vec, sig)
                plt.vlines([int(stretch_point_TF/2)], ymin=sig.min(), ymax=sig.max(), colors='r')
                plt.hlines([median_inspi_trough], xmin=0, xmax=stretch_point_TF, colors='b')
                plt.scatter([oc_top_i], sig[oc_top_i], color='r')
                plt.show()

        _df_oc = pd.DataFrame({'sujet' : [sujet]*len(oc_ratio), 'cond' : [cond]*len(oc_ratio), 'cycle_i' : np.arange(len(oc_ratio)).tolist(), 'oc_ratio' : oc_ratio, 'oc_val' : oc_val})

        if _respfeature.shape[0] != _df_oc.shape[0]:

            raise ValueError('!!! NOT SAME CYCLE NUMBER !!!')

        df_oc.append(_df_oc)

    df_oc = pd.concat(df_oc)

    #### plot
    df_plot = pd.melt(df_oc, id_vars=[_col for _col in df_oc.columns if _col not in ['oc_ratio', 'oc_val']], value_vars=['oc_ratio', 'oc_val'], var_name="metric_type", value_name="val",)
    g = sns.catplot(df_plot, kind='strip', x='cond', y='val', col='metric_type', jitter=True, sharey=False)
    plt.suptitle(f"{sujet}")
    # plt.show()


    #### save
    os.chdir(os.path.join(path_results, 'respi', 'oc_ratio'))
    g.savefig(f"{sujet}_oc.png")
    plt.close('all')

    os.chdir(os.path.join(path_precompute, 'RESP', 'respfeatures')) 
    df_oc.to_excel(f"{sujet}_df_oc.xlsx")

    

#sujet = sujet_list[0]
def extract_oc_size_MECACO2(sujet):

    print(sujet)

    #### extract
    df_oc = []

    #cond = cond_list_interaction[2]
    for cond in cond_list_interaction:

        os.chdir(os.path.join(path_precompute, 'RESP', 'MECACO2INTER'))
        resp_stretch = np.load(f"{sujet}_oc_{cond}_stretch_resp_MECACO2INTER.npy")    
        
        oc_ratio = []
        oc_val = []

        for cycle_i in range(resp_stretch.shape[0]): 

            sig = resp_stretch[cycle_i]

            _oc_ratio, oc_top_val = get_oc_values(sig)
            
            oc_ratio.append(_oc_ratio)
            oc_val.append(oc_top_val)

            # if _oc_ratio > 10 and cond == 'rsp_chl':
            #     raise

            if debug:

                time_vec = np.arange(sig.size)[:-1]
                plt.plot(time_vec, sig[:-1])
                plt.vlines([start_inspi_asc, start_inspi_dec], ymin=sig.min(), ymax=sig.max(), colors='g')
                plt.hlines([median_inspi_trough], xmin=0, xmax=stretch_point_TF, colors='b')
                plt.scatter([oc_top_i], sig[oc_top_i], color='r')
                plt.plot(time_vec,sig_diff)
                plt.show()

        if debug:

            for cycle_i in range(resp_stretch.shape[0]):

                plt.plot(resp_stretch[cycle_i])
            
            plt.show()

            for cycle_i in range(resp_stretch.shape[0]):

                plt.plot(np.diff(resp_stretch[cycle_i]))
            
            plt.show()

            cycle_i = 0

            for cycle_i in range(resp_stretch.shape[0]): 

                sig = resp_stretch[cycle_i][:-1]
                sig_diff = np.diff(resp_stretch[cycle_i])
                start_inspi_dec = np.argmin(sig_diff[:int(stretch_point_TF/4)])
                start_inspi_asc = np.argmax(sig_diff[int(stretch_point_TF/4):int(stretch_point_TF/2)]) + int(stretch_point_TF/4)

                oc_top_i = np.where(sig_diff[start_inspi_dec:] > 0)[0][0] + start_inspi_dec

                median_inspi_trough = np.median(sig[start_inspi_dec:start_inspi_asc])

                time_vec = np.arange(sig.size)
                plt.plot(time_vec, sig)
                plt.vlines([start_inspi_asc, start_inspi_dec], ymin=sig.min(), ymax=sig.max(), colors='g')
                plt.hlines([median_inspi_trough], xmin=0, xmax=stretch_point_TF, colors='b')
                plt.scatter([oc_top_i], sig[oc_top_i], color='r')
                plt.plot(time_vec,sig_diff)
                plt.show()

        _df_oc = pd.DataFrame({'sujet' : [sujet]*len(oc_ratio), 'cond' : [cond]*len(oc_ratio), 'cycle_i' : np.arange(len(oc_ratio)).tolist(), 'oc_ratio' : oc_ratio, 'oc_val' : oc_val})

        df_oc.append(_df_oc)

    df_oc = pd.concat(df_oc)

    #### plot
    df_plot = pd.melt(df_oc, id_vars=[_col for _col in df_oc.columns if _col not in ['oc_ratio', 'oc_val']], value_vars=['oc_ratio', 'oc_val'], var_name="metric_type", value_name="val",)
    g = sns.catplot(df_plot, kind='strip', x='cond', y='val', col='metric_type', jitter=True, sharey=False)
    plt.suptitle(f"{sujet}")
    # plt.show()


    #### save
    os.chdir(os.path.join(path_results, 'respi', 'oc_ratio'))
    g.savefig(f"{sujet}_MECACO2_oc.png")
    plt.close('all')

    os.chdir(os.path.join(path_precompute, 'RESP', 'respfeatures')) 
    df_oc.to_excel(f"{sujet}_MECACO2_df_oc.xlsx")

    


def generate_df_R():

    #### load resp features
    rf_metrics_raw = ['cycle_duration', 'inspi_duration', 'expi_duration', 'cycle_freq', 'inspi_volume',
       'expi_volume', 'total_amplitude', 'inspi_amplitude', 'expi_amplitude',
       'total_volume']
    phase_cycle_list = ['whole', 'inspi', 'expi']
    mapping_pase_cycle_list = {'whole' : ['cycle_duration', 'cycle_freq', 'total_amplitude', 'total_volume', 'oc_ratio', 'oc_val'], 
                               'inspi' : ['inspi_duration', 'inspi_cycle_freq', 'inspi_amplitude', 'inspi_volume', 'oc_ratio', 'oc_val'], 
                               'expi' : ['expi_duration', 'expi_cycle_freq', 'expi_amplitude', 'expi_volume', 'oc_ratio', 'oc_val']}
    rf_metric_allphase_cycle = ['duration', 'cycle_freq', 'amplitude', 'volume', 'oc_ratio', 'oc_val']

    R_keep_OC_metric = ['inspi_volume', 'expi_volume', 'total_amplitude', 'inspi_amplitude', 'expi_amplitude', 'total_volume', 'oc_ratio', 'oc_val']
    
    df_rf_allcond = []
    df_rf_allcond_OC = []

    #sujet = sujet_list[0]
    for sujet in sujet_list:

        # load respiration features
        os.chdir(os.path.join(path_precompute, 'RESP', 'respfeatures')) 
        _oc_respfeatures = pd.read_excel(f'{sujet}_df_oc.xlsx')

        _respfeatures = pd.read_excel(f'{sujet}_respfeatures_cleaned_label.xlsx')
        _respfeatures_pre = pd.read_excel(f'{sujet}_respfeatures_cleaned_label_pre.xlsx')
        
        for cond_i, cond in enumerate(conditions):

            resp_pre  = _respfeatures_pre.query(f"cond == '{cond}'").drop(columns=['Unnamed: 0'])[rf_metrics_raw]
            resp_post = _respfeatures.query(f"cond == '{cond}'").drop(columns=['Unnamed: 0'])[rf_metrics_raw]

            resp_pre[['oc_ratio', 'oc_val']] = _oc_respfeatures.query(f"cond == '{cond}'")[['oc_ratio', 'oc_val']]
            resp_post[['oc_ratio', 'oc_val']] = _oc_respfeatures.query(f"cond == '{cond}'")[['oc_ratio', 'oc_val']]

            resp_pre['inspi_cycle_freq'], resp_pre['expi_cycle_freq'] = 1/resp_pre['inspi_duration'], 1/resp_pre['expi_duration']
            resp_post['inspi_cycle_freq'], resp_post['expi_cycle_freq'] = 1/resp_post['inspi_duration'], 1/resp_post['expi_duration']
            
            #### allmetric
            for _phase_cycle in phase_cycle_list:

                _df_map_metric_pre = resp_pre[mapping_pase_cycle_list[_phase_cycle]].copy()
                _df_map_metric_pre['phase'] = [_phase_cycle] * _df_map_metric_pre.shape[0]
                _df_map_metric_pre['pre_post'] = ['pre'] * _df_map_metric_pre.shape[0]
                _df_map_metric_pre['cond'] = [cond] * _df_map_metric_pre.shape[0]
                _df_map_metric_pre['sujet'] = [sujet] * _df_map_metric_pre.shape[0]

                _df_map_metric_post = resp_post[mapping_pase_cycle_list[_phase_cycle]].copy()
                _df_map_metric_post['phase'] = [_phase_cycle] * _df_map_metric_post.shape[0]
                _df_map_metric_post['pre_post'] = ['post'] * _df_map_metric_post.shape[0]
                _df_map_metric_post['cond'] = [cond] * _df_map_metric_post.shape[0]
                _df_map_metric_post['sujet'] = [sujet] * _df_map_metric_post.shape[0]

                if _phase_cycle == 'expi':
                    _rep_dict_col = {col : col[5:] for col in _df_map_metric_pre.columns if col.find(_phase_cycle) != -1}
                if _phase_cycle == 'inspi':
                    _rep_dict_col = {col : col[6:] for col in _df_map_metric_pre.columns if col.find(_phase_cycle) != -1}
                if _phase_cycle == 'whole':
                    _rep_dict_col = {'cycle_duration' : 'duration', 'total_amplitude' : 'amplitude', 'total_volume' : 'volume'}

                _df_map_metric_pre = _df_map_metric_pre.rename(columns=_rep_dict_col)
                _df_map_metric_post = _df_map_metric_post.rename(columns=_rep_dict_col)

                _df_map_prepost = pd.concat([_df_map_metric_pre, _df_map_metric_post], axis=0)

                df_rf_allcond.append(_df_map_prepost)

            #### OC metric
            resp_pre_short = resp_pre[R_keep_OC_metric]
            _rep_dict_col = {col : f"pre_{col}" for col in resp_pre_short.columns}
            _df_map_metric_pre = resp_pre_short.rename(columns=_rep_dict_col)

            resp_post_short = resp_post[R_keep_OC_metric]
            _rep_dict_col = {col : f"post_{col}" for col in resp_post_short.columns}
            _df_map_metric_post = resp_post_short.rename(columns=_rep_dict_col)

            _df_R_OC = pd.concat([_df_map_metric_pre, _df_map_metric_post], axis=1)
            _df_R_OC = _df_R_OC.drop(columns=['pre_oc_ratio', 'pre_oc_val'])

            _df_R_OC['cond'] = [cond] * _df_R_OC.shape[0]
            _df_R_OC['sujet'] = [sujet] * _df_R_OC.shape[0]
            
            df_rf_allcond_OC.append(_df_R_OC)

    df_rf_allcond = pd.concat(df_rf_allcond)
    df_rf_allcond_OC = pd.concat(df_rf_allcond_OC)

    filename = f"df_R_RFonly_allmetric.xlsx"
    filepath = os.path.join(path_precompute, 'RESP', 'df_R', filename)

    df_rf_allcond.to_excel(filepath)

    filename = f"df_R_RFonly_selOC.xlsx"
    filepath = os.path.join(path_precompute, 'RESP', 'df_R', filename)

    df_rf_allcond_OC.to_excel(filepath)








################################
######## EXECUTE ########
################################


if __name__ == '__main__':

    #### main analysis
    for sujet in sujet_list:

        extract_oc_size_general(sujet)

    generate_df_R()

    ####
    for sujet in sujet_list_interaction_analysis:

        extract_oc_size_MECACO2(sujet)

                        