



from n00_config_params import *
from n00bis_config_analysis_functions import *
from n01_manip_data import *







################################
######## PRECOMPUTE ########
################################

#sujet, cond = 'NS086_02', 'rsp_ctrl'
def precompute_respfeatures(sujet, cond):

    #### verify if already computed
    os.chdir(os.path.join(path_precompute, 'RESP'))

    if os.path.exists(f'{sujet}_{cond}_df_respfeatures_pre.xlsx'):
        print(f'{sujet} {cond} ALREADY COMPUTED', flush=True)
        return

    print(f'RESP PRECOMPUTE {sujet} {cond}', flush=True)

    #### get params
    data_allcond, resp_allcond, chanlist, localist = get_data_sujet(sujet)
    data, resp = data_allcond[cond], resp_allcond[cond]
    del data_allcond, resp_allcond

    time_vec_sec = np.arange(resp[0].size)/srate - (resp[0].size/srate)/2

    #### inspect
    if debug:

        _rast = data[0]
        for cycle_i in range(_rast.shape[0]):
            plt.plot(_rast[cycle_i])
        plt.show()

    _half_time = np.where(np.isclose(time_vec_sec, 0, atol=1e-3))[0][0]

        #### generate cycles stretch
    cycles_pre = []
    cycles_post = []

    for epoch_i in range(resp.shape[0]):

        _resp = np.array(resp[epoch_i], copy=True)
        _resp_clean = physio.preprocess(_resp, srate, band=25., btype='lowpass', ftype='bessel', order=5, normalize=False)
        _resp_clean_smooth = physio.smooth_signal(_resp_clean, srate, win_shape='gaussian', sigma_ms=40.0)
        
        try:
            _cycles = physio.detect_respiration_cycles(_resp_clean_smooth, srate, baseline_mode='median')
        except:
            continue

        if debug:

            fig_verif, ax = plt.subplots(figsize=(18, 10))
            ax.plot(_resp)
            ax.scatter(_cycles[:,0], _resp[_cycles[:,0]], color='g', label='inspi_selected')
            ax.scatter(_cycles[:,-1], _resp[_cycles[:,-1]], color='g', label='inspi_selected')
            ax.scatter(_cycles[:,1], _resp[_cycles[:,1]], color='c', label='expi_selected', marker='s')
            plt.legend()
            fig_verif.show()   

        if (_cycles[:,1] <= _half_time).sum():

            _cycles_pre_sel = _cycles[:,1] <= _half_time
            _cycles_pre = _cycles[_cycles_pre_sel]
            _cycles_pre = _cycles_pre[-1]
            _cycles_pre[-1] = np.where(np.isclose(time_vec_sec, 0, atol=1e-3))[0][0]

            cycles_pre.append([_cycles_pre, epoch_i])

        if (_cycles[:,1] >= _half_time).sum():
            
            _cycles_post_sel = _cycles[:,1] >= _half_time
            _cycles_post = _cycles[_cycles_post_sel]    
            _cycles_post[0,0] = np.where(np.isclose(time_vec_sec, 0, atol=1e-3))[0][0]
            _cycles_post = _cycles_post[0]

            cycles_post.append([_cycles_post, epoch_i])

    os.chdir(os.path.join(path_precompute, 'RESP')) 
    epoch_sel_pre = np.array([_cycle_i for _, _cycle_i in cycles_pre])
    epoch_sel_post = np.array([_cycle_i for _, _cycle_i in cycles_post])
    np.save(f'{sujet}_{cond}_epoch_sel_pre.xlsx', epoch_sel_pre) 
    np.save(f'{sujet}_{cond}_epoch_sel_post.xlsx', epoch_sel_post) 

        #### pre
    cycle_stretch_pre = []

    for cycle_i, _cycle in enumerate(cycles_pre):

        if cycle_i == 0:
            first_point = _cycle[0][0]
            cycle_stretch_pre.append(_cycle[0] - first_point)
        else:
            _previous_end = cycle_stretch_pre[cycle_i-1][-1]
            first_point = _cycle[0][0]
            cycle_stretch_pre.append(_cycle[0] - first_point + _previous_end)

    cycle_stretch_pre = np.array(cycle_stretch_pre)

        #### post
    cycle_stretch_post = []

    for cycle_i, _cycle in enumerate(cycles_post):

        if cycle_i == 0:
            cycle_stretch_post.append(_cycle[0] - _half_time)
        else:
            _previous_end = cycle_stretch_post[cycle_i-1][-1]
            cycle_stretch_post.append(_cycle[0] - _half_time + _previous_end)

    cycle_stretch_post = np.array(cycle_stretch_post)

        #### extract values
    _resp_pre = np.zeros((cycle_stretch_pre[-1,-1]))
    _resp_post = np.zeros((cycle_stretch_post[-1,-1]))
        
        ### pre
    print('pre')
    for cycle_i in range(cycle_stretch_pre.shape[0]):

        print_advancement(cycle_i, cycle_stretch_pre.shape[0], steps=[25,50,75])

        _start, _stop, _epoch = cycles_pre[cycle_i][0][0], cycles_pre[cycle_i][0][-1], cycles_pre[cycle_i][-1]
        _resp_pre[cycle_stretch_pre[cycle_i][0]:cycle_stretch_pre[cycle_i][-1]] = resp[_epoch][..., _start:_stop]

        #### post
    print('post')
    for cycle_i in range(cycle_stretch_post.shape[0]):

        print_advancement(cycle_i, cycle_stretch_pre.shape[0], steps=[25,50,75])

        _start, _stop, _epoch = cycles_post[cycle_i][0][0], cycles_post[cycle_i][0][-1], cycles_post[cycle_i][-1]
        _resp_post[cycle_stretch_post[cycle_i][0]:cycle_stretch_post[cycle_i][-1]] = resp[_epoch][..., _start:_stop]
            
    if debug:
        plt.plot(_resp_pre)
        plt.show()

        plt.plot(_resp_post)
        plt.show()
        
    #### compute respfeatures
    cycle_stretch_pre[-1, -1] = cycle_stretch_pre[-1, -1]-1
    cycle_stretch_post[-1, -1] = cycle_stretch_post[-1, -1]-1

    resp_features_pre = physio.compute_respiration_cycle_features(_resp_pre, srate, cycle_stretch_pre, baseline=None)    
    resp_features_post = physio.compute_respiration_cycle_features(_resp_post, srate, cycle_stretch_post, baseline=None)    

    #### save stretch tf
    print(f'SAVE', flush=True)
    os.chdir(os.path.join(path_precompute, 'RESP')) 
    resp_features_pre.to_excel(f'{sujet}_{cond}_df_respfeatures_pre.xlsx')    
    resp_features_post.to_excel(f'{sujet}_{cond}_df_respfeatures_post.xlsx')    

    






################################
######## EXECUTE ########
################################


if __name__ == '__main__':


    for cond in conditions:
    
        for sujet in sujet_list_allcond[cond]:

            precompute_respfeatures(sujet, cond)



                        