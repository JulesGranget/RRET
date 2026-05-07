

import joblib


from A_config.n01_O_params import *
from A_config.n02_O_analysis_functions import *
from A_config.n03_X_manip_data import *







################################
######## PRECOMPUTE TF ########
################################

#sujet, cond, norm_param = sujet_list[0], 'rsp_ctrl', 'rscore'
def precompute_tf_allconv(sujet, cond, norm_param):

    #### verify if already computed
    os.chdir(os.path.join(path_precompute, 'TF', 'epoch'))

    if os.path.exists(f'{sujet}_{cond}_df_Pxx.xlsx'):
        print(f'{sujet} {cond} ALREADY COMPUTED', flush=True)
        return

    print(f'TF PRECOMPUTE {sujet} {cond}', flush=True)

    #### get params
    data_allcond, resp_allcond, chanlist, localist = get_data_sujet(sujet)
    data, resp = data_allcond[cond], resp_allcond[cond]
    del data_allcond, resp_allcond

    time_vec_sec = np.arange(resp[0].size)/srate - (resp[0].size/srate)/2
    chan_list, loca_list = get_chanlist(sujet)

    #### inspect
    if debug:

        _rast = data[0]
        min, max = np.percentile(_tf, 1), np.percentile(_tf, 99)
        for cycle_i in range(_rast.shape[0]):
            plt.plot(_rast[cycle_i])
        plt.show()
    
    #### select wavelet parameters
    wavelets = get_wavelets()

    #### compute
    cycle_counts = resp.shape[0]

    print('CONV', flush=True)

    os.chdir(path_memmap)
    tf_allconv = np.memmap(f'memmap_{sujet}_{cond}_tf_conv.npy', mode="w+", shape=(len(chanlist), cycle_counts, nfrex, data.shape[-1]), dtype=np.float32)
            
    #### conv
    #chan_i = 0
    # for chan_i, _ in enumerate(chanlist):
    def extract_chan_conv(chan_i):

        print_advancement(chan_i, data.shape[0], steps=[25, 50, 75])

        for cycle_i in range(cycle_counts):

            for fi in range(nfrex):
                
                tf_allconv[chan_i, cycle_i, fi,:] = abs(scipy.signal.fftconvolve(data[chan_i,cycle_i,:], wavelets[fi,:], 'same'))**2

    joblib.Parallel(n_jobs = n_core, prefer = 'processes')(joblib.delayed(extract_chan_conv)(chan_i) for chan_i, _ in enumerate(chanlist))

    #### extract and norm
    print(f'NORM params', flush=True)

    os.chdir(path_memmap)
    rscore_param = np.memmap(f'memmap_{sujet}_{cond}_rscore_param.npy', mode="w+", shape=(2, len(chanlist), nfrex), dtype=np.float32)

    if norm_param == 'rscore':

        # for chan_i, chan in enumerate(chanlist):
        def extract_chan_rscore_params(chan_i): 

            print_advancement(chan_i, data.shape[0], steps=[25, 50, 75])

            for fi in range(nfrex):
                
                _med = np.median(tf_allconv[chan_i, :, fi,:])
                _mad = scipy.stats.median_abs_deviation(tf_allconv[chan_i, :, fi,:].reshape(-1), axis=0)

                if debug:

                    _cycles = tf_allconv[chan_i, :, fi,:]

                    for cycle_i in range(_cycles.shape[0]):
                        plt.plot(_cycles[cycle_i])
                    plt.show()

                    min, max = np.percentile(_cycles, 1), np.percentile(_cycles, 99)
                    plt.pcolormesh(_cycles, vmin=min, vmax=max)
                    plt.colorbar()
                    plt.show()

                rscore_param[0,chan_i,fi] = _mad
                rscore_param[1,chan_i,fi] = _med

        joblib.Parallel(n_jobs = n_core, prefer = 'processes')(joblib.delayed(extract_chan_rscore_params)(chan_i) for chan_i, _ in enumerate(chanlist))

        print(f'NORM', flush=True)

        os.chdir(path_memmap)
        tf_allconv_norm = np.memmap(f'memmap_{sujet}_{cond}_tf_conv_norm.npy', mode="w+", shape=(len(chanlist), cycle_counts, nfrex, data.shape[-1]), dtype=np.float32)

        # for chan_i, _ in enumerate(chanlist):
        def norm_chan_conv(chan_i):

            print_advancement(chan_i, data.shape[0], steps=[25, 50, 75])

            for cycle_i in range(cycle_counts):

                _mat = tf_allconv[chan_i, cycle_i]
                _mat_norm = np.zeros(_mat.shape)

                if debug:

                    plt.pcolormesh(_mat)
                    plt.colorbar()
                    plt.show()

                    plt.plot(np.median(_mat[:10], axis=0))
                    plt.show()

                for fi in range(nfrex):

                    _mad, _med = rscore_param[0,chan_i,fi], rscore_param[1,chan_i,fi]

                    _x = _mat[fi]
                    _mat_norm[fi] = (_x-_med) * 0.6745 / _mad

                if debug:

                    plt.pcolormesh(_mat_norm)
                    plt.colorbar()
                    plt.show()

                    _mat_norm.min()
                    _mat_norm.max()

                    plt.hist(_mat_norm.reshape(-1))
                    plt.show()
                
                tf_allconv_norm[chan_i,cycle_i] = _mat_norm

        joblib.Parallel(n_jobs = n_core, prefer = 'processes')(joblib.delayed(norm_chan_conv)(chan_i) for chan_i, _ in enumerate(chanlist))

    if norm_param == None:

        print(f'NORM', flush=True)

        tf_allconv_norm[:] = tf_allconv

    if debug:

        for chan_i in range(10):
            _tf = tf_allconv_norm[chan_i, 0]
            min, max = np.percentile(_tf, 1), np.percentile(_tf, 99)
            plt.pcolormesh(_tf, vmin=min, vmax=max)
            plt.colorbar()
            plt.show()

        plt.hist(_tf.reshape(-1))
        plt.show()

        _tf.reshape(-1).max()

    #### construct dataset
    print('CONSTRUCT DATASET')
    os.chdir(path_memmap)
    resp_memmap = np.memmap(f'memmap_{sujet}_{cond}_resp.npy', mode="w+", shape=resp.shape, dtype=np.float32)
    resp_memmap[:] = resp

    # cycle_count_pre = 0
    # cycle_count_post = 0

    _half_time = np.where(np.isclose(time_vec_sec, 0, atol=1e-3))[0][0]
    thresh_pre_post = _half_time - srate*1.5

        #### generate cycles stretch
    cycles_pre = []
    cycles_post = []

    for epoch_i in range(tf_allconv_norm.shape[1]):

        _resp = np.array(resp_memmap[epoch_i], copy=True)
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
    os.chdir(path_memmap)
    _tf_pre = np.memmap(f'memmap_{sujet}_{cond}_tf_pre.npy', mode="w+", shape=(len(chan_list), nfrex, cycle_stretch_pre[-1,-1]), dtype=np.float32)
    _tf_post = np.memmap(f'memmap_{sujet}_{cond}_tf_post.npy', mode="w+", shape=(len(chan_list), nfrex, cycle_stretch_post[-1,-1]), dtype=np.float32)
    _resp_pre = np.zeros((cycle_stretch_pre[-1,-1]))
    _resp_post = np.zeros((cycle_stretch_post[-1,-1]))
        
        ### pre
    print('pre')
    for cycle_i in range(cycle_stretch_pre.shape[0]):

        print_advancement(cycle_i, cycle_stretch_pre.shape[0], steps=[25,50,75])

        _start, _stop, _epoch = cycles_pre[cycle_i][0][0], cycles_pre[cycle_i][0][-1], cycles_pre[cycle_i][-1]
        _resp_pre[cycle_stretch_pre[cycle_i][0]:cycle_stretch_pre[cycle_i][-1]] = resp_memmap[_epoch][..., _start:_stop]
        _tf_pre[:,:,cycle_stretch_pre[cycle_i][0]:cycle_stretch_pre[cycle_i][-1]] = tf_allconv_norm[:, _epoch, :, _start:_stop]

        #### post
    print('post')
    for cycle_i in range(cycle_stretch_post.shape[0]):

        print_advancement(cycle_i, cycle_stretch_pre.shape[0], steps=[25,50,75])

        _start, _stop, _epoch = cycles_post[cycle_i][0][0], cycles_post[cycle_i][0][-1], cycles_post[cycle_i][-1]
        _resp_post[cycle_stretch_post[cycle_i][0]:cycle_stretch_post[cycle_i][-1]] = resp_memmap[_epoch][..., _start:_stop]
        _tf_post[:,:,cycle_stretch_post[cycle_i][0]:cycle_stretch_post[cycle_i][-1]] = tf_allconv_norm[:, _epoch, :, _start:_stop]
            
    if debug:
        plt.plot(_resp_pre)
        plt.show()

        plt.plot(_resp_post)
        plt.show()
        
    #### stretch
    print('STRETCH')
    cycle_stretch_pre[-1, -1] = cycle_stretch_pre[-1, -1]-1
    cycle_stretch_post[-1, -1] = cycle_stretch_post[-1, -1]-1

    resp_features_pre = physio.compute_respiration_cycle_features(_resp_pre, srate, cycle_stretch_pre, baseline=None)    
    resp_features_post = physio.compute_respiration_cycle_features(_resp_post, srate, cycle_stretch_post, baseline=None)    

    tf_allchan_stretch_pre = []
    tf_allchan_stretch_post = []

    for chan_i, _ in enumerate(chan_list):

        tf_allchan_stretch_pre.append(stretch_data_tf(resp_features_pre, stretch_point_TF, _tf_pre[chan_i], srate)[0].astype(np.float32))
        tf_allchan_stretch_post.append(stretch_data_tf(resp_features_post, stretch_point_TF, _tf_post[chan_i], srate)[0].astype(np.float32))

    tf_allchan_stretch_pre = np.array(tf_allchan_stretch_pre)
    tf_allchan_stretch_post = np.array(tf_allchan_stretch_post)

    resp_stretch_pre = stretch_data(resp_features_pre, stretch_point_TF, _resp_pre, srate)[0]
    resp_stretch_post = stretch_data(resp_features_post, stretch_point_TF, _resp_post, srate)[0]

    if debug:
        for cycle_i in range(resp_stretch_pre.shape[0]):
            plt.plot(resp_stretch_pre[cycle_i])
        plt.show()

        for cycle_i in range(resp_stretch_post.shape[0]):
            plt.plot(resp_stretch_post[cycle_i])
        plt.show()

        plt.pcolormesh(np.median(tf_allchan_stretch_pre[0], axis=[0]))
        plt.show()

        plt.pcolormesh(np.median(tf_allchan_stretch_post[0], axis=[0]))
        plt.show()

    #### save stretch tf
    print(f'SAVE TF STRETCH', flush=True)
    os.chdir(os.path.join(path_precompute, 'TF', 'epoch'))
    np.save(f'{sujet}_{cond}_tf_allchan_stretch_pre.npy', tf_allchan_stretch_pre)    
    np.save(f'{sujet}_{cond}_tf_allchan_stretch_post.npy', tf_allchan_stretch_post)   

    os.chdir(os.path.join(path_precompute, 'RESP', 'epochs')) 
    np.save(f'{sujet}_{cond}_pre_stretch_resp.npy', resp_stretch_pre)    
    np.save(f'{sujet}_{cond}_post_stretch_resp.npy', resp_stretch_post)    

    #### extract power
    print(f'EXTRACT POWER', flush=True)

    df_Pxx = pd.DataFrame()
    pre_phase_list = ['pre_inspi', 'pre_expi']
    post_phase_list = ['inspi', 'expi']

    for chan_i in range(tf_allconv_norm.shape[0]):
    # def extract_Pxx_chan(chan_i):

        print_advancement(chan_i, data.shape[0], steps=[25, 50, 75])

        # _df_chan = pd.DataFrame()

        chan_name = chanlist[chan_i]
        chan_loca = localist[chan_i]

        for epoch_i in range(tf_allchan_stretch_pre.shape[1]):

            #### pre
            _tf = tf_allchan_stretch_pre[chan_i, epoch_i]

            for band, freq in freq_band_dict.items():

                frex_sel = (frex >= freq[0]) & (frex <= freq[-1])

                inspi_sel = np.arange(stretch_point_TF) <= stretch_point_TF/2
                expi_sel = np.arange(stretch_point_TF) > stretch_point_TF/2
                Pxx_pre_inspi_sel = np.median(_tf[frex_sel,:][:,inspi_sel])
                Pxx_pre_expi_sel = np.median(_tf[frex_sel,:][:,expi_sel])

                _df =   pd.DataFrame({'sujet' : [sujet]*len(pre_phase_list), 'cond' : [cond]*len(pre_phase_list), 'chan' : [chan_name]*len(pre_phase_list), 'ROI' : [chan_loca]*len(pre_phase_list),
                                    'band' : [band]*len(pre_phase_list), 'phase' : pre_phase_list, 'cycle' : [epoch_i]*len(pre_phase_list), 'Pxx' : [Pxx_pre_inspi_sel, Pxx_pre_expi_sel]})
                
                df_Pxx = pd.concat([df_Pxx, _df])
                # _df_chan = pd.concat([_df_chan, _df])

        for epoch_i in range(tf_allchan_stretch_post.shape[1]):

            #### post
            _tf = tf_allchan_stretch_post[chan_i, epoch_i]

            for band, freq in freq_band_dict.items():

                frex_sel = (frex >= freq[0]) & (frex <= freq[-1])

                inspi_sel = np.arange(stretch_point_TF) <= stretch_point_TF/2
                expi_sel = np.arange(stretch_point_TF) > stretch_point_TF/2
                Pxx_post_inspi_sel = np.median(_tf[frex_sel,:][:,inspi_sel])
                Pxx_post_expi_sel = np.median(_tf[frex_sel,:][:,expi_sel])

                _df =   pd.DataFrame({'sujet' : [sujet]*len(post_phase_list), 'cond' : [cond]*len(post_phase_list), 'chan' : [chan_name]*len(post_phase_list), 'ROI' : [chan_loca]*len(post_phase_list),
                                    'band' : [band]*len(post_phase_list), 'phase' : post_phase_list, 'cycle' : [epoch_i]*len(post_phase_list), 'Pxx' : [Pxx_post_inspi_sel, Pxx_post_expi_sel]})
                
                df_Pxx = pd.concat([df_Pxx, _df])
                # _df_chan = pd.concat([_df_chan, _df])

        # return _df_chan

    # df_allchan = joblib.Parallel(n_jobs = n_core, prefer = 'processes')(joblib.delayed(extract_Pxx_chan)(chan_i) for chan_i, _ in enumerate(chanlist))

    # for chan_i, _ in enumerate(chanlist):

    #     extract_Pxx_chan(chan_i)

    #### SAVE
    print('SAVE EXTRACT PXX', flush=True)
    os.chdir(os.path.join(path_precompute, 'TF', 'epoch'))
    df_Pxx.to_excel(f'{sujet}_{cond}_df_Pxx.xlsx')    

    #### remove
    os.chdir(path_memmap)
    os.remove(f'memmap_{sujet}_{cond}_tf_conv.npy')
    os.remove(f'memmap_{sujet}_{cond}_rscore_param.npy')
    os.remove(f'memmap_{sujet}_{cond}_tf_conv_norm.npy')
    os.remove(f'memmap_{sujet}_{cond}_resp.npy')
    os.remove(f'memmap_{sujet}_{cond}_tf_pre.npy')
    os.remove(f'memmap_{sujet}_{cond}_tf_post.npy')

    del data, resp, tf_allconv, tf_allconv_norm, _tf_pre, _tf_post

    print('done', flush=True)






################################
######## EXECUTE ########
################################


if __name__ == '__main__':

    param_list = []
    norm_param = 'rscore'

    for cond in conditions:
    
        for sujet in sujet_list_allcond[cond]:

            # param_list.append([sujet, cond, norm_param])
            precompute_tf_allconv(sujet, cond, norm_param)

    # joblib.Parallel(n_jobs = n_core, prefer = 'processes')(joblib.delayed(precompute_tf_allconv)(sujet, cond, norm_param) for sujet, cond in param_list[:3])



    #### identify nan

    if debug:

        df_nan = pd.DataFrame()
            
        for sujet in sujet_list_allcond[cond]:

            #### verify if already computed
            os.chdir(os.path.join(path_precompute, 'TF', 'epoch'))

            print(f'{sujet} {cond}', flush=True)

            #### get params
            data_allcond, resp_allcond, chanlist, localist = get_data_sujet(sujet)

            for cond in conditions:

                data, resp = data_allcond[cond], resp_allcond[cond]

                if np.isnan(data).sum() != 0:

                    time_vec_sec = np.arange(resp[0].size)/srate - (resp[0].size/srate)/2
                    nan_vec = np.isnan(data[0,0])
                    nan_vec[0] = False
                    _start, _stop = np.where(np.diff(nan_vec))[0]
                    _start, _stop = time_vec_sec[_start], time_vec_sec[_stop]
                    _df = pd.DataFrame({'sujet' : [sujet], 'cond' : [cond], 'nan' : [np.isnan(data).sum()], 'start' : [_start], 'stop' : [_stop]})
                    df_nan = pd.concat([df_nan, _df])

                else:

                    _df = pd.DataFrame({'sujet' : [sujet], 'cond' : [cond], 'nan' : [np.isnan(data).sum()], 'start' : [0], 'stop' : [0]})
                    df_nan = pd.concat([df_nan, _df])

        # just NS196 having nan



                        