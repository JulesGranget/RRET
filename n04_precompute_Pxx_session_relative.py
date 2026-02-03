

import joblib


from n00_config_params import *
from n00bis_config_analysis_functions import *
from n01_manip_data import *





################################
######## PRECOMPUTE TF ########
################################

#sujet, norm_param = sujet_list[6], 'rscore'
def precompute_tf_allconv_from_epoch(sujet, norm_param):

    #### verify if already computed
    os.chdir(os.path.join(path_precompute, 'TF', 'session'))

    if os.path.exists(f'{sujet}_xr_Pxx.nc'):
        print(f'{sujet} ALREADY COMPUTED', flush=True)
        return

    print(f'TF PRECOMPUTE {sujet}', flush=True)

    #### get params
    data, resp, chanlist, localist = get_data_sujet_fullsig(sujet)

    time_vec_sec = np.arange(resp[0].size)/srate - (resp[0].size/srate)/2
    chan_list, loca_list = get_chanlist(sujet)

    _half_time = np.where(np.isclose(time_vec_sec, 0, atol=1e-3))[0][0]

    respi_cycle_search_time = [-2, 6] #sec

    os.chdir(os.path.join(path_prep))
    df_resp_cycle = pd.read_excel(f"{sujet}_breathInfo.xlsx")

    #### generate linear signals
    print('CONSTRUCT LINEAR SIG')
    resp_linear = []
    data_linear = []
    inspi_starts = []
    _len_resp = 0

    nan_epochs = []

    for epoch_i in range(resp.shape[0]):

        print_advancement(epoch_i, resp.shape[0], [25,50,75])

        if np.isnan(data[:,epoch_i]).sum() != 0 or np.isnan(resp[epoch_i]).sum() != 0:
            nan_epochs.append(epoch_i)
            continue

        if debug:
            plt.plot(resp[epoch_i])
            plt.show()

            plt.pcolormesh(data[:,epoch_i])
            plt.show()

        if epoch_i == 0:

            resp_linear.append(resp[epoch_i])
            data_linear.append(data[:,epoch_i])
            inspi_starts.append(_half_time)
            _len_resp += resp[epoch_i].size

        else:

            _resp_m1 = resp[epoch_i-1]
            _resp_1 = resp[epoch_i]

            _congruence = []
            for i in range(_resp_m1.size):
                _congruence.append((_resp_m1 == np.roll(_resp_1, i)).sum())
            _congruence = np.array(_congruence)

            _i_roll = np.argmax(_congruence)

            if _i_roll == 0 or srate*1 > _congruence[_i_roll]:
                _i_roll_tot = 0
            else:
                _i_roll_tot = _resp_m1.size - _i_roll

            if debug:
                plt.plot(np.arange(_resp_m1.size), _resp_m1)
                plt.plot(np.arange(np.concat([np.zeros(int(_resp_m1.size/2)), _resp_1]).size), np.concat([np.zeros(int(_resp_m1.size/2)), _resp_1]))
                plt.show()

                plt.plot(_congruence)
                plt.show()

                plt.scatter(np.arange(_resp_m1.size), _resp_m1)
                plt.scatter(np.arange(np.concat([np.zeros(_i_roll), _resp_1]).size), np.concat([np.zeros(_i_roll), _resp_1]))
                plt.show()

                plt.plot(_resp_m1)
                plt.plot(np.concat([np.zeros(_resp_m1.size), _resp_1[_i_roll_tot:]]))
                plt.show()

            _sig_add_resp = _resp_1[_i_roll_tot:]
            _sig_add_data = data[:,epoch_i,_i_roll_tot:]

            resp_linear.append(_sig_add_resp)
            data_linear.append(_sig_add_data)

            _len_resp += _resp_1[_i_roll_tot:].size
            _inspi_start = _len_resp - _half_time 

            if debug:
                plt.plot(_resp_m1)
                plt.vlines(inspi_starts[epoch_i-1], ymin=_resp_m1.min(), ymax=_resp_m1.max(), color='b')
                plt.plot(np.concat([np.zeros(_i_roll), _resp_1]))
                plt.vlines(_inspi_start, ymin=_resp_m1.min(), ymax=_resp_m1.max(), color='r')
                plt.show()

            inspi_starts.append(_inspi_start)
            
        if debug:

            plt.plot(resp_linear)
            plt.show()

    resp_linear = np.concatenate(resp_linear, axis=0)
    data_linear = np.concatenate(data_linear, axis=1)
    inspi_starts = np.array(inspi_starts)

    if len(nan_epochs) != 0:
        df_resp_cycle = df_resp_cycle.drop(index=nan_epochs)

    if np.isnan(resp_linear).sum() != 0 or np.isnan(data_linear).sum() != 0:
        raise ValueError('NAN IN DATA')

    if debug:

        plt.plot(resp_linear)
        plt.vlines(inspi_starts, ymin=resp_linear.min(), ymax=resp_linear.max(), color='r')
        plt.show()

        np.where(np.isnan(data_linear[0]))
        plt.plot(data_linear[0])
        plt.show()
    
    #### select wavelet parameters
    wavelets = get_wavelets()

    #### compute
    print('CONV', flush=True)

    os.chdir(path_memmap)
    tf_allconv = np.memmap(f'memmap_{sujet}_tf_conv.npy', mode="w+", shape=(len(chanlist), nfrex, data_linear.shape[-1]), dtype=np.float32)
            
    #chan_i = 0
    # for chan_i, _ in enumerate(chanlist):
    def extract_chan_conv(chan_i):

        print_advancement(chan_i, len(chanlist), steps=[25, 50, 75])

        for fi in range(nfrex):
            
            tf_allconv[chan_i, fi] = abs(scipy.signal.fftconvolve(data_linear[chan_i], wavelets[fi], 'same'))**2

    joblib.Parallel(n_jobs = n_core, prefer = 'processes')(joblib.delayed(extract_chan_conv)(chan_i) for chan_i, _ in enumerate(chanlist))

    #### extract and norm
    print(f'NORM PARAMS', flush=True)

    if norm_param == 'rscore':

        os.chdir(path_memmap)
        rscore_param = np.memmap(f'memmap_{sujet}_rscore_param.npy', mode="w+", shape=(2, len(chanlist), nfrex), dtype=np.float32)

        # for chan_i, chan in enumerate(chanlist):
        def extract_chan_rscore_params(chan_i): 

            print_advancement(chan_i, data.shape[0], steps=[25, 50, 75])

            for fi in range(nfrex):
                
                _med = np.median(tf_allconv[chan_i, fi,:])
                _mad = scipy.stats.median_abs_deviation(tf_allconv[chan_i, fi,:].reshape(-1), axis=0)

                if debug:

                    _cycles = tf_allconv[chan_i, fi,:]

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
        tf_allconv_norm = np.memmap(f'memmap_{sujet}_tf_conv_norm.npy', mode="w+", shape=(len(chanlist), nfrex, data_linear.shape[-1]), dtype=np.float32)

        def norm_chan_conv(chan_i):

            print_advancement(chan_i, len(chanlist), steps=[25, 50, 75])

            mat = tf_allconv[chan_i]

            mad = rscore_param[0, chan_i, :]
            med = rscore_param[1, chan_i, :]

            mat_norm = (mat - med[:, None]) * 0.6745 / mad[:, None]

            if debug:

                time_plot = 120*srate

                plt.pcolormesh(mat[:,:time_plot])
                plt.colorbar()
                plt.show()

                plt.plot(np.median(mat[:10], axis=0))
                plt.show()

                vmin, vmax = np.percentile(mat_norm[:,:time_plot].reshape(-1), 1), np.percentile(mat_norm[:,:time_plot].reshape(-1), 99)
                plt.pcolormesh(mat_norm[:,:time_plot], vmin=vmin, vmax=vmax)
                plt.colorbar()
                plt.show()

                plt.hist(mat_norm.ravel(), bins=200)
                plt.show()

            tf_allconv_norm[chan_i] = mat_norm.astype(np.float32, copy=False)

        joblib.Parallel(n_jobs=n_core, prefer="threads")(joblib.delayed(norm_chan_conv)(chan_i) for chan_i in range(len(chanlist)))

    if norm_param == None:

        print(f'NO NORM', flush=True)

        tf_allconv_norm[:] = tf_allconv

    if debug:

        for chan_i in range(10):
            _tf = tf_allconv_norm[chan_i]
            min, max = np.percentile(_tf, 1), np.percentile(_tf, 99)
            plt.pcolormesh(_tf, vmin=min, vmax=max)
            plt.colorbar()
            plt.show()

        plt.hist(_tf.reshape(-1))
        plt.show()

        _tf.reshape(-1).max()

    #### construct dataset
    print('CONSTRUCT CYCLES STRETCH')

    _half_time = np.where(np.isclose(time_vec_sec, 0, atol=1e-3))[0][0]

    expi_starts = []

    for _inspi_i, _inspi_val in enumerate(inspi_starts):

        # if _inspi_i == 14:
        #     raise

        _dummy_cycles = np.array([[np.abs(respi_cycle_search_time[0]*srate), np.abs(respi_cycle_search_time[0]*srate) + 0.5*srate, np.abs(respi_cycle_search_time[0]*srate) + 1*srate]])

        _resp = np.array(resp_linear[_inspi_val+respi_cycle_search_time[0]*srate:_inspi_val+respi_cycle_search_time[1]*srate], copy=True)
        _resp_clean = physio.preprocess(_resp, srate, band=25., btype='lowpass', ftype='bessel', order=5, normalize=False)
        _resp_clean_smooth = physio.smooth_signal(_resp_clean, srate, win_shape='gaussian', sigma_ms=40.0)

        if debug:
            plt.plot(_resp_clean_smooth)
            plt.show()
        
        try:
            _cycles = physio.detect_respiration_cycles(_resp_clean_smooth, srate, baseline_mode='median')
        except:
            _cycles = _dummy_cycles

        if debug:

            fig_verif, ax = plt.subplots(figsize=(18, 10))
            ax.plot(_resp)
            ax.scatter(_cycles[:,0], _resp[_cycles[:,0]], color='g', label='inspi_selected')
            ax.scatter(_cycles[:,-1], _resp[_cycles[:,-1]], color='g', label='inspi_selected')
            ax.scatter(_cycles[:,1], _resp[_cycles[:,1]], color='c', label='expi_selected', marker='s')
            plt.vlines([np.abs(respi_cycle_search_time[0]*srate), inspi_starts[_inspi_i+1]-_inspi_val+np.abs(respi_cycle_search_time[0]*srate)], ymin=_resp.min(), ymax=_resp.max(), color='r')
            ax.title('RAW')
            plt.legend()
            fig_verif.show()   

        if _cycles.shape[0] > 1:
            _inspi_sel_clean = _cycles[:,1] > np.abs(respi_cycle_search_time[0]*srate)
            if _inspi_sel_clean.sum() != 0:
                _cycle_clean = _cycles[_inspi_sel_clean]
            else:
                raise
        else:
            if _cycles[:,1] > np.abs(respi_cycle_search_time[0]*srate):
                _cycle_clean = _cycles
            else:
                _cycle_clean = _dummy_cycles

        _expi = _inspi_val+respi_cycle_search_time[0]*srate + _cycle_clean[0,1]

        if _inspi_i != inspi_starts.shape[0]-1:
            if _expi > inspi_starts[_inspi_i+1]:
                _expi = inspi_starts[_inspi_i] + (inspi_starts[_inspi_i+1] - inspi_starts[_inspi_i])/2

        if _inspi_i != 0:
            if expi_starts[_inspi_i-1] > _inspi_val:
                raise

        if _inspi_i != 0 and _inspi_val <= inspi_starts[_inspi_i-1]:
            raise

        if debug:

            fig_verif, ax = plt.subplots(figsize=(18, 10))
            ax.plot(_resp)
            ax.scatter(_cycles[:,0], _resp[_cycles[:,0]], color='g', label='inspi_selected')
            ax.scatter(_cycles[:,-1], _resp[_cycles[:,-1]], color='g', label='inspi_selected')
            ax.scatter(_cycles[:,1], _resp[_cycles[:,1]], color='c', label='expi_selected', marker='s')
            plt.vlines([np.abs(respi_cycle_search_time[0]*srate), inspi_starts[_inspi_i+1]-_inspi_val+np.abs(respi_cycle_search_time[0]*srate)], ymin=_resp.min(), ymax=_resp.max(), color='r')
            ax.title('RAW')
            plt.legend()
            fig_verif.show()   

            fig_verif, ax = plt.subplots(figsize=(18, 10))
            ax.plot(_resp)
            ax.scatter([np.abs(respi_cycle_search_time[0]*srate)], _resp[np.abs(respi_cycle_search_time[0]*srate)], color='g', label='inspi_selected')
            ax.scatter([int(np.abs(respi_cycle_search_time[0]*srate)+_expi-_inspi_val)], _resp[int(np.abs(respi_cycle_search_time[0]*srate)+_expi-_inspi_val)], color='r', label='expi_selected')
            ax.title('SELECTED')
            plt.legend()
            fig_verif.show()   

            plt.plot(resp_linear[:expi_starts[-1]+10*srate])
            plt.scatter(inspi_starts[:len(expi_starts)+1], resp_linear[inspi_starts[:len(expi_starts)+1]], label='inspi', color='g')
            plt.scatter(expi_starts, resp_linear[expi_starts], label='expi', color='r')
            plt.legend()
            plt.show()

        expi_starts.append(_expi)

    expi_starts = np.array(expi_starts)

    if debug :

        plt.plot(resp_linear)
        plt.vlines(inspi_starts, ymin=resp_linear.min(), ymax=resp_linear.max(), color='r')
        plt.vlines(expi_starts, ymin=resp_linear.min(), ymax=resp_linear.max(), color='b')
        plt.show()

    next_inspi = np.append(inspi_starts, resp_linear.size-1)[1:]
    cycles = np.concatenate([inspi_starts.reshape(-1,1), expi_starts.reshape(-1,1), next_inspi.reshape(-1,1)], axis=1)

    if debug:

        plt.scatter(range(cycles[:,0].size), cycles[:,0], label='inspi')
        plt.scatter(range(cycles[:,0].size), cycles[:,1], label='expi')
        plt.scatter(range(cycles[:,0].size), cycles[:,2], label='next_inspi')
        plt.legend()
        plt.show()
        
    #### stretch
    print('STRETCH')

    resp_features = physio.compute_respiration_cycle_features(resp_linear, srate, cycles.astype('int'), baseline=None)    

    tf_allchan_stretch = []

    for chan_i, _ in enumerate(chan_list):

        print_advancement(chan_i, len(chanlist), [25,50,75])

        tf_allchan_stretch.append(stretch_data_tf(resp_features, stretch_point_TF, tf_allconv_norm[chan_i], srate)[0].astype(np.float32))

    tf_allchan_stretch = np.array(tf_allchan_stretch)

    resp_stretch = stretch_data(resp_features, stretch_point_TF, resp_linear, srate)[0]

    if debug:
        for cycle_i in range(resp_stretch.shape[0]):
            plt.plot(resp_stretch[cycle_i])
        plt.show()

        plt.pcolormesh(np.median(tf_allchan_stretch[0], axis=[0]))
        plt.show()

    #### extract pre/post/condition
    if debug:
        plt.hist(resp_features['cycle_duration'], bins=50)
        plt.show()

        plt.hist(resp_features['inspi_duration'], bins=50)
        plt.show()

    extract_good_cycle_jules = (resp_features['cycle_duration'] >= exclude_cycle_duration_jules).values & (resp_features['inspi_duration'] >= exclude_inspi_duration_jules).values
    extract_good_cycle_sam = (df_resp_cycle['isGoodTrial'] == 1).values & (df_resp_cycle['isGoodBreath'] == 1).values 
    extract_good_all = extract_good_cycle_jules & extract_good_cycle_sam

    if debug:
        plt.plot(resp_linear)
        plt.vlines(resp_features[extract_good_all]['inspi_index'], ymin=resp_linear.min(), ymax=resp_linear.max(), color='g')
        plt.vlines(resp_features[extract_good_all]['expi_index'], ymin=resp_linear.min(), ymax=resp_linear.max(), color='r')
        plt.show()

    df_resp_cycle_cleaned = df_resp_cycle[extract_good_all]
    tf_allchan_stretch_cleaned = tf_allchan_stretch[:,extract_good_all]
    resp_stretch_cleaned = resp_stretch[extract_good_all]
    resp_features_cleaned = resp_features[extract_good_all]

    cond_sel_dict = {}
    cond_sel_dict_pre = {}

    try:

        for cond in conditions:
            if cond == 'rsp_ctrl':
                _cond_sel = (df_resp_cycle_cleaned['isControl'] == 1).values & (df_resp_cycle_cleaned['occlusionType'] == 0).values & (df_resp_cycle_cleaned['challengeO2conc'] == 21).values 
            elif cond == 'rsp_chl':
                _cond_sel = (df_resp_cycle_cleaned['isControl'] == 0).values & (df_resp_cycle_cleaned['isChallenge'] == 1).values & (df_resp_cycle_cleaned['occlusionType'] == 0).values & (df_resp_cycle_cleaned['challengeO2conc'] == 21).values 
            elif cond == 'oc_ctrl':
                _cond_sel = (df_resp_cycle_cleaned['isControl'] == 1).values & (df_resp_cycle_cleaned['occlusionType'] == 2).values & (df_resp_cycle_cleaned['challengeO2conc'] == 21).values 
            elif cond == 'oc_chl':
                _cond_sel = (df_resp_cycle_cleaned['isControl'] == 0).values & (df_resp_cycle_cleaned['isChallenge'] == 1).values & (df_resp_cycle_cleaned['occlusionType'] == 2).values & (df_resp_cycle_cleaned['challengeO2conc'] == 21).values 

            cond_sel_dict[cond] = np.where(_cond_sel)[0]
            cond_sel_dict_pre[cond] = cond_sel_dict[cond] -1

    except:

        for cond in conditions:
            if cond == 'rsp_ctrl':
                _cond_sel = (df_resp_cycle_cleaned['isControl'] == 1).values & (df_resp_cycle_cleaned['occlusionType'] == 0).values
            elif cond == 'rsp_chl':
                _cond_sel = (df_resp_cycle_cleaned['isControl'] == 0).values & (df_resp_cycle_cleaned['isChallenge'] == 1).values & (df_resp_cycle_cleaned['occlusionType'] == 0).values
            elif cond == 'oc_ctrl':
                _cond_sel = (df_resp_cycle_cleaned['isControl'] == 1).values & (df_resp_cycle_cleaned['occlusionType'] == 2).values
            elif cond == 'oc_chl':
                _cond_sel = (df_resp_cycle_cleaned['isControl'] == 0).values & (df_resp_cycle_cleaned['isChallenge'] == 1).values & (df_resp_cycle_cleaned['occlusionType'] == 2).values 

            cond_sel_dict[cond] = np.where(_cond_sel)[0]
            cond_sel_dict_pre[cond] = cond_sel_dict[cond] -1

    if debug:

        for cond in conditions:

            resp_plot = resp_stretch_cleaned[cond_sel_dict[cond]]
            for cycle_i in range(resp_plot.shape[0]):
                plt.plot(resp_plot[cycle_i])
            plt.title(f"{cond}, n_cycle:{resp_plot.shape[0]}")
            plt.show()

            chan_i = 0
            tf_plot = tf_allchan_stretch_cleaned[chan_i,cond_sel_dict[cond]]
            plt.pcolormesh(np.median(tf_plot, axis=0))
            plt.title(f"{cond}, n_cycle:{resp_plot.shape[0]}")
            plt.show()
    
    #### save stretch tf
    print(f'SAVE TF STRETCH', flush=True)
    os.chdir(os.path.join(path_precompute, 'TF', 'session'))
    for cond in conditions:
        np.save(f'{sujet}_{cond}_tf_allchan_stretch_post.npy', tf_allchan_stretch_cleaned[:,cond_sel_dict[cond]])    
        np.save(f'{sujet}_{cond}_tf_allchan_stretch_pre.npy', tf_allchan_stretch_cleaned[:,cond_sel_dict_pre[cond]])    

    os.chdir(os.path.join(path_precompute, 'RESP', 'session')) 
    for cond in conditions:
        np.save(f'{sujet}_{cond}_stretch_resp_post.npy', resp_stretch_cleaned[cond_sel_dict[cond]])
        np.save(f'{sujet}_{cond}_stretch_resp_pre.npy', resp_stretch_cleaned[cond_sel_dict_pre[cond]])

    os.chdir(os.path.join(path_precompute, 'RESP', 'respfeatures')) 
    resp_features_cleaned.to_excel(f'{sujet}_respfeatures_cleaned.xlsx')
    df_resp_cycle_cleaned.to_excel(f"{sujet}_cycles_info_cleaned.xlsx")

    #### extract power
    print("EXTRACT POWER", flush=True)

    idx = np.arange(stretch_point_TF)
    inspi_sel = idx <= stretch_point_TF / 2
    expi_sel  = idx >  stretch_point_TF / 2

    bands = list(freq_band_dict.keys())
    pre_posts = ["pre", "post"]
    phases = ["inspi", "expi"]
    chans = list(chan_list)
    conds = list(conditions)

    band_frex_sel = {band: (frex >= freq[0]) & (frex <= freq[-1])
                    for band, freq in freq_band_dict.items()}

    max_cycles = np.max([len(cond_sel_dict[c]) for c in conds])

    Pxx = np.full( (len(conds), len(chans), len(bands), len(pre_posts), len(phases), max_cycles), np.nan, dtype=np.float32 )

    cond_to_i = {c: i for i, c in enumerate(conds)}
    band_to_i = {b: i for i, b in enumerate(bands)}
    prepost_to_i = {"pre": 0, "post": 1}
    phase_to_i = {"inspi": 0, "expi": 1}

    for chan_i, chan_name in enumerate(chans):

        print_advancement(chan_i, len(chans), steps=[25, 50, 75])

        for cond in conds:

            c_i = cond_to_i[cond]

            for cycle_i, cycle_sel_df in enumerate(cond_sel_dict[cond]):

                # indices in your TF array
                pre_post_indices = {"pre": cycle_sel_df - 1, "post": cycle_sel_df}

                for pre_post_sel, pre_post_sel_i in pre_post_indices.items():

                    pp_i = prepost_to_i[pre_post_sel]

                    _tf = tf_allchan_stretch_cleaned[chan_i, pre_post_sel_i] 

                    for band, frex_sel in band_frex_sel.items():
                        
                        b_i = band_to_i[band]

                        tf_band = _tf[frex_sel, :]  

                        pxx_inspi = np.median(tf_band[:, inspi_sel])
                        pxx_expi  = np.median(tf_band[:, expi_sel])

                        Pxx[c_i, chan_i, b_i, pp_i, phase_to_i["inspi"], cycle_i] = pxx_inspi
                        Pxx[c_i, chan_i, b_i, pp_i, phase_to_i["expi"],  cycle_i] = pxx_expi

    da_Pxx = xr.DataArray(
        Pxx,
        dims=("cond", "chan", "band", "pre_post", "phase", "cycle"),
        coords={
            "cond": conds,
            "chan": chans,
            "band": bands,
            "pre_post": pre_posts,
            "phase": phases,
            "cycle": np.arange(max_cycles),
            "ROI": ("chan", localist),   
            "sujet": sujet               
        },
        name="Pxx"
    )


    #### SAVE
    print('SAVE EXTRACT PXX', flush=True)
    os.chdir(os.path.join(path_precompute, 'TF', 'session'))
    da_Pxx.to_netcdf(f'{sujet}_xr_Pxx.nc')    

    #### remove
    os.chdir(path_memmap)
    os.remove(f'memmap_{sujet}_tf_conv.npy')
    os.remove(f'memmap_{sujet}_rscore_param.npy')
    os.remove(f'memmap_{sujet}_tf_conv_norm.npy')

    del data, resp, tf_allconv, tf_allconv_norm

    print('done', flush=True)






#sujet, norm_param = sujet_list[0], 'rscore'
def precompute_tf_allconv_allsession(sujet, norm_param):

    #### verify if already computed
    os.chdir(os.path.join(path_precompute, 'TF', 'session'))

    if os.path.exists(f'{sujet}_df_Pxx.xlsx'):
        print(f'{sujet} ALREADY COMPUTED', flush=True)
        return

    print(f'TF PRECOMPUTE {sujet}', flush=True)

    #### get params
    data, respi, inspi_starts_i, expi_starts_i, df_resp_cycle = get_data_sujet_wholesession(sujet)
    chanlist, localist = get_chanlist(sujet)
    n_trial = len(respi)

    n_cycle_vec = []
    for trial_i in range(n_trial):
        n_cycle_vec.append(len(inspi_starts_i[trial_i]))
    n_cycle_vec = np.array(n_cycle_vec)
    n_cycle = n_cycle_vec.sum()

    phase_vec = np.arange(stretch_point_TF)
    nchan = data[0].shape[0]

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
    os.chdir(path_memmap)
    tf_allconv_stretch = np.memmap(f'memmap_{sujet}_tf_conv_stretch.npy', mode="w+", shape=(nchan, n_cycle, nfrex, phase_vec.shape[0]), dtype=np.float32)

    for trial_i in range(n_trial):

        _data_trial, _respi_trial = data[trial_i], respi[trial_i]
        _inspi_start_trial, _expi_starts_trial = inspi_starts_i[trial_i], expi_starts_i[trial_i]

        if debug:

            plt.plot(_respi_trial)
            plt.vlines(_inspi_start_trial, ymin=_respi_trial.min(), ymax=_respi_trial.max(), color='r')
            plt.vlines(_expi_starts_trial, ymin=_respi_trial.min(), ymax=_respi_trial.max(), color='b')
            plt.show()

        #### conv
        print(f"TRIAL {trial_i+1} conv")

        os.chdir(path_memmap)
        tf_allconv_trial = np.memmap(f'memmap_{sujet}_tf_conv_trial.npy', mode="w+", shape=(nchan, nfrex, _data_trial.shape[-1]), dtype=np.float32)
        #chan_i = 0
        # for chan_i, _ in enumerate(chanlist):
        def extract_chan_conv(chan_i):

            # print_advancement(chan_i, nchan, steps=[25, 50, 75])

            for fi in range(nfrex):
                
                tf_allconv_trial[chan_i,fi,:] = abs(scipy.signal.fftconvolve(_data_trial[chan_i], wavelets[fi,:], 'same'))**2

        joblib.Parallel(n_jobs = n_core, prefer = 'processes')(joblib.delayed(extract_chan_conv)(chan_i) for chan_i in range(nchan))

        #### extract and norm
        print(f"TRIAL {trial_i+1} norm")

        os.chdir(path_memmap)
        rscore_param = np.memmap(f'memmap_{sujet}_rscore_param.npy', mode="w+", shape=(2, nchan, nfrex), dtype=np.float32)

        if norm_param == 'rscore':

            # for chan_i, chan in enumerate(chanlist):
            def extract_chan_rscore_params(chan_i): 

                # print_advancement(chan_i, nchan, steps=[25, 50, 75])

                for fi in range(nfrex):
                    
                    _med = np.median(tf_allconv_trial[chan_i,fi])
                    _mad = scipy.stats.median_abs_deviation(tf_allconv_trial[chan_i,fi], axis=0)

                    if debug:

                        _cycles = tf_allconv_trial[chan_i, :, fi,:]

                        for cycle_i in range(_cycles.shape[0]):
                            plt.plot(_cycles[cycle_i])
                        plt.show()

                        min, max = np.percentile(_cycles, 1), np.percentile(_cycles, 99)
                        plt.pcolormesh(_cycles, vmin=min, vmax=max)
                        plt.colorbar()
                        plt.show()

                    rscore_param[0,chan_i,fi] = _mad
                    rscore_param[1,chan_i,fi] = _med

            joblib.Parallel(n_jobs = n_core, prefer = 'processes')(joblib.delayed(extract_chan_rscore_params)(chan_i) for chan_i in range(nchan))

            os.chdir(path_memmap)
            tf_allconv_norm = np.memmap(f'memmap_{sujet}_tf_conv_norm.npy', mode="w+", shape=(nchan, nfrex, _data_trial.shape[-1]), dtype=np.float32)

            # for chan_i, _ in enumerate(chanlist):
            def norm_chan_conv(chan_i):

                # print_advancement(chan_i, nchan, steps=[25, 50, 75])

                _mat = tf_allconv_trial[chan_i]
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
                
                tf_allconv_norm[chan_i] = _mat_norm

            joblib.Parallel(n_jobs = n_core, prefer = 'processes')(joblib.delayed(norm_chan_conv)(chan_i) for chan_i in range(nchan))

        if norm_param == None:

            print(f'NORM', flush=True)

            tf_allconv_norm[:] = tf_allconv_stretch

        if debug:

            for chan_i in range(10):
                _tf = tf_allconv_norm[chan_i]
                min, max = np.percentile(_tf, 1), np.percentile(_tf, 99)
                plt.pcolormesh(_tf, vmin=min, vmax=max)
                plt.colorbar()
                plt.show()

            plt.hist(_tf.reshape(-1))
            plt.show()

            _tf.reshape(-1).max()
            
        #### stretch
        print(f"TRIAL {trial_i+1} stretch")
    
        _cycle_trial = np.concatenate([_inspi_start_trial.reshape(-1,1), _expi_starts_trial.reshape(-1,1)], axis=1)
        _next_inspi = np.append(_inspi_start_trial[1:], np.array([_respi_trial.shape[0]-1]), axis=0).reshape(-1,1)
        _cycle_trial = np.concatenate([_cycle_trial, _next_inspi], axis=1)

        _resp_features = physio.compute_respiration_cycle_features(_respi_trial, srate, _cycle_trial, baseline=None)    

        tf_allchan_stretch = []

        for chan_i in range(nchan):

            tf_allchan_stretch.append(stretch_data_tf(_resp_features, stretch_point_TF, tf_allconv_norm[chan_i], srate)[0].astype(np.float32))

        if trial_i == 0:
            _insert_start, _insert_stop = 0, n_cycle_vec[trial_i]
        else:
            _insert_start, _insert_stop = n_cycle_vec[:trial_i].sum(), n_cycle_vec[:trial_i].sum() + n_cycle_vec[trial_i]

        tf_allchan_stretch = np.array(tf_allchan_stretch)
        tf_allconv_stretch[:,_insert_start:_insert_stop] = tf_allchan_stretch

        if debug:
            chan_i = 0
            for cycle_i in range(tf_allchan_stretch[chan_i].shape[0]):
                plt.pcolormesh(tf_allchan_stretch[chan_i,cycle_i])
                plt.show()

            plt.pcolormesh(np.median(np.array(tf_allchan_stretch[chan_i]), axis=[0]))
            plt.show()

    #### extract condition
    extract_cond = {}

    df_resp_cycle_corr = df_resp_cycle.query(f"isGoodTrial == 1 and isGoodBreath == 1 and challengeO2conc == 21")

    for cond in conditions:

        if cond == 'rsp_ctrl':
            _sel_cycles = df_resp_cycle_corr.query(f"isControl == 1 and isChallenge == 0 and occlusionType in [0, 1] and challengeLoadMagExh == 0 and challengeLoadMagInh == 0").index

        elif cond == 'rsp_chl':
            _sel_cycles = df_resp_cycle_corr.query(f"isControl == 0 and isChallenge == 1 and occlusionType in [0, 1] and challengeLoadMagExh == 0 and challengeLoadMagInh == 0").index

        elif cond == 'oc_ctrl':
            _sel_cycles = df_resp_cycle_corr.query(f"isControl == 1 and isChallenge == 0 and occlusionType == 2").index

        elif cond == 'oc_chl':
            _sel_cycles = df_resp_cycle_corr.query(f"isControl == 0 and isChallenge == 1 and occlusionType == 2").index

        extract_cond[cond] = _sel_cycles

    #### save stretch tf
    print(f'SAVE TF STRETCH', flush=True)
    os.chdir(os.path.join(path_precompute, 'TF', 'session'))
        
    for cond in conditions:
        np.save(f'{sujet}_{cond}_tf_allchan_stretch.npy', tf_allconv_stretch[:,extract_cond[cond]])    

    #### extract power
    print(f'EXTRACT POWER', flush=True)

    df_Pxx = pd.DataFrame()
    phase_list = ['inspi', 'expi']

    for cond in conditions:

        print(cond)

        _tf_cond = tf_allconv_stretch[:,extract_cond[cond]]

        for chan_i in range(_tf_cond.shape[0]):
        # def extract_Pxx_chan(chan_i):

            # print_advancement(chan_i, _tf_cond.shape[0], steps=[25, 50, 75])

            # _df_chan = pd.DataFrame()

            chan_name = chanlist[chan_i]
            chan_loca = localist[chan_i]

            for cycle_i in range(_tf_cond.shape[1]):

                #### pre
                _tf = _tf_cond[chan_i, cycle_i]

                for band, freq in freq_band_dict.items():

                    frex_sel = (frex >= freq[0]) & (frex <= freq[-1])

                    inspi_sel = np.arange(stretch_point_TF) <= stretch_point_TF/2
                    expi_sel = np.arange(stretch_point_TF) > stretch_point_TF/2
                    # Pxx_pre_inspi_sel = np.median(_tf[frex_sel,:][:,inspi_sel])
                    # Pxx_pre_expi_sel = np.median(_tf[frex_sel,:][:,expi_sel])
                    Pxx_inspi_sel = np.median(_tf[frex_sel,:][:,inspi_sel])
                    Pxx_expi_sel = np.median(_tf[frex_sel,:][:,expi_sel])

                    _df =   pd.DataFrame({'sujet' : [sujet]*len(phase_list), 'cond' : [cond]*len(phase_list), 'chan' : [chan_name]*len(phase_list), 'ROI' : [chan_loca]*len(phase_list),
                                        'band' : [band]*len(phase_list), 'phase' : phase_list, 'cycle' : [cycle_i]*len(phase_list), 'Pxx' : [Pxx_inspi_sel, Pxx_expi_sel]})
                    
                    df_Pxx = pd.concat([df_Pxx, _df])
                    # _df_chan = pd.concat([_df_chan, _df])

            # for epoch_i in range(tf_allchan_stretch_post.shape[1]):

            #     #### post
            #     _tf = tf_allchan_stretch_post[chan_i, epoch_i]

            #     for band, freq in freq_band_dict.items():

            #         frex_sel = (frex >= freq[0]) & (frex <= freq[-1])

            #         inspi_sel = np.arange(stretch_point_TF) <= stretch_point_TF/2
            #         expi_sel = np.arange(stretch_point_TF) > stretch_point_TF/2
            #         Pxx_post_inspi_sel = np.median(_tf[frex_sel,:][:,inspi_sel])
            #         Pxx_post_expi_sel = np.median(_tf[frex_sel,:][:,expi_sel])

            #         _df =   pd.DataFrame({'sujet' : [sujet]*len(post_phase_list), 'cond' : [cond]*len(post_phase_list), 'chan' : [chan_name]*len(post_phase_list), 'ROI' : [chan_loca]*len(post_phase_list),
            #                             'band' : [band]*len(post_phase_list), 'phase' : post_phase_list, 'cycle' : [epoch_i]*len(post_phase_list), 'Pxx' : [Pxx_post_inspi_sel, Pxx_post_expi_sel]})
                    
            #         df_Pxx = pd.concat([df_Pxx, _df])
                    # _df_chan = pd.concat([_df_chan, _df])

            # return _df_chan

    # df_allchan = joblib.Parallel(n_jobs = n_core, prefer = 'processes')(joblib.delayed(extract_Pxx_chan)(chan_i) for chan_i, _ in enumerate(chanlist))

    # for chan_i, _ in enumerate(chanlist):

    #     extract_Pxx_chan(chan_i)

    #### SAVE
    print('SAVE EXTRACT PXX', flush=True)
    os.chdir(os.path.join(path_precompute, 'TF', 'session'))

    df_Pxx.to_excel(f'{sujet}_df_Pxx.xlsx')    

    #### remove
    os.chdir(path_memmap)
    os.remove(f'memmap_{sujet}_tf_conv_stretch.npy')
    os.remove(f'memmap_{sujet}_tf_conv_trial.npy')
    os.remove(f'memmap_{sujet}_rscore_param.npy')
    os.remove(f'memmap_{sujet}_tf_conv_norm.npy')

    del data, respi, tf_allconv_stretch, tf_allconv_trial, rscore_param, tf_allconv_norm

    print('done', flush=True)






################################
######## EXECUTE ########
################################


if __name__ == '__main__':

    # param_list = []
    norm_param = 'rscore'
    
    #sujet = sujet_list[0]
    for sujet in sujet_list:

        # param_list.append([sujet, norm_param])
        precompute_tf_allconv_from_epoch(sujet, norm_param)

    # joblib.Parallel(n_jobs = n_core, prefer = 'processes')(joblib.delayed(precompute_tf_allconv)(sujet, norm_param) for sujet in param_list[:3])



    



                        