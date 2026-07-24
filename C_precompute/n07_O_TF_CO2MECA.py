

import joblib
import plotly

from A_config.n01_O_params import *
from A_config.n02_O_analysis_functions import *
from A_config.n03_X_manip_data import *









################################
######## PRECOMPUTE TF ########
################################

#sujet, norm_param = sujet_list_interaction_analysis[2], 'rscore'
def precompute_tf_allconv_from_epoch(sujet, norm_param):

    #### verify if already computed
    path_export_TF = os.path.join(path_precompute, 'TF', 'MECACO2_INTER') 

    # if os.path.exists(os.path.join(path_export_TF, f"{sujet}_noc_ctrl_tf_allchan_stretch_MECACO2INTER.npy")):
    #     print(f'{sujet} ALREADY COMPUTED', flush=True)
    #     return

    print(f'TF PRECOMPUTE {sujet}', flush=True)

    #### get params
    data, resp, chanlist, localist = get_data_sujet_fullsig(sujet)

    time_vec_sec = np.arange(resp[0].size)/srate - (resp[0].size/srate)/2
    chan_list, loca_list = get_chanlist(sujet)

    _half_time = np.where(np.isclose(time_vec_sec, 0, atol=1e-3))[0][0]

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

        sig = data_linear[0]
        plt.plot(sig)
        plt.vlines(inspi_starts, ymin=sig.min(), ymax=sig.max(), color='r')
        plt.show()

        correct_sujet_inspi_start = np.array([369787, 371347, 373417,
        375477, 377727, 379957, 382047, 383867, 385992, 388252, 390372,
        392637, 394702, 396822, 398642, 400922, 402927, 408928])
        inspi_starts_corrected = np.array([_inspi for _inspi in inspi_starts if _inspi not in correct_sujet_inspi_start])

        sig = data_linear[0]
        plt.plot(sig)
        plt.vlines(inspi_starts_corrected, ymin=sig.min(), ymax=sig.max(), color='r')
        plt.show()

        np.where(np.isnan(data_linear[0]))
        plt.plot(data_linear[0])
        plt.show()

    #### correct for artifacts
    if sujet == 'NS211':

        correct_sujet_inspi_start = np.array([369787, 371347, 373417,
        375477, 377727, 379957, 382047, 383867, 385992, 388252, 390372,
        392637, 394702, 396822, 398642, 400922, 402927, 408928])
        inspi_starts_corrected = np.array([[_inspi_i, _inspi] for _inspi_i, _inspi in enumerate(inspi_starts) if _inspi not in correct_sujet_inspi_start])

        df_resp_cycle = df_resp_cycle.iloc[inspi_starts_corrected[:,0]]

        inspi_starts = inspi_starts_corrected[:,1]

    
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

    #### free space
    os.chdir(path_memmap)
    os.remove(f'memmap_{sujet}_tf_conv.npy')
    del tf_allconv

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

    #### prepare resp_feature df
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

    ##### keep only O2 concentration 21 
    if 'challengeO2conc' not in df_resp_cycle_cleaned.columns: 
        df_resp_cycle_cleaned['challengeO2conc'] = np.ones(df_resp_cycle_cleaned.shape[0]) * 21

    mask_O2 = (df_resp_cycle_cleaned['challengeO2conc'] == 21).values 
    df_resp_cycle_cleaned = df_resp_cycle_cleaned[mask_O2]
    tf_allchan_stretch_cleaned = tf_allchan_stretch_cleaned[:,mask_O2]
    resp_stretch_cleaned = resp_stretch_cleaned[mask_O2]
    resp_features_cleaned = resp_features_cleaned[mask_O2]

    #### extract conditions
    cond_mask = {}
    oc_list = ['oc', 'noc']

    for oc_cond in oc_list:

        cond_mask[oc_cond] = {}

        if oc_cond == 'oc':
            _cond_sel_OC = (df_resp_cycle_cleaned['occlusionType'] == 2).values
        elif oc_cond == 'noc':
            _cond_sel_OC = (df_resp_cycle_cleaned['occlusionType'] == 0).values

        for cond in cond_list_interaction:

            if cond == 'ctrl':
                _cond_sel = (df_resp_cycle_cleaned['isControl'] == 1).values & (df_resp_cycle_cleaned['challengeLoadMag'] == 0).values & (df_resp_cycle_cleaned['challengeCO2conc'] == 0).values
            elif cond == 'CO2':
                _cond_sel = (df_resp_cycle_cleaned['isControl'] == 0).values & (df_resp_cycle_cleaned['challengeLoadMag'] == 0).values & (df_resp_cycle_cleaned['challengeCO2conc'] != 0).values
            elif cond == 'MECA':
                _cond_sel = (df_resp_cycle_cleaned['isControl'] == 0).values & (df_resp_cycle_cleaned['challengeLoadMag'] != 0).values & (df_resp_cycle_cleaned['challengeCO2conc'] == 0).values
            elif cond == 'BOTH':
                _cond_sel = (df_resp_cycle_cleaned['isControl'] == 0).values & (df_resp_cycle_cleaned['challengeLoadMag'] != 0).values & (df_resp_cycle_cleaned['challengeCO2conc'] != 0).values

            cond_mask[oc_cond][cond] = _cond_sel & _cond_sel_OC

    if debug:

        for cond in cond_list_interaction:

            print(f"{cond} {cond_mask[oc_cond][cond].sum()}")

        for cond in cond_list_interaction:

            resp_plot = resp_stretch_cleaned[cond_mask[oc_cond][cond]]
            for cycle_i in range(resp_plot.shape[0]):
                plt.plot(resp_plot[cycle_i])
            plt.title(f"{cond}, n_cycle:{resp_plot.shape[0]}")
            plt.show()

            chan_i = 0
            tf_plot = tf_allchan_stretch_cleaned[chan_i,cond_mask[oc_cond][cond]]
            plt.pcolormesh(np.median(tf_plot, axis=0))
            plt.title(f"{cond}, n_cycle:{resp_plot.shape[0]}")
            plt.show()

    #### add cond label to respfeature and save
    df_respfeature_label = []
    df_psycho = []
    oc_cond = 'noc'

    for cond in cond_list_interaction:

        _df_psycho = df_resp_cycle_cleaned[cond_mask[oc_cond][cond]][['trial', 'currentTrialCount', 'trialUnpleasantness', 'trialAnxiety', 'trialAnxietySource']]
        _df_psycho['cycle_i'] = np.arange(_df_psycho.shape[0])
        _df_psycho['sujet'] = [sujet] * _df_psycho.shape[0]
        _df_psycho['cond'] = [cond] * _df_psycho.shape[0]
        df_psycho.append(_df_psycho)

        _df = resp_features_cleaned.iloc[cond_mask[oc_cond][cond]].copy()
        _df['sujet'] = [sujet] * _df.shape[0]
        _df['cond'] = [cond] * _df.shape[0]

        df_respfeature_label.append(_df)

    df_respfeature_label = pd.concat(df_respfeature_label)
    df_psycho = pd.concat(df_psycho)

    path_export = os.path.join(path_precompute, 'RESP', 'respfeatures')
    df_respfeature_label.to_excel(os.path.join(path_export, f'{sujet}_respfeatures_cleaned_label_MECACO2INTER.xlsx'))
    df_resp_cycle_cleaned.to_excel(os.path.join(path_export, f"{sujet}_cycles_info_cleaned_MECACO2INTER.xlsx"))

    path_export_psycho = os.path.join(path_precompute, 'PSYCHO', 'df_export')
    df_psycho.to_excel(os.path.join(path_export_psycho, f'{sujet}_df_psycho.xlsx'))
    
    #### save stretch tf
    print(f'SAVE TF STRETCH', flush=True)
    for oc_cond in oc_list:
        for cond in cond_list_interaction:
            np.save(os.path.join(path_export_TF, f'{sujet}_{oc_cond}_{cond}_tf_allchan_stretch_MECACO2INTER.npy'), tf_allchan_stretch_cleaned[:,cond_mask[oc_cond][cond]])

    path_export_resp = os.path.join(path_precompute, 'RESP', 'MECACO2INTER')
    for oc_cond in oc_list:
        for cond in cond_list_interaction:
            np.save(os.path.join(path_export_resp, f'{sujet}_{oc_cond}_{cond}_stretch_resp_MECACO2INTER.npy'), resp_stretch_cleaned[cond_mask[oc_cond][cond]])

    #### remove
    os.chdir(path_memmap)
    os.remove(f'memmap_{sujet}_rscore_param.npy')
    os.remove(f'memmap_{sujet}_tf_conv_norm.npy')

    del data, resp, tf_allconv_norm









################################
######## EXTRACT POWER ########
################################


def extract_power(sujet):

    print("EXTRACT POWER", flush=True)

    #### verify if already computed
    path_export_xr = os.path.join(path_precompute, 'TF', 'MECACO2_INTER')

    if os.path.exists(os.path.join(path_export_xr, f'{sujet}_xr_Pxx_MECACO2INTER.nc')):
        print(f'{sujet} ALREADY COMPUTED', flush=True)
        return

    #### params
    oc_list = ['noc', 'oc']
    chanlist, localist = get_chanlist(sujet)

    idx = np.arange(stretch_point_TF)
    inspi_sel = idx <= stretch_point_TF / 2
    expi_sel  = idx >  stretch_point_TF / 2

    bands = list(freq_band_dict.keys())
    phases = ["inspi", "expi"]

    band_frex_sel = {band: (frex >= freq[0]) & (frex <= freq[-1])
                    for band, freq in freq_band_dict.items()}
    
    #### load
    path_load_data = os.path.join(path_precompute, 'TF', 'MECACO2_INTER')

    tf_allcond = {}

    for oc_cond in oc_list:

        tf_allcond[oc_cond] = {}

        for cond in cond_list_interaction:

            _mat = np.load(os.path.join(path_load_data, f'{sujet}_{oc_cond}_{cond}_tf_allchan_stretch_MECACO2INTER.npy'))
            tf_allcond[oc_cond][cond] = _mat

    #### extract
    max_list = []
    for oc_cond_i, oc_cond in enumerate(oc_list):

        for cond in cond_list_interaction:

            max_list.append(tf_allcond[oc_cond][cond].shape[1])
    
    max_cycles = np.max(np.array(max_list))

    Pxx = np.full( (len(oc_list), len(cond_list_interaction), len(chanlist), len(bands), len(phases), max_cycles), np.nan, dtype=np.float32 )

    for oc_cond_i, oc_cond in enumerate(oc_list):

        print(f"{oc_cond} extract")

        for chan_i, chan_name in enumerate(chanlist):

            print_advancement(chan_i, len(chanlist), steps=[25, 50, 75])

            for cond_i, cond in enumerate(cond_list_interaction):

                ncycle = tf_allcond[oc_cond][cond].shape[1]

                for cycle_i in range(ncycle):

                    for band_i, (band, frex_sel) in enumerate(band_frex_sel.items()):

                        _tf = tf_allcond[oc_cond][cond][chan_i, cycle_i][frex_sel, :]  

                        pxx_inspi = np.median(_tf[:, inspi_sel])
                        pxx_expi  = np.median(_tf[:, expi_sel])

                        Pxx[oc_cond_i, cond_i, chan_i, band_i, 0, cycle_i] = pxx_inspi
                        Pxx[oc_cond_i, cond_i, chan_i, band_i, 1,  cycle_i] = pxx_expi

    da_Pxx = xr.DataArray(
        Pxx,
        dims=("oc_cond", "cond", "chan", "band", "phase", "cycle"),
        coords={
            "oc_cond" : oc_list,
            "cond": cond_list_interaction,
            "chan": chanlist,
            "band": bands,
            "phase": phases,
            "cycle": np.arange(max_cycles),
            "ROI": ("chan", localist),   
            "sujet": sujet               
        },
        name="Pxx"
    )

    #### SAVE
    print('SAVE EXTRACT PXX', flush=True)
    da_Pxx.to_netcdf(os.path.join(path_export_xr, f'{sujet}_xr_Pxx_MECACO2INTER.nc'))   

    print('done', flush=True)






################################
######## EXECUTE ########
################################


if __name__ == '__main__':

    #### PRECOMPUTE
    norm_param = 'rscore'
    
    #sujet = sujet_list_interaction_analysis[0]
    for sujet in sujet_list_interaction_analysis:

        precompute_tf_allconv_from_epoch(sujet, norm_param)
        extract_power(sujet)

    


                        