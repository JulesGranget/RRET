

import joblib
import plotly

from A_config.n01_O_params import *
from A_config.n02_O_analysis_functions import *
from A_config.n03_X_manip_data import *






################################
######## PRECOMPUTE TF ########
################################

#sujet = sujet_list[0]
def extract_press_AF_resp(sujet):

    #### verify if already computed
    os.chdir(os.path.join(path_precompute, 'TF', 'session'))

    # if os.path.exists(f"{sujet}_rsp_ctrl_tf_allchan_stretch_post.npy"):
    #     print(f'{sujet} ALREADY COMPUTED', flush=True)
    #     return

    print(f'RESP {sujet}', flush=True)

    #### get params
    data, press, AF, co2, chanlist, localist = get_data_sujet_resp_fullsig(sujet)

    time_vec_sec = np.arange(press[0].size)/srate - (press[0].size/srate)/2
    chan_list, loca_list = get_chanlist(sujet)

    _half_time = np.where(np.isclose(time_vec_sec, 0, atol=1e-3))[0][0]

    os.chdir(os.path.join(path_prep))
    df_resp_cycle = pd.read_excel(f"{sujet}_breathInfo.xlsx")

    #### generate linear signals
    print('CONSTRUCT LINEAR SIG')
    press_linear = []
    AF_linear = []
    co2_linear = []
    data_linear = []
    inspi_starts = []
    _len_resp = 0

    nan_epochs = []

    for epoch_i in range(press.shape[0]):

        print_advancement(epoch_i, press.shape[0], [25,50,75])

        if np.isnan(data[:,epoch_i]).sum() != 0 or np.isnan(press[epoch_i]).sum() != 0:
            nan_epochs.append(epoch_i)
            continue

        if debug:
            plt.plot(press[epoch_i])
            plt.show()

            plt.plot(co2[epoch_i])
            plt.show()

            plt.pcolormesh(data[:,epoch_i])
            plt.show()

        if epoch_i == 0:

            press_linear.append(press[epoch_i])
            AF_linear.append(AF[epoch_i])
            co2_linear.append(co2[epoch_i])
            data_linear.append(data[:,epoch_i])
            inspi_starts.append(_half_time)
            _len_resp += press[epoch_i].size

        else:

            _resp_m1 = press[epoch_i-1]
            _resp_1 = press[epoch_i]
            _AF_1 = AF[epoch_i]
            _co2_1 = co2[epoch_i]

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

            _sig_add_press = _resp_1[_i_roll_tot:]
            _sig_add_AF = _AF_1[_i_roll_tot:]
            _sig_add_co2 = _co2_1[_i_roll_tot:]
            _sig_add_data = data[:,epoch_i,_i_roll_tot:]

            press_linear.append(_sig_add_press)
            AF_linear.append(_sig_add_AF)
            co2_linear.append(_sig_add_co2)
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

            plt.plot(press_linear)
            plt.show()

    press_linear = np.concatenate(press_linear, axis=0)
    AF_linear = np.concatenate(AF_linear, axis=0)
    co2_linear = np.concatenate(co2_linear, axis=0)
    data_linear = np.concatenate(data_linear, axis=1)
    inspi_starts = np.array(inspi_starts)

    if len(nan_epochs) != 0:
        df_resp_cycle = df_resp_cycle.drop(index=nan_epochs)

    if np.isnan(press_linear).sum() != 0 or np.isnan(data_linear).sum() != 0:
        raise ValueError('NAN IN DATA')

    #### norm
    press_linear = scipy.stats.zscore(press_linear)
    AF_linear = scipy.stats.zscore(AF_linear)

    if debug:

        plt.plot(press_linear)
        plt.vlines(inspi_starts, ymin=press_linear.min(), ymax=press_linear.max(), color='r')
        plt.show()

        plt.plot(co2_linear)
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

    #### construct dataset
    print('CONSTRUCT CYCLES STRETCH')

    _half_time = np.where(np.isclose(time_vec_sec, 0, atol=1e-3))[0][0]

    expi_starts = []

    for _inspi_i, _inspi_val in enumerate(inspi_starts):

        # if _inspi_i == 14:
        #     raise

        _dummy_cycles = np.array([[np.abs(respi_cycle_search_time[0]*srate), np.abs(respi_cycle_search_time[0]*srate) + 0.5*srate, np.abs(respi_cycle_search_time[0]*srate) + 1*srate]])

        _resp = np.array(press_linear[_inspi_val+respi_cycle_search_time[0]*srate:_inspi_val+respi_cycle_search_time[1]*srate], copy=True)
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

            plt.plot(press_linear[:expi_starts[-1]+10*srate])
            plt.scatter(inspi_starts[:len(expi_starts)+1], press_linear[inspi_starts[:len(expi_starts)+1]], label='inspi', color='g')
            plt.scatter(expi_starts, press_linear[expi_starts], label='expi', color='r')
            plt.legend()
            plt.show()

        expi_starts.append(_expi)

    expi_starts = np.array(expi_starts)

    if debug :

        plt.plot(press_linear)
        plt.vlines(inspi_starts, ymin=press_linear.min(), ymax=press_linear.max(), color='r')
        plt.vlines(expi_starts, ymin=press_linear.min(), ymax=press_linear.max(), color='b')
        plt.show()

    next_inspi = np.append(inspi_starts, press_linear.size-1)[1:]
    cycles = np.concatenate([inspi_starts.reshape(-1,1), expi_starts.reshape(-1,1), next_inspi.reshape(-1,1)], axis=1)

    if debug:

        plt.scatter(range(cycles[:,0].size), cycles[:,0], label='inspi')
        plt.scatter(range(cycles[:,0].size), cycles[:,1], label='expi')
        plt.scatter(range(cycles[:,0].size), cycles[:,2], label='next_inspi')
        plt.legend()
        plt.show()
        
    #### stretch
    print('STRETCH')

    resp_features_press = physio.compute_respiration_cycle_features(press_linear, srate, cycles.astype('int'), baseline=None) 
    resp_features_AF = physio.compute_respiration_cycle_features(AF_linear, srate, cycles.astype('int'), baseline=None)    

    press_stretch = stretch_data(resp_features_press, stretch_point_TF, press_linear, srate)[0]
    AF_stretch = stretch_data(resp_features_AF, stretch_point_TF, AF_linear, srate)[0]
    co2_stretch = stretch_data(resp_features_press, stretch_point_TF, co2_linear, srate)[0]

    if debug:
        for cycle_i in range(press_stretch.shape[0]):
            plt.plot(press_stretch[cycle_i])
        plt.show()

        for cycle_i in range(co2_stretch.shape[0]):
            plt.plot(co2_stretch[cycle_i])
        plt.show()

    #### extract pre/post/condition
    if debug:
        plt.hist(resp_features_press['cycle_duration'], bins=50)
        plt.show()

        plt.hist(resp_features_press['inspi_duration'], bins=50)
        plt.show()

    extract_good_cycle_jules = (resp_features_press['cycle_duration'] >= exclude_cycle_duration_jules).values & (resp_features_press['inspi_duration'] >= exclude_inspi_duration_jules).values
    extract_good_cycle_sam = (df_resp_cycle['isGoodTrial'] == 1).values & (df_resp_cycle['isGoodBreath'] == 1).values 
    extract_good_all = extract_good_cycle_jules & extract_good_cycle_sam

    if debug:
        plt.plot(press_linear)
        plt.vlines(resp_features_press[extract_good_all]['inspi_index'], ymin=press_linear.min(), ymax=press_linear.max(), color='g')
        plt.vlines(resp_features_press[extract_good_all]['expi_index'], ymin=press_linear.min(), ymax=press_linear.max(), color='r')
        plt.show()

    df_resp_cycle_cleaned = df_resp_cycle[extract_good_all]
    press_stretch_cleaned = press_stretch[extract_good_all]
    AF_stretch_cleaned = AF_stretch[extract_good_all]
    co2_stretch_cleaned = co2_stretch[extract_good_all]
    press_respfeatures_cleaned = resp_features_press[extract_good_all]
    AF_respfeatures_cleaned = resp_features_AF[extract_good_all]

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

            resp_plot = press_respfeatures_cleaned[cond_sel_dict[cond]]
            for cycle_i in range(resp_plot.shape[0]):
                plt.plot(resp_plot[cycle_i])
            plt.title(f"{cond}, n_cycle:{resp_plot.shape[0]}")
            plt.show()


    #### add cond label to respfeature
    df_respfeature_label = []

    for cond in conditions:

        _df = press_respfeatures_cleaned.iloc[cond_sel_dict[cond]].copy()
        _df['sujet'] = [sujet] * _df.shape[0]
        _df['cond'] = [cond] * _df.shape[0]
        _df['etCO2'] = co2_stretch_cleaned[cond_sel_dict[cond]].max(axis=1)
        _df['resptype'] = ['press'] * _df.shape[0]
        df_respfeature_label.append(_df)

        _df = AF_respfeatures_cleaned.iloc[cond_sel_dict[cond]].copy()
        _df['sujet'] = [sujet] * _df.shape[0]
        _df['cond'] = [cond] * _df.shape[0]
        _df['etCO2'] = co2_stretch_cleaned[cond_sel_dict[cond]].max(axis=1)
        _df['resptype'] = ['AF'] * _df.shape[0]
        df_respfeature_label.append(_df)

    df_respfeature_label = pd.concat(df_respfeature_label)

    if debug:

        plt.plot(scipy.stats.zscore(press_linear), label='press')
        plt.plot(scipy.stats.zscore(AF_linear), label='AF')
        plt.legend()
        plt.show()

        plt.plot(press_stretch_cleaned.mean(axis=0), label='press')
        plt.plot(AF_stretch_cleaned.mean(axis=0), label='AF')
        plt.legend()
        plt.show()

    #### cond stretch extract
    press_stretch_cleaned_cond, AF_stretch_cleaned_cond = {}, {}
    for cond in conditions:

        press_stretch_cleaned_cond[cond] = press_stretch_cleaned[cond_sel_dict[cond]]
        AF_stretch_cleaned_cond[cond] = AF_stretch_cleaned[cond_sel_dict[cond]]

    return df_respfeature_label, press_linear, press_stretch_cleaned_cond, AF_linear, AF_stretch_cleaned_cond










################################
######## EXECUTE ########
################################


if __name__ == '__main__':

    #### load
    df_respfeature_label, press_linear, press_stretch_cleaned, AF_linear, AF_stretch_cleaned = [], [], [], [], []

    sujet_list_press_AF = ['NS131_02', 'NS203', 'NS205', 'NS206', 'NS208', 'NS211', 'NS213', 'NS215'] #NS202 no eflow
    
    #sujet = sujet_list[1]
    for sujet in sujet_list_press_AF:

        _df_respfeature_label, _press_linear, _press_stretch_cleaned, _AF_linear, _AF_stretch_cleaned = extract_press_AF_resp(sujet)
        df_respfeature_label.append(_df_respfeature_label)
        press_linear.append(_press_linear)
        press_stretch_cleaned.append(_press_stretch_cleaned)
        AF_linear.append(_AF_linear)
        AF_stretch_cleaned.append(_AF_stretch_cleaned)

    df_respfeature_label = pd.concat(df_respfeature_label)

    #### plot
    path_export_fig = os.path.join(path_results, 'respi', 'AF_vs_press')
    df_respfeature_label.to_excel(os.path.join(path_export_fig, f"respfeatures_allsujet_AFvsPRESS.xlsx"))

    metrics = ['cycle_duration', 'inspi_duration', 'expi_duration', 'cycle_freq', 'cycle_ratio', 'inspi_volume',
               'expi_volume', 'total_amplitude', 'total_volume', 'etCO2']

    for metric in metrics:

        sns.catplot(data=df_respfeature_label, kind='strip', x='cond', y=metric, hue='resptype', col='sujet')
        plt.suptitle(metric)
        # plt.show()
        plt.savefig(os.path.join(path_export_fig, f"{metric}_AFvsPRESS.png"))

    plt.close('all')

    fig, axs = plt.subplots(ncols=4)

    for col, cond in enumerate(conditions):

        ax = axs[col]

        for sujet_i, sujet in enumerate(sujet_list_press_AF):

            ax.plot(np.median(press_stretch_cleaned[sujet_i][cond], axis=0), color='b')
            ax.plot(np.median(AF_stretch_cleaned[sujet_i][cond], axis=0), color='r')

        ax.set_title(cond)

    # plt.show()
    fig.savefig(os.path.join(path_export_fig, f"MEAN_respi_AFvsPRESS.png"))
    plt.close('all')

    #### stats
    import pingouin as pg
    
    metrics = ['cycle_duration', 'inspi_duration', 'expi_duration', 'cycle_freq', 'cycle_ratio', 'inspi_volume',
                   'expi_volume', 'total_amplitude', 'total_volume', 'etCO2']

    df_stats_test = df_respfeature_label.query(f"resptype == 'AF'")
    df_stats_resp = []

    for metric in metrics:

        aov = pg.rm_anova(data=df_stats_test[['sujet', 'cond', metric]], dv=metric, within="cond", subject="sujet", detailed=True)
        posthoc = pg.pairwise_tests(data=df_respfeature_label, dv=metric, within="cond", subject="sujet", parametric=False, padjust="holm")
        posthoc['pre_hoc'] = ['rmANOVA'] * posthoc.shape[0]
        posthoc['pre_hoc_pval'] = [aov["p_unc"].iloc[0]] * posthoc.shape[0]
        posthoc['metric'] = [metric] * posthoc.shape[0]
        df_stats_resp.append(posthoc)
    
    df_stats_resp = pd.concat(df_stats_resp)
    df_stats_resp = df_stats_resp[['metric', 'pre_hoc_pval', 'Contrast', 'A', 'B', 'Paired', 'Parametric', 'W_val', 'alternative', 'p_unc', 'p_corr', 'p_adjust', 'hedges', 'pre_hoc']]

    




                        