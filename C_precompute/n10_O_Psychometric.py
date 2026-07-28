


from A_config.n01_O_params import *
from A_config.n02_O_analysis_functions import *
from A_config.n03_X_manip_data import *







################################
######## EXTRACT STRESS ########
################################

def get_df_unpl():

    #### config
    col_to_extract = ['sujet', 'cond', 'trialUnpleasantness', 'trialAnxiety', 'trialAnxietySource']
    path_import_df_psycho = os.path.join(path_precompute, 'PSYCHO', 'df_export')
    path_export_res = os.path.join(path_results, 'Psycho')

    #### load_df
    df_psycho = pd.concat([pd.read_excel(os.path.join(path_import_df_psycho, f'{sujet}_df_psycho_MECACO2BOTH.xlsx')) for sujet in sujet_list_interaction_analysis])

    #### trial num
    df_plot = df_psycho[["sujet", "cond", "trial"]].drop_duplicates().reset_index(drop=True).groupby(['sujet', 'cond']).count().reset_index()

    sns.catplot(data=df_plot, kind='bar', x='sujet', y='trial', hue='cond', hue_order=['ctrl', 'MECA', 'CO2', 'BOTH'])
    # plt.show()
    plt.savefig(os.path.join(path_export_res, f"trial_num.png"))

    #### plot
    df_plot = df_psycho.groupby(['sujet', 'cond', 'trial']).median().reset_index()[col_to_extract]
    df_plot = pd.melt(df_psycho, id_vars=['sujet', 'cond'], value_vars=['trialUnpleasantness', 'trialAnxiety', 'trialAnxietySource'], var_name='psychometric', value_name='val')
        
    sns.catplot(data=df_plot, kind='bar', x='sujet', y='val', hue='cond', col='psychometric', hue_order=['ctrl', 'MECA', 'CO2', 'BOTH'])
    # plt.show()
    plt.savefig(os.path.join(path_export_res, f"psychometric_values.png"))
        
    plt.close('all')        






########################################
######## PSYCHOMETRIC ANALYSIS ########
########################################


def extract_df_psycho_main_analysis(sujet):

    print(sujet)

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

            resp_linear.append(_sig_add_resp)

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
    inspi_starts = np.array(inspi_starts)

    if len(nan_epochs) != 0:
        df_resp_cycle = df_resp_cycle.drop(index=nan_epochs)

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

    cond_sel_dict = {}

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


    #### add cond label to respfeature and save
    df_psycho = []

    for cond in conditions:

        _df_psycho = df_resp_cycle_cleaned.iloc[cond_sel_dict[cond]]
        _df_psycho['sujet'] = [sujet] * _df_psycho.shape[0]
        _df_psycho['cond'] = [cond] * _df_psycho.shape[0]

        df_psycho.append(_df_psycho)

    df_psycho = pd.concat(df_psycho)

    os.chdir(os.path.join(path_precompute, 'PSYCHO', 'df_export')) 
    df_psycho.to_excel(f"{sujet}_df_psycho.xlsx")
    
    


    
    
    
    





################################
######## EXECUTE ########
################################


if __name__ == '__main__':
    
    get_df_unpl()

    for sujet in sujet_list:
    
        extract_df_psycho_main_analysis(sujet)

    


                        