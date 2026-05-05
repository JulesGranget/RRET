

import joblib
import plotly

from n00_config_O_params import *
from n00bis_config_O_analysis_functions import *
from n00ter_X_manip_data import *







################################
######## EXTRACT ECG ########
################################


#sujet = sujet_list_ecg[0]
def get_resp_ecg_sujet(sujet):

    print(f"open data {sujet}")
    os.chdir(os.path.join(path_prep))

    # time vec : 0 start of the inspiration

    # trial : trial number
    # occlusionType	:  0 no occlusion, 1 control occulusion, 2 occlusion
    # isControl : 1 True, 0 False
    # isChallenge : 1 True, 0 False
    # isRelief : 1 True, 0 False
    # trialUnpleasantness : 0 - 10 after each trial
    # trialAnxiety : 0 - 10 after each trial
    # trialAnxietySource  : 4 breathing effort, lack of air, neither or both
    # isGoodTrial : 0 bad 1 good 	
    # isGoodBreath : 0 bad 1 good

    try:
        data = scipy.io.loadmat(f"RRBO_{sujet}.mat")
        print('open loadmat')
        print(data.keys())
        print(data['RRBO'][0,0].dtype.names)
        data['RRBO'][0,0]['pressure']

        data_epoch = np.transpose(data['RRBO'][0,0]['eeg'][:], (0,2,1)) #chan, breaths, time points 
        resp_epoch = np.transpose(data['RRBO'][0,0]['pressure'][:], (1,0))
        ecg_epoch_raw = np.transpose(data['RRBO'][0,0]['emg'][:], (1,0))
        
        chanlist = np.array([_chan[0][0] for _chan in data['RRBO'][0,0]['labels']])
        localist = np.array([_chan[0][0] for _chan in data['RRBO'][0,0]['labelsFSurf']])

    except:
        # print('open h5py')
        with h5py.File(f"RRBO_{sujet}.mat", "r") as f:

            print(list(f.keys()))

            data = f["RRBO"]
            print(list(data.keys()))

            data_epoch = np.transpose(data['eeg'][:], (2,0,1))
            resp_epoch = data['pressure'][:]
            ecg_epoch_raw = data['emg'][:]

            chanlist = np.array(["".join(chr(c) for c in f[_ref][:].flatten()) for _ref in data['labels'][0]])
            localist = np.array(["".join(chr(c) for c in f[_ref][:].flatten()) for _ref in data['labelsFSurf'][0]])

    ecg_epoch_raw = ecg_epoch_raw[:,:,ecg_sujet_pin_choice[sujet]]
    ecg_epoch = ecg_epoch_raw.reshape(ecg_epoch_raw.shape[0], -1)

    return resp_epoch, ecg_epoch




################################
######## EXPORT DF ECG ########
################################


def export_df_ecg():

    df_ecg = []

    for sujet in sujet_list_ecg:

        resp_epoch, ecg_epoch = get_resp_ecg_sujet(sujet)

        if debug:

            plt.plot(np.concat(ecg_epoch))
            plt.show()

            plt.plot(physio.preprocess(np.concat(ecg_epoch), srate))
            plt.show()

        os.chdir(os.path.join(path_precompute, 'RESP', 'respfeatures')) 
        # df_respfeature_label.to_excel(f'{sujet}_respfeatures_cleaned_label.xlsx')
        # df_respfeature_label_pre.to_excel(f'{sujet}_respfeatures_cleaned_label_pre.xlsx')
        df_resp_cycle_raw = pd.read_excel(f"{sujet}_cycles_info_cleaned.xlsx")
        df_breath_info_cond = df_resp_cycle_raw.query(f"isRelief == 0")

        trial_i_list = df_breath_info_cond['trial'].unique()

        for trial_i_add, trial_i_sel in enumerate(trial_i_list):

            df_trial = df_breath_info_cond.query(f"trial == {trial_i_sel}")

            if df_trial['isChallenge'].sum() != 0:
                trial_cond = 'CHL'
                df_trial_cond = df_trial.query(f"isChallenge == 1")
            else:
                trial_cond = 'RB'
                df_trial_cond = df_trial.query(f"isControl == 1")
            
            trial_unpleasantness = df_trial['trialUnpleasantness'].unique()[0]
            trial_anxiety = df_trial['trialAnxiety'].unique()[0]
            
            first_i, last_i = df_trial_cond['breathIdx'].values[0], df_trial_cond['breathIdx'].values[-1]

            if first_i == last_i:
                linear_data = ecg_epoch[first_i]
            else:
                linear_data = np.concat(ecg_epoch[first_i:last_i])

            if np.isnan(linear_data).sum() != 0:
                linear_data[np.isnan(linear_data)] = np.median(linear_data)

            clean_ecg = physio.preprocess(linear_data, srate)

            if clean_ecg.size/srate < ecg_extract_time * 2:

                _df_start = pd.DataFrame({'sujet' : [sujet], 'trial_num' : [trial_i_sel], 'cond' : [trial_cond], 'time_seg' : ['start'], 'HR_count' : [0], 
                                        'ecg_length' : [clean_ecg.size/srate], 'included' : [False], 'unpleasantness' : [trial_unpleasantness], 'anxiety' : [trial_anxiety]})
                _df_end = pd.DataFrame({'sujet' : [sujet], 'trial_num' : [trial_i_sel], 'cond' : [trial_cond], 'time_seg' : ['end'], 'HR_count' : [0], 
                                        'ecg_length' : [clean_ecg.size/srate], 'included' : [False], 'unpleasantness' : [trial_unpleasantness], 'anxiety' : [trial_anxiety]})
                df_ecg.append(_df_start)
                df_ecg.append(_df_end)
                
            else:

                start_sel_vec = np.arange(clean_ecg.size)/srate <= ecg_extract_time
                end_sel_vec = np.arange(clean_ecg.size)/srate >= (clean_ecg.size/srate - ecg_extract_time) 

                extract_ecg_phase = {'start' : start_sel_vec, 'end' : end_sel_vec}

                #phase_name, phase_sel_vec = 'end', end_sel_vec
                for phase_name, phase_sel_vec in extract_ecg_phase.items():

                    ecg_phase_sel = clean_ecg[extract_ecg_phase[phase_name]]
                    raw_ecg_peak = physio.detect_peak(ecg_phase_sel, srate)
                    ecg_R_peaks = physio.clean_ecg_peak(ecg_phase_sel, srate, raw_ecg_peak)
                    df_ecg_R_peaks = pd.DataFrame({'peak_index' : ecg_R_peaks, 'peak_time' : ecg_R_peaks/srate})

                    if df_ecg_R_peaks.shape[0] < ecg_lower_limit_nHR:
                        df_add = pd.DataFrame({'sujet' : [sujet], 'trial_num' : [trial_i_sel], 'cond' : [trial_cond], 'time_seg' : [phase_name], 'HR_count' : [df_ecg_R_peaks.shape[0]], 
                                        'ecg_length' : [clean_ecg.size/srate], 'included' : [False], 'unpleasantness' : [trial_unpleasantness], 'anxiety' : [trial_anxiety]})
                        df_ecg.append(df_add)

                    else:
                        peak_ms = df_ecg_R_peaks['peak_time'].values * 1000.
                        delta_ms = np.diff(peak_ms)
                        metrics = pd.Series(dtype=float)
                        metrics['HRV_Mean'] = np.nanmean(delta_ms)
                        metrics['HRV_SD'] = np.nanstd(delta_ms)
                        metrics['HRV_Median'], metrics['HRV_Mad'] = physio.compute_median_mad(delta_ms[~np.isnan(delta_ms)])
                        metrics['HRV_CV'] = metrics['HRV_SD'] / metrics['HRV_Mean']
                        metrics['HRV_MCV'] = metrics['HRV_Mad'] / metrics['HRV_Median']
                        metrics['HRV_Asymmetry'] = metrics['HRV_Median'] - metrics['HRV_Mean']
                        metrics['HRV_RMSSD'] = np.sqrt(np.nanmean(np.diff(delta_ms)**2))

                        df = metrics.reset_index()
                        df = df.pivot_table(columns='index', values=0)

                        df_add = pd.DataFrame({'sujet' : [sujet], 'trial_num' : [trial_i_sel], 'cond' : [trial_cond], 'time_seg' : [phase_name], 'HR_count' : [df_ecg_R_peaks.shape[0]], 
                                            'ecg_length' : [clean_ecg.size/srate], 'included' : [True], 'unpleasantness' : [trial_unpleasantness], 'anxiety' : [trial_anxiety]})
                        df_add = pd.concat([df_add, df], axis=1)
                        df_ecg.append(df_add)

                    if debug:
                        plt.plot(ecg_phase_sel)
                        plt.vlines(df_ecg_R_peaks['peak_index'], ymin=ecg_phase_sel.min(), ymax=ecg_phase_sel.max(), color='r')
                        plt.show()

    df_ecg = pd.concat(df_ecg)  

    #### export  
    path_df_ecg_export = os.path.join(path_precompute, 'ECG', f"df_ecg.xlsx")
    df_ecg.to_excel(path_df_ecg_export)





################################
######## EXPORT ECG SIG ########
################################



def export_chunk_sig_ecg():

    phase_ecg_list = ['start', 'end']
    data_xr = []

    #sujet_i, sujet = 0, sujet_list_ecg[0]
    for sujet_i, sujet in enumerate(sujet_list_ecg):

        resp_epoch, ecg_epoch = get_resp_ecg_sujet(sujet)

        if debug:

            plt.plot(np.concat(ecg_epoch))
            plt.show()

            plt.plot(physio.preprocess(np.concat(ecg_epoch), srate))
            plt.show()

        os.chdir(os.path.join(path_precompute, 'RESP', 'respfeatures')) 
        # df_respfeature_label.to_excel(f'{sujet}_respfeatures_cleaned_label.xlsx')
        # df_respfeature_label_pre.to_excel(f'{sujet}_respfeatures_cleaned_label_pre.xlsx')
        df_resp_cycle_raw = pd.read_excel(f"{sujet}_cycles_info_cleaned.xlsx")

        len_extract = srate*ecg_chunk_time*len(phase_ecg_list)

        trial_i_list = df_resp_cycle_raw['trial'].unique()

        cond_alltrial_sujet = []
        data_alltrial_sujet = []

        for trial_i_add, trial_i_sel in enumerate(trial_i_list):

            # print(trial_i_sel)

            data_trial = np.zeros((len(phase_ecg_list), srate*ecg_chunk_time*len(phase_ecg_list)))

            df_trial = df_resp_cycle_raw.query(f"trial == {trial_i_sel}")

            if df_trial['isChallenge'].sum() != 0:
                trial_cond = 'CHL'
            else:
                trial_cond = 'RB'
            
            first_i, last_i = df_trial['breathIdx'].values[0], df_trial['breathIdx'].values[-1]

            if first_i == last_i:
                linear_data = ecg_epoch[first_i]
            else:
                linear_data = np.concat(ecg_epoch[first_i:last_i])

            if np.isnan(linear_data).sum() != 0:
                linear_data[np.isnan(linear_data)] = np.median(linear_data)

            clean_ecg = physio.preprocess(linear_data, srate)

            if trial_cond == 'RB':

                extract_vec_sel_strat = np.arange(clean_ecg.size)[:len_extract]
                extract_vec_sel_end = np.arange(clean_ecg.size)[clean_ecg.size - len_extract:]

            elif trial_cond == 'CHL':

                transition_cond = np.where(np.abs(np.diff(df_trial["isChallenge"])))[0] + 1

                if transition_cond.size < 2:
                    continue

                len_epoch = ecg_epoch.shape[-1]
                transition_cond_i = transition_cond * len_epoch
                extract_vec_sel_strat = np.arange(clean_ecg.size)[transition_cond_i[0] - int(len_extract/2):transition_cond_i[0] + int(len_extract/2)]
                extract_vec_sel_end = np.arange(clean_ecg.size)[transition_cond_i[1] - int(len_extract/2):transition_cond_i[1] + int(len_extract/2)]

            extract_ecg_start = linear_data[extract_vec_sel_strat]
            extract_ecg_end = linear_data[extract_vec_sel_end]

            #### construct iHR START
            extract_ecg_start_clean = physio.preprocess(extract_ecg_start, srate)
            raw_ecg_peak = physio.detect_peak(extract_ecg_start_clean, srate)
            ecg_R_peaks = physio.clean_ecg_peak(extract_ecg_start_clean, srate, raw_ecg_peak)

            if ecg_R_peaks.size < 10:
                continue

            if debug:
                plt.plot(extract_ecg_start)
                plt.vlines(ecg_R_peaks, ymin=extract_ecg_start.min(), ymax=extract_ecg_start.max(), color='r')
                plt.show()

            RRI = np.insert(np.diff(ecg_R_peaks), 0, np.median(np.diff(ecg_R_peaks)))
            RRI_start = np.zeros((extract_ecg_start.size))

            for peak_i, peak_i_val in enumerate(ecg_R_peaks):
                if peak_i == 0:
                    RRI_start[:peak_i_val] = RRI[0]
                elif peak_i == ecg_R_peaks.size-1:
                    RRI_start[ecg_R_peaks[peak_i-1]:peak_i_val] = RRI[peak_i]
                    RRI_start[peak_i_val:] = RRI[-1]
                else:
                    RRI_start[ecg_R_peaks[peak_i-1]:peak_i_val] = RRI[peak_i]

            if debug:
                plt.plot(scipy.stats.zscore(extract_ecg_start))
                plt.plot(scipy.stats.zscore(RRI_start))
                plt.vlines(ecg_R_peaks, ymin=scipy.stats.zscore(extract_ecg_start).min(), ymax=scipy.stats.zscore(extract_ecg_start).max(), color='r')
                plt.show()

            #### construct iHR END
            extract_ecg_end_clean = physio.preprocess(extract_ecg_end, srate)
            raw_ecg_peak = physio.detect_peak(extract_ecg_end_clean, srate)
            ecg_R_peaks = physio.clean_ecg_peak(extract_ecg_end_clean, srate, raw_ecg_peak)

            if ecg_R_peaks.size < 10:
                continue

            if debug:
                plt.plot(extract_ecg_end)
                plt.vlines(ecg_R_peaks, ymin=extract_ecg_end.min(), ymax=extract_ecg_end.max(), color='r')
                plt.show()

            RRI = np.insert(np.diff(ecg_R_peaks), 0, np.median(np.diff(ecg_R_peaks)))
            RRI_end = np.zeros((extract_ecg_end.size))

            for peak_i, peak_i_val in enumerate(ecg_R_peaks):
                if peak_i == 0:
                    RRI_end[:peak_i_val] = RRI[peak_i]
                elif peak_i == ecg_R_peaks.size-1:
                    RRI_end[ecg_R_peaks[peak_i-1]:peak_i_val] = RRI[peak_i]
                    RRI_end[peak_i_val:] = RRI[peak_i]
                else:
                    RRI_end[ecg_R_peaks[peak_i-1]:peak_i_val] = RRI[peak_i]

            if debug:
                plt.plot(scipy.stats.zscore(extract_ecg_end))
                plt.plot(scipy.stats.zscore(RRI_end))
                plt.vlines(ecg_R_peaks, ymin=scipy.stats.zscore(extract_ecg_end).min(), ymax=scipy.stats.zscore(extract_ecg_end).max(), color='r')
                plt.show()

                plt.plot(RRI_end)
                plt.vlines(ecg_R_peaks, ymin=scipy.stats.zscore(extract_ecg_end).min(), ymax=scipy.stats.zscore(extract_ecg_end).max(), color='r')
                plt.show()

                plt.plot(RRI)
                plt.show()

            #### load
            data_trial[0] = RRI_start
            data_trial[1] = RRI_end

            if debug:

                plt.plot(RRI_start)
                plt.plot(RRI_end)
                plt.show()

            data_alltrial_sujet.append(data_trial)
            cond_alltrial_sujet.append(trial_cond)

        data_alltrial_sujet = np.stack(data_alltrial_sujet, axis=0)
        cond_alltrial_sujet = np.array(cond_alltrial_sujet)

        RB_median = np.median(data_alltrial_sujet[cond_alltrial_sujet == 'RB'], axis=0)
        CHL_median = np.median(data_alltrial_sujet[cond_alltrial_sujet == 'CHL'], axis=0)
        data_median_trial = np.stack([RB_median, CHL_median])

        data_xr.append(data_median_trial)

    data_xr = np.stack(data_xr, axis=0)
    coords_xr = {'patient' : sujet_list_ecg, 'cond' : ['RB', 'CHL'], 'phase' : phase_ecg_list, 'time' : np.arange(-ecg_chunk_time, ecg_chunk_time, 1/srate)}

    xr_ecg = xr.DataArray(data_xr, dims=coords_xr.keys(), coords=coords_xr)

    #### export  
    path_xr_ecg_export = os.path.join(path_precompute, 'ECG', f"ecg_xr.nc")
    xr_ecg.to_netcdf(path_xr_ecg_export)








################################
######## EXECUTE ########
################################


if __name__ == '__main__':

    export_df_ecg()
    export_chunk_sig_ecg()





                        