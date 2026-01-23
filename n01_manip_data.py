

import h5py


from n00_config_params import *
from n00bis_config_analysis_functions import *



debug = False



"""

I ) Jose says increases in pressure during the challenge part. 
He wants to see if this increase is linked to a specific increase in gamma band activity.
The objective is to take succesive respiratory cycle during challenge that have high pressure difference and correlate it to gamma responses.
The hypothesis behid : maybe changes in gamma activities make patient aware in their breathing and make them breath with more pressure.

II ) During the occlusion there is an increase in pressure. We want to see if during the challenge these increase are different.
We take pressure during occlusion during baseline and rest, take the increase and then see if its different and if the gamma power is different. 

"""





################################
######## LOAD DATA ########
################################



#sujet = sujet_list[0]
def get_data_sujet_fullsig(sujet):

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

            chanlist = np.array(["".join(chr(c) for c in f[_ref][:].flatten()) for _ref in data['labels'][0]])
            localist = np.array(["".join(chr(c) for c in f[_ref][:].flatten()) for _ref in data['labelsFSurf'][0]])

    return data_epoch, resp_epoch, chanlist, localist




#sujet = sujet_list[0]
def get_data_sujet(sujet):

    print(f"open data {sujet}")
    os.chdir(os.path.join(path_prep))

    df_resp_cycle = pd.read_excel(f"{sujet}_breathInfo.xlsx")

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

    extract_good_cycle = (df_resp_cycle['isGoodTrial'] == 1).values & (df_resp_cycle['isGoodBreath'] == 1).values 

    try:
        data = scipy.io.loadmat(f"RRBO_{sujet}.mat")
        print('open loadmat')
        print(data.keys())
        print(data['RRBO'][0,0].dtype.names)
        data['RRBO'][0,0]['pressure']

        data_epochs_raw = np.transpose(data['RRBO'][0,0]['eeg'][:], (0,2,1)) #chan, breaths, time points 
        resp_epochs_raw = np.transpose(data['RRBO'][0,0]['pressure'][:], (1,0))
        
        data_epoch = data_epochs_raw[:,extract_good_cycle] #chan, time points, breaths
        resp_epoch = resp_epochs_raw[extract_good_cycle]
        chanlist = np.array([_chan[0][0] for _chan in data['RRBO'][0,0]['labels']])
        resp_time = data['RRBO'][0,0]['time'][:].reshape(-1)
        localist = np.array([_chan[0][0] for _chan in data['RRBO'][0,0]['labelsFSurf']])

    except:
        # print('open h5py')
        with h5py.File(f"RRBO_{sujet}.mat", "r") as f:

            print(list(f.keys()))

            data = f["RRBO"]
            print(list(data.keys()))

            data_epochs_raw = np.transpose(data['eeg'][:], (2,0,1))
            resp_epochs_raw = data['pressure'][:]

            data_epoch = data_epochs_raw[:,extract_good_cycle] #breaths, time points, chan
            resp_epoch = resp_epochs_raw[extract_good_cycle]
            chanlist = np.array(["".join(chr(c) for c in f[_ref][:].flatten()) for _ref in data['labels'][0]])
            resp_time = data['time'][:].reshape(-1)
            localist = np.array(["".join(chr(c) for c in f[_ref][:].flatten()) for _ref in data['labelsFSurf'][0]])

    os.chdir(os.path.join(path_prep))
    df_resp_cycle_cleaned = df_resp_cycle[extract_good_cycle]

    resp_control_sel = (df_resp_cycle_cleaned['occlusionType'] == 0).values & (df_resp_cycle_cleaned['isControl'] == 1).values 
    resp_chall_sel = (df_resp_cycle_cleaned['occlusionType'] == 0).values & (df_resp_cycle_cleaned['isChallenge'] == 1).values 
    oc_control_sel = (df_resp_cycle_cleaned['occlusionType'] == 2).values & (df_resp_cycle_cleaned['isControl'] == 1).values 
    oc_chall_sel = (df_resp_cycle_cleaned['occlusionType'] == 2).values & (df_resp_cycle_cleaned['isChallenge'] == 1).values 

    data_allcond = {'rsp_ctrl' : data_epoch[:,resp_control_sel], 'rsp_chl' : data_epoch[:,resp_chall_sel],
                                'oc_ctrl' : data_epoch[:,oc_control_sel], 'oc_chl' : data_epoch[:,oc_chall_sel]}
    
    resp_allcond = {'rsp_ctrl' : resp_epoch[resp_control_sel], 'rsp_chl' : resp_epoch[resp_chall_sel],
                                'oc_ctrl' : resp_epoch[oc_control_sel], 'oc_chl' : resp_epoch[oc_chall_sel]}


    return data_allcond, resp_allcond, chanlist, localist






#sujet = sujet_list[1]
def get_data_sujet_wholesession(sujet):

    print(f"open data {sujet}")
    os.chdir(os.path.join(path_prep))

    #### trial id
    df_resp_cycle = pd.read_excel(f"{sujet}_breathInfo.xlsx")

    filename_list = [_file for _file in os.listdir() if _file.find(sujet) != -1 and _file.find('RRET') != -1]
    filename_list.sort()

    data_chunked = []
    inspi_starts_i_chunked = []
    expi_starts_i_chunked = []
    respi_chunked = []

    trial_tot_count = 0
    trial_count_list = []
    trial_to_keep_totcount = []
    trial_to_keep_infile = []

    for filename_i, filename in enumerate(filename_list):

        print(f"{sujet} open trial{filename_i+1}")

        #### open
        try:
            data_load = scipy.io.loadmat(filename)
            print('open loadmat')

            n_trial = data_load["ecog"]["epochs"][0][0]['eeg'][0,0].shape[0]

            vec_nan = []

            data = []
            respi = []
            inspi_starts = []
            expi_starts = []

            for trial_i in range(n_trial):

                _sig_data = data_load["ecog"]["epochs"][0][0]['eeg'][0,0][trial_i][0]
                if np.isnan(_sig_data).sum() != 0:
                    _sel_exclude_nan = np.isnan(_sig_data).sum(axis=0) != 0
                    data.append(_sig_data)
                    vec_nan.append(_sel_exclude_nan)
                else:
                    data.append(_sig_data)

                _sig = data_load["ecog"]["epochs"][0][0]['inhonsets'][0,0][trial_i][0]
                inspi_starts.append(_sig.reshape(-1))

                _sig = data_load["ecog"]["epochs"][0][0]['inhoffsets'][0,0][trial_i][0]
                expi_starts.append(_sig.reshape(-1))

                _sig = data_load["ecog"]["epochs"][0][0]['iflow'][0,0][trial_i][0]
                iflow = _sig.reshape(-1)

                _sig = data_load["ecog"]["epochs"][0][0]['eflow'][0,0][trial_i][0]
                eflow = _sig.reshape(-1)

                _respi = iflow + eflow
                
                respi.append(_respi)

            timevec = np.array(data_load["ecog"]["epochs"][0][0]['time'][0,0]).reshape(-1)

        except:
            # print('open h5py')
            with h5py.File(filename, "r") as f:

                print('open h5py')

                # print(list(f.keys()))
                # print(list(f["ecog"].keys()))
                # print(list(f["ecog"]["epochs"]))
                # print(list(f["ecog"]["epochs"]['inhonsets'][:]))
                # print(list(f["ecog"]["epochs"]['time'][:]))
                # print(list(f["ecog"]["params"]))

                n_trial = f["ecog"]["epochs"]['eeg'][:].shape[-1]

                vec_nan = []

                data = []
                respi = []
                inspi_starts = []
                expi_starts = []

                for trial_i in range(n_trial):

                    ref = f["ecog"]["epochs"]['eeg'][0][trial_i]
                    data.append(np.transpose(f[ref][()], (1,0)))

                    ref = f["ecog"]["epochs"]['inhonsets'][0][trial_i]
                    inspi_starts.append(f[ref][()].reshape(-1))

                    ref = f["ecog"]["epochs"]['inhoffsets'][0][trial_i]
                    expi_starts.append(f[ref][()].reshape(-1))

                    ref = f["ecog"]["epochs"]['iflow'][()][0][trial_i]
                    iflow = f[ref][()].reshape(-1)

                    ref = f["ecog"]["epochs"]['eflow'][()][0][trial_i]
                    eflow = f[ref][()].reshape(-1) * -1

                    _respi = iflow + eflow
                    
                    respi.append(_respi)

                timevec = np.array(f["ecog"]["epochs"]['time'][:]).reshape(-1)

        data = np.array(data)
        respi = np.array(respi)

        if len(vec_nan) != 0:

            for _vec_nan_i, _vec_nan in enumerate(vec_nan):

                if _vec_nan_i == 0:

                    vec_nan_final = _vec_nan

                else:

                    vec_nan_final = vec_nan_final | _vec_nan

            data = data[:,:,vec_nan_final]
            respi = respi[:,vec_nan_final]

        inspi_starts_i = []
        expi_starts_i = []

        for trial_i in range(n_trial):

            inspi_starts_i.append(np.searchsorted(timevec, inspi_starts[trial_i]))
            expi_starts_i.append(np.searchsorted(timevec, expi_starts[trial_i]))

        #### verif ccle count
        for trial_i in range(n_trial):

            _n_cycle = df_resp_cycle.query(f"trial == {trial_i+1+trial_tot_count}").shape[0]

            if inspi_starts_i[trial_i].size == _n_cycle:

                trial_to_keep_totcount.append(trial_i+1+trial_tot_count)
                trial_to_keep_infile.append(trial_i)

        #### chunk
        for trial_i in range(n_trial):

            _start = inspi_starts_i[trial_i][0] - chunk_time * srate
            _stop = expi_starts_i[trial_i][-1] + chunk_time * srate
            data_chunked.append(data[trial_i,:,_start:_stop])
            respi_chunked.append(respi[trial_i][_start:_stop])

            inspi_starts_i_chunked.append(inspi_starts_i[trial_i] - _start)
            expi_starts_i_chunked.append(expi_starts_i[trial_i] - _start)

        trial_tot_count += respi.shape[0]
        trial_count_list.append(np.arange(n_trial))

        if debug:

            trial_i = 0
            for trial_i in range(respi.shape[0]):
                plt.plot(respi[trial_i])
                plt.vlines(inspi_starts_i[trial_i], ymin=respi.min(), ymax=respi.max(), colors='r')
                plt.vlines(expi_starts_i[trial_i], ymin=respi.min(), ymax=respi.max(), colors='b')
                plt.show()

            for trial_i in range(respi.shape[0]):
                plt.plot(timevec, scipy.stats.zscore(respi[trial_i]) + 5*trial_i)
            plt.show()

            trial_i = 0
            for trial_i in range(respi.shape[0]):
                _sig_respi = respi_chunked[trial_i]
                plt.plot(_sig_respi)
                plt.vlines(inspi_starts_i_chunked[trial_i], ymin=_sig_respi.min(), ymax=_sig_respi.max(), colors='r')
                plt.vlines(expi_starts_i_chunked[trial_i], ymin=_sig_respi.min(), ymax=_sig_respi.max(), colors='b')
                plt.show()

    if len(data_chunked) != df_resp_cycle['trial'].unique().size:
        print(f"!!! {sujet} trial num not the same in data and df_resp_cycle !!!")
    
    #### keep only correct
    df_resp_cycle_corr = df_resp_cycle.query(f"trial in {trial_to_keep_totcount}")
    df_resp_cycle_corr.to_excel(f"{sujet}_breathInfo_corr.xlsx")
    data_chunked = [data_chunked[i] for i in np.array(trial_to_keep_totcount)-1]
    respi_chunked = [respi_chunked[i] for i in np.array(trial_to_keep_totcount)-1]
    inspi_starts_i_chunked = [inspi_starts_i_chunked[i] for i in np.array(trial_to_keep_totcount)-1]
    expi_starts_i_chunked = [expi_starts_i_chunked[i] for i in np.array(trial_to_keep_totcount)-1]

    return data_chunked, respi_chunked, inspi_starts_i_chunked, expi_starts_i_chunked, df_resp_cycle_corr
    
    




def explore_data():

    df_info_data = pd.DataFrame()
    
    #sujet_i, sujet = 0, 'NS131_02'
    for sujet_i, sujet in enumerate(sujet_list[2:]):

        os.chdir(os.path.join(path_prep))

        df_resp_cycle = pd.read_excel(f"{sujet}_breathInfo.xlsx")

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

        extract_good_cycle = (df_resp_cycle['isGoodTrial'] == 1).values & (df_resp_cycle['isGoodBreath'] == 1).values 

        try:
            data = scipy.io.loadmat(f"RRBO_{sujet}.mat")
            print(data.keys())
            print(data['RRBO'][0,0].dtype.names)
            data['RRBO'][0,0]['pressure']

            data_epoch = data['RRBO'][0,0]['eeg'][:][:,:,extract_good_cycle] #chan, time points, breaths
            resp_epoch = data['RRBO'][0,0]['pressure'][:][:,extract_good_cycle]
            chanlist = np.array([_chan[0][0] for _chan in data['RRBO'][0,0]['labels']])
            resp_time = data['RRBO'][0,0]['time'][:].reshape(-1)
            loca_list = np.array([_chan[0][0] for _chan in data['RRBO'][0,0]['labelsFSurf']])

        except:
            with h5py.File(f"RRBO_{sujet}.mat", "r") as f:

                print(list(f.keys()))

                data = f["RRBO"]
                print(list(data.keys()))

                data_epoch = data['eeg'][:]
                resp_epoch = data['pressure'][:]

                data_epoch = data_epoch[extract_good_cycle] #chan, time points, breaths
                resp_epoch = resp_epoch[extract_good_cycle]
                chanlist = np.array(["".join(chr(c) for c in f[_ref][:].flatten()) for _ref in data['labels'][0]])
                resp_time = data['time'][:].reshape(-1)
                loca_list = np.array(["".join(chr(c) for c in f[_ref][:].flatten()) for _ref in data['labelsFSurf'][0]])

        df_resp_cycle_cleaned = df_resp_cycle[extract_good_cycle]

        resp_control_sel = (df_resp_cycle_cleaned['occlusionType'] == 0).values & (df_resp_cycle_cleaned['isControl'] == 1).values 
        resp_chall_sel = (df_resp_cycle_cleaned['occlusionType'] == 0).values & (df_resp_cycle_cleaned['isChallenge'] == 1).values 
        oc_control_sel = (df_resp_cycle_cleaned['occlusionType'] == 2).values & (df_resp_cycle_cleaned['isControl'] == 1).values 
        oc_chall_sel = (df_resp_cycle_cleaned['occlusionType'] == 2).values & (df_resp_cycle_cleaned['isChallenge'] == 1).values 

        if debug:

            print(f"resp_control : {resp_control_sel.sum()}, resp_chall : {resp_chall_sel.sum()}, oc_control : {oc_control_sel.sum()}, oc_chall : {oc_chall_sel.sum()}, ")

        df_export_sujet_cond_count = pd.DataFrame({'sujet' : [sujet], 'resp_control' : [resp_control_sel.sum()], 'resp_chall' : [resp_chall_sel.sum()], 
                                                   'oc_control' : [oc_control_sel.sum()], 'oc_chall' : [oc_chall_sel.sum()]})
        
        os.chdir(os.path.join(path_results, 'respi', 'count_cycles'))
        df_export_sujet_cond_count.to_excel(f"{sujet}_count_cycles.xlsx")

        #### organize data
        data_allcond = {'rsp_ctrl' : resp_epoch[resp_control_sel], 'rsp_chl' : resp_epoch[resp_chall_sel],
                                'oc_ctrl' : resp_epoch[oc_control_sel], 'oc_chl' : resp_epoch[oc_chall_sel]}
        
        #### plot all cycles
        os.chdir(os.path.join(path_results, 'respi', 'epochs'))

        time_vec = np.linspace(-4, 4, resp_epoch.shape[-1])
        
        data_plot = resp_epoch[resp_control_sel]
        fig_resp_control, ax = plt.subplots()
        for cycle_i in range(data_plot.shape[0]):  
            plt.plot(time_vec, data_plot[cycle_i])
        plt.title(f"{sujet} : resp_control / n : {resp_control_sel.sum()}")
        # plt.show()
        fig_resp_control.savefig(f"{sujet}_allcycles_resp_control.jpeg")

        data_plot = resp_epoch[resp_chall_sel]
        fig_resp_chall, ax = plt.subplots()
        for cycle_i in range(data_plot.shape[0]):  
            plt.plot(time_vec, data_plot[cycle_i])
        plt.title(f"{sujet} : resp_chall / n : {resp_chall_sel.sum()}")
        # plt.show()
        fig_resp_chall.savefig(f"{sujet}_allcycles_resp_chall.jpeg")

        data_plot = resp_epoch[oc_control_sel]
        fig_oc_control, ax = plt.subplots()
        for cycle_i in range(data_plot.shape[0]):  
            plt.plot(time_vec, data_plot[cycle_i])
        plt.title(f"{sujet} : oc_control / n : {oc_control_sel.sum()}")
        # plt.show()
        fig_oc_control.savefig(f"{sujet}_allcycles_oc_control.jpeg")

        data_plot = resp_epoch[oc_chall_sel]
        fig_oc_chall, ax = plt.subplots()
        for cycle_i in range(data_plot.shape[0]):  
            plt.plot(time_vec, data_plot[cycle_i])
        plt.title(f"{sujet} : oc_chall / n : {oc_chall_sel.sum()}")
        # plt.show()
        fig_oc_chall.savefig(f"{sujet}_allcycles_oc_chall.jpeg")

        plt.close('all')

        #### plot mean

        time_vec = np.linspace(-4, 4, resp_epoch.shape[-1])

        fig_summary_mean, ax = plt.subplots()
        data_plot = resp_epoch[resp_control_sel]
        plt.plot(time_vec, data_plot.mean(axis=0), label=f"resp_control ({resp_control_sel.sum()})", linestyle='--', color='r')
        data_plot = resp_epoch[resp_chall_sel]
        plt.plot(time_vec, data_plot.mean(axis=0), label=f"resp_chall ({resp_chall_sel.sum()})", linestyle='-', color='r')
        data_plot = resp_epoch[oc_control_sel]
        plt.plot(time_vec, data_plot.mean(axis=0), label=f"oc_control ({oc_control_sel.sum()})", linestyle='--', color='b')
        data_plot = resp_epoch[oc_chall_sel]
        plt.plot(time_vec, data_plot.mean(axis=0), label=f"oc_chall ({oc_chall_sel.sum()})", linestyle='-', color='b')
        plt.title(f"{sujet}")
        plt.legend()
        # plt.show()

        fig_summary_mean.savefig(f"{sujet}_mean_cycles.jpeg")

        plt.close('all')

    #### extract count cycle allsujet
    df_count_cycle_allsujet = pd.DataFrame()
    
    #sujet_i, sujet = 0, 'NS131_02'
    for sujet_i, sujet in enumerate(sujet_list): 

        os.chdir(os.path.join(path_results, 'respi', 'count_cycles'))
        _df_add = pd.read_excel(f"{sujet}_count_cycles.xlsx")
        df_count_cycle_allsujet = pd.concat([df_count_cycle_allsujet, _df_add])
        df_count_cycle_allsujet = df_count_cycle_allsujet.drop(columns=['Unnamed: 0'])

    df_count_cycle_allsujet.to_excel('ALLSUJET_count_cycle.xlsx')

    #### extract chanlist

    #sujet_i, sujet = 0, 'NS131_02'
    for sujet_i, sujet in enumerate(sujet_list): 

        print(sujet)

        data_allcond, resp_allcond, chanlist, localist = get_data_sujet(sujet)

        os.chdir(os.path.join(path_precompute, 'chanlist'))
        np.save(f"{sujet}_chanlist.npy", chanlist)
        np.save(f"{sujet}_localist.npy", localist)


    #### extract pressure for all

    diff_allcond_allsujet = {}
    pressure_allcond_allsujet = {}

    #sujet_i, sujet = 0, 'NS131_02'
    for sujet_i, sujet in enumerate(sujet_list): 

        print(sujet)

        data_allcond, resp_allcond, chanlist, localist = get_data_sujet(sujet)

        diff_allcond = {}
        pressure_allcond = {}

        for cond in conditions:

            diff_vec = []
            press_vec = []
            _resp = resp_allcond[cond]

            #cycle_i = 4
            for cycle_i in range(_resp.shape[0]):

                _respi_sig = _resp[cycle_i]
                _inv_respi_sig = _respi_sig * -1
                _peaks, _ = scipy.signal.find_peaks(_inv_respi_sig, width=srate*0.1, prominence=_inv_respi_sig.std()*0.1, distance=srate*1)

                if (_peaks < _respi_sig.size/2).sum() == 0:
                    continue

                if debug:

                    plt.plot(_respi_sig)
                    plt.show()

                    plt.plot(_respi_sig)
                    plt.scatter(_peaks, _respi_sig[_peaks], color='r')
                    plt.show()

                _peak_pre_sel_i = np.where(np.diff(_peaks > _respi_sig.size/2))[0][0]
                _peak_post_sel_i = _peak_pre_sel_i + 1
                _peak_pre_val = _respi_sig[_peaks[_peak_pre_sel_i]]
                _peak_post_val = _respi_sig[_peaks[_peak_post_sel_i]]
                _diff_val = _peak_post_val - _peak_pre_val

                diff_vec.append(_diff_val)
                press_vec.append(_peak_post_val)

            diff_allcond[cond] = np.array(diff_vec)
            pressure_allcond[cond] = np.array(press_vec)

        if debug:

            plt.plot(np.arange(pressure_allcond['rsp_ctrl'].size), pressure_allcond['rsp_ctrl'], label='rsp_ctrl')
            plt.plot(np.arange(pressure_allcond['rsp_chl'].size), pressure_allcond['rsp_chl'], label='rsp_chl')
            plt.legend()  
            plt.show()

            plt.plot(np.arange(pressure_allcond['rsp_ctrl'].size), pressure_allcond['rsp_ctrl'], label='rsp_ctrl')
            plt.plot(np.arange(pressure_allcond['oc_ctrl'].size), pressure_allcond['oc_ctrl'], label='oc_ctrl')
            plt.legend()  
            plt.show()

            plt.plot(np.arange(pressure_allcond['oc_ctrl'].size), pressure_allcond['oc_ctrl'], label='oc_ctrl')
            plt.plot(np.arange(pressure_allcond['oc_chl'].size), pressure_allcond['oc_chl'], label='oc_chl')
            plt.legend()  
            plt.show()

            plt.plot(np.arange(diff_allcond['oc_ctrl'].size), diff_allcond['oc_ctrl'], label='oc_ctrl')
            plt.plot(np.arange(diff_allcond['oc_chl'].size), diff_allcond['oc_chl'], label='oc_chl')
            plt.legend()  
            plt.show()

            plt.hist(diff_allcond['rsp_ctrl'], label='rsp_ctrl', alpha=0.5, bins=100, range=[-50,50])
            plt.hist(diff_allcond['rsp_chl'], label='rsp_chl', alpha=0.5, bins=100, range=[-50,50])

            plt.hist(diff_allcond['oc_ctrl'], label='oc_ctrl', alpha=0.5, bins=100, range=[-50,50])
            plt.hist(diff_allcond['oc_chl'], label='oc_chl', alpha=0.5, bins=100, range=[-50,50])

            plt.legend()        
            plt.show()

        diff_allcond_allsujet[sujet] = diff_allcond
        pressure_allcond_allsujet[sujet] = pressure_allcond
    

    if debug:

        for sujet_i, sujet in enumerate(sujet_list):
            plt.plot(np.arange(pressure_allcond_allsujet[sujet]['rsp_ctrl'].size), scipy.stats.zscore(pressure_allcond_allsujet[sujet]['rsp_ctrl']) + sujet_i * 3, label=sujet)
        plt.title('rsp_ctrl')
        plt.legend()  
        plt.show()

        for sujet_i, sujet in enumerate(sujet_list):
            plt.plot(np.arange(pressure_allcond_allsujet[sujet]['rsp_chl'].size), scipy.stats.zscore(pressure_allcond_allsujet[sujet]['rsp_chl']) + sujet_i * 3, label=sujet)
        plt.title('rsp_chl')
        plt.legend()  
        plt.show()

        for sujet_i, sujet in enumerate(sujet_list):
            plt.plot(np.arange(pressure_allcond_allsujet[sujet]['oc_ctrl'].size), scipy.stats.zscore(pressure_allcond_allsujet[sujet]['oc_ctrl']) + sujet_i * 3, label=sujet)
        plt.title('oc_ctrl')
        plt.legend()  
        plt.show()

        for sujet_i, sujet in enumerate(sujet_list):
            plt.plot(np.arange(pressure_allcond_allsujet[sujet]['oc_chl'].size), scipy.stats.zscore(pressure_allcond_allsujet[sujet]['oc_chl']) + sujet_i * 3, label=sujet)
        plt.title('oc_chl')
        plt.legend()  
        plt.show()

        for sujet_i, sujet in enumerate(sujet_list_allcond['oc_ctrl']):
            plt.scatter(np.arange(pressure_allcond_allsujet[sujet]['oc_ctrl'].size), np.sort(scipy.stats.zscore(pressure_allcond_allsujet[sujet]['oc_ctrl']))[::-1], label=sujet)
        plt.title('oc_ctrl')
        plt.legend()  
        plt.show()

        for sujet_i, sujet in enumerate(sujet_list_allcond['oc_chl']):
            plt.scatter(np.arange(pressure_allcond_allsujet[sujet]['oc_chl'].size), np.sort(scipy.stats.zscore(pressure_allcond_allsujet[sujet]['oc_chl']))[::-1], label=sujet)
        plt.title('oc_chl')
        plt.legend()  
        plt.show()







