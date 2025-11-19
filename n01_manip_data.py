

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




conditions = ['rsp_ctrl', 'rsp_chl', 'oc_ctrl', 'oc_chl']


################################
######## LOAD DATA ########
################################


def get_data_sujet(sujet):

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
            localist = np.array(["".join(chr(c) for c in f[_ref][:].flatten()) for _ref in data['labelsFSurf'][0]])

    return data_epoch, resp_epoch, chanlist, localist







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
        np.concatenate((resp_epoch[resp_control_sel].reshape(1,resp_control_sel.sum(),-1), resp_epoch[resp_chall_sel].reshape(1,resp_chall_sel.sum(),-1)), axis=0)


        
        #### plot all cycles
        os.chdir(os.path.join(path_results, 'respi'))
        
        data_plot = resp_epoch[resp_control_sel]
        fig_resp_control, ax = plt.subplots()
        for cycle_i in range(data_plot.shape[0]):  
            plt.plot(data_plot[cycle_i])
        plt.title(f"{sujet} : resp_control / n : {resp_control_sel.sum()}")
        # plt.show()
        fig_resp_control.savefig(f"{sujet}_allcycles_resp_control.jpeg")

        data_plot = resp_epoch[resp_chall_sel]
        fig_resp_chall, ax = plt.subplots()
        for cycle_i in range(data_plot.shape[0]):  
            plt.plot(data_plot[cycle_i])
        plt.title(f"{sujet} : resp_chall / n : {resp_chall_sel.sum()}")
        # plt.show()
        fig_resp_chall.savefig(f"{sujet}_allcycles_resp_chall.jpeg")

        data_plot = resp_epoch[oc_control_sel]
        fig_oc_control, ax = plt.subplots()
        for cycle_i in range(data_plot.shape[0]):  
            plt.plot(data_plot[cycle_i])
        plt.title(f"{sujet} : oc_control / n : {oc_control_sel.sum()}")
        # plt.show()
        fig_oc_control.savefig(f"{sujet}_allcycles_oc_control.jpeg")

        data_plot = resp_epoch[oc_chall_sel]
        fig_oc_chall, ax = plt.subplots()
        for cycle_i in range(data_plot.shape[0]):  
            plt.plot(data_plot[cycle_i])
        plt.title(f"{sujet} : oc_chall / n : {oc_chall_sel.sum()}")
        # plt.show()
        fig_oc_chall.savefig(f"{sujet}_allcycles_oc_chall.jpeg")

        plt.close('all')

        #### plot mean

        fig_summary_mean, ax = plt.subplots()
        data_plot = resp_epoch[resp_control_sel]
        plt.plot(data_plot.mean(axis=0), label=f"resp_control ({resp_control_sel.sum()})", linestyle='--', color='r')
        data_plot = resp_epoch[resp_chall_sel]
        plt.plot(data_plot.mean(axis=0), label=f"resp_chall ({resp_chall_sel.sum()})", linestyle='-', color='r')
        data_plot = resp_epoch[oc_control_sel]
        plt.plot(data_plot.mean(axis=0), label=f"oc_control ({oc_control_sel.sum()})", linestyle='--', color='b')
        data_plot = resp_epoch[oc_chall_sel]
        plt.plot(data_plot.mean(axis=0), label=f"oc_chall ({oc_chall_sel.sum()})", linestyle='-', color='b')
        plt.title(f"{sujet}")
        plt.legend()
        # plt.show()

        fig_summary_mean.savefig(f"{sujet}_mean_cycles.jpeg")

        plt.close('all')

    #### extract count cycle allsujet
    df_count_cycle_allsujet = pd.DataFrame()
    
    #sujet_i, sujet = 0, 'NS131_02'
    for sujet_i, sujet in enumerate(sujet_list[2:]): 

        os.chdir(os.path.join(path_results, 'respi', 'count_cycles'))
        _df_add = pd.read_excel(f"{sujet}_count_cycles.xlsx")
        df_count_cycle_allsujet = pd.concat([df_count_cycle_allsujet, _df_add])
        df_count_cycle_allsujet = df_count_cycle_allsujet.drop(columns=['Unnamed: 0'])

    df_count_cycle_allsujet.to_excel('ALLSUJET_count_cycle.xlsx')


    #### extract pressure for all
    
    diff_allcond = {}

    for cond in conditions:

        diff_vec = []

        for cycle_i in range(resp_epoch.shape[0]):

            _respi_sig = resp_epoch[cycle_i]
            _inv_respi_sig = _respi_sig * -1
            _peaks, _ = scipy.signal.find_peaks(_inv_respi_sig, width=srate*0.1, prominence=_inv_respi_sig.std()*0.1, distance=srate*1)

            if debug:

                plt.plot(resp_epoch[cycle_i])
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

        diff_allcond[cond] = np.array(diff_vec)

    if debug:

        plt.plot(diff_vec)
        plt.show()
            











