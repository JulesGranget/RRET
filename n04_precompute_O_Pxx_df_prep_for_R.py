

import joblib


from n00_config_O_params import *
from n00bis_config_O_analysis_functions import *
from n00ter_X_manip_data import *

import joblib
import statsmodels.api as sm
from io import StringIO
from statsmodels.tools.sm_exceptions import ValueWarning
import warnings





################################
######## EXPORT DF ########
################################

def export_Pxx_df():
        
    for band in freq_band_dict:

        print(band)
        
        #### load
        os.chdir(os.path.join(path_precompute, 'TF', 'session'))

        df_band_Pxx_allsujet = []

        for sujet in sujet_list:

            _xr = xr.open_dataarray(f"{sujet}_xr_Pxx.nc").sel(pre_post='post', band=band).rename({'phase' : 'phase_cycle'})
            _xr_whole = _xr.median('phase_cycle').expand_dims(phase_cycle=['whole'])
            _xr_allphase = xr.concat([_xr, _xr_whole], dim='phase_cycle')

            df_band_Pxx_allsujet.append(_xr_allphase.to_dataframe().reset_index().dropna())

        df_band_Pxx_allsujet = pd.concat(df_band_Pxx_allsujet)
        ROI_corr = modify_loca_name(df_band_Pxx_allsujet['ROI'].values)[:,0]
        df_band_Pxx_allsujet['ROI'] = ROI_corr

        df_band_Pxx_allsujet[['resp', 'state']] = df_band_Pxx_allsujet['cond'].str.split('_', expand=True).rename(columns={0: 'resp', 1: 'state'})
        df_band_Pxx_allsujet = df_band_Pxx_allsujet.drop(columns=['cond'])

        #### export
        os.chdir(os.path.join(path_precompute, 'TF', 'session', 'df_R'))
        df_band_Pxx_allsujet.to_excel(f"df_R_{band}.xlsx")

    








################################
######## PERM ########
################################



def spearman_permutation_test(X, Y, n_perm_spearman, percentile_thresh_list):

    rho_obs, _ = scipy.stats.spearmanr(X, Y)

    #### ensure that there is no double
    Y_raw_i = np.arange(Y.size)
    seen_perm = set()     
    rand_perm_i = []       

    for i in range(n_perm_spearman):

        _Y_perm_i = np.random.permutation(Y_raw_i)
        _perm_key = tuple(_Y_perm_i)

        if _perm_key not in seen_perm:
            seen_perm.add(_perm_key)
            rand_perm_i.append(_Y_perm_i)

    #### compute perm
    null_dist = []
    for _perm_i in rand_perm_i:

        _rho, _ = scipy.stats.spearmanr(X, Y[_perm_i])
        null_dist.append(_rho)

    null_dist = np.array(null_dist)

    lower, upper = np.percentile(null_dist, percentile_thresh_list)
    signi = rho_obs < lower or rho_obs > upper

    return rho_obs, signi













################################
######## REG ########
################################




def get_reg_allsujet():

    if os.path.exists(os.path.join(path_precompute, 'TF', 'session', 'df_reg', f"df_reg.xlsx")):
        
        print('ALREADY COMPUTED')
    
    #### load Pxx
    das = []

    for sujet in sujet_list:
        
        fp = os.path.join(path_precompute, "TF", "session", f"{sujet}_xr_Pxx.nc")
        da = xr.load_dataarray(fp)
        da = da.expand_dims(sujet=[sujet])

        das.append(da)

    da_Pxx_allsujet_phase_cycle = xr.concat(das, dim="sujet")
    da_Pxx_allsujet_phase_cycle = da_Pxx_allsujet_phase_cycle.rename({'phase' : 'phase_cycle'})
    da_whole = da_Pxx_allsujet_phase_cycle.median('phase_cycle').expand_dims(phase_cycle=["whole"])
    da_Pxx_allphase_cycle = xr.concat([da_Pxx_allsujet_phase_cycle, da_whole],dim="phase_cycle")
    
    #### load resp features
    rf_metrics_raw = ['cycle_duration', 'inspi_duration', 'expi_duration', 'cycle_freq', 'inspi_volume',
       'expi_volume', 'total_amplitude', 'inspi_amplitude', 'expi_amplitude',
       'total_volume']
    phase_cycle_list = ['whole', 'inspi', 'expi']
    mapping_pase_cycle_list = {'whole' : ['cycle_duration', 'cycle_freq', 'total_amplitude', 'total_volume', 'oc_ratio', 'oc_val'], 
                               'inspi' : ['inspi_duration', 'inspi_cycle_freq', 'inspi_amplitude', 'inspi_volume', 'oc_ratio', 'oc_val'], 
                               'expi' : ['expi_duration', 'expi_cycle_freq', 'expi_amplitude', 'expi_volume', 'oc_ratio', 'oc_val']}
    rf_metric_allphase_cycle = ['duration', 'cycle_freq', 'amplitude', 'volume', 'oc_ratio', 'oc_val']
    
    resp_vars = {}

    for sujet in sujet_list:

        # load respiration features
        os.chdir(os.path.join(path_precompute, 'RESP', 'respfeatures')) 
        _oc_respfeatures = pd.read_excel(f'{sujet}_df_oc.xlsx')

        _respfeatures = pd.read_excel(f'{sujet}_respfeatures_cleaned_label.xlsx')
        _respfeatures_pre = pd.read_excel(f'{sujet}_respfeatures_cleaned_label_pre.xlsx')

        _cond_da = {}
        
        for cond_i, cond in enumerate(conditions):

            resp_pre  = _respfeatures_pre.query(f"cond == '{cond}'").drop(columns=['Unnamed: 0'])[rf_metrics_raw]
            resp_post = _respfeatures.query(f"cond == '{cond}'").drop(columns=['Unnamed: 0'])[rf_metrics_raw]

            resp_pre[['oc_ratio', 'oc_val']] = _oc_respfeatures.query(f"cond == '{cond}'")[['oc_ratio', 'oc_val']]
            resp_post[['oc_ratio', 'oc_val']] = _oc_respfeatures.query(f"cond == '{cond}'")[['oc_ratio', 'oc_val']]

            resp_pre['inspi_cycle_freq'], resp_pre['expi_cycle_freq'] = 1/resp_pre['inspi_duration'], 1/resp_pre['expi_duration']
            resp_post['inspi_cycle_freq'], resp_post['expi_cycle_freq'] = 1/resp_post['inspi_duration'], 1/resp_post['expi_duration']
            
            allmetric_resp_pre = []
            allmetric_resp_post = []

            for _phase_cycle in phase_cycle_list:

                allmetric_resp_pre.append(resp_pre[mapping_pase_cycle_list[_phase_cycle]].to_numpy())
                allmetric_resp_post.append(resp_post[mapping_pase_cycle_list[_phase_cycle]].to_numpy())

            # convert to numpy
            pre_vals  = np.stack(allmetric_resp_pre)
            post_vals = np.stack(allmetric_resp_post)

            ds_resp = xr.Dataset(
                data_vars={
                    "respfeatures": (("phase_group", "phase_cycle", "cycle", "feature"),
                                    np.stack([pre_vals, post_vals], axis=0))
                },
                coords={
                    "phase_group": ["pre", "post"],
                    "phase_cycle": phase_cycle_list,
                    "cycle": np.arange(pre_vals.shape[1]),
                    "feature": rf_metric_allphase_cycle,
                    "sujet": sujet,
                    "cond": cond,
                }
            ).expand_dims(("sujet", "cond"))

            # collect
            _cond_da[cond] = ds_resp

        resp_vars[sujet] = xr.concat(list(_cond_da.values()), dim="cond").sortby("cond")

    ds_respfeatures_all = xr.concat(list(resp_vars.values()), dim="sujet").sortby("sujet")
    ds_respfeatures_all = ds_respfeatures_all.rename({"phase_group": "pre_post"})

    #### merge
    ds_whole = xr.merge([da_Pxx_allphase_cycle, ds_respfeatures_all])

    #### inspect
    if debug:

        time_phase_list = ['pre', 'post']

        for sujet in sujet_list:

            chanlist_sujet = get_chanlist(sujet)[0]

            for chan_i, chan in enumerate(chanlist_sujet):

                for cond in conditions:

                    for band in freq_band_dict:

                        for phase_protocol in time_phase_list:

                            for cycle_phase in phase_cycle_list:

                                for feature in rf_metric_allphase_cycle:

                                    _ds_Pxx = ds_whole.sel(sujet=sujet, chan=chan, cond=cond, band=band, pre_post='pre', phase_cycle=cycle_phase, feature=feature)['Pxx']
                                    X_pre = _ds_Pxx.values[~np.isnan(_ds_Pxx.values)]

                                    _ds_Pxx = ds_whole.sel(sujet=sujet, chan=chan, cond=cond, band=band, pre_post='post', phase_cycle=cycle_phase, feature=feature)['Pxx']
                                    X_post = _ds_Pxx.values[~np.isnan(_ds_Pxx.values)]

                                    _ds_Pxx = ds_whole.sel(sujet=sujet, chan=chan, cond=cond, band=band, pre_post='pre', phase_cycle=cycle_phase, feature='oc_ratio')['Pxx']
                                    X_pre = _ds_Pxx.values[~np.isnan(_ds_Pxx.values)]

                                    _ds_Pxx = ds_whole.sel(sujet=sujet, chan=chan, cond=cond, band=band, pre_post='post', phase_cycle=cycle_phase, feature='oc_ratio')['Pxx']
                                    X_post = _ds_Pxx.values[~np.isnan(_ds_Pxx.values)]

                                    _df_rf = ds_whole.sel(sujet=sujet, chan=chan, cond=cond, band=band, pre_post='pre', phase_cycle=cycle_phase, feature='oc_ratio')['respfeatures']
                                    Y_pre = _df_rf.values[~np.isnan(_df_rf.values)]

                                    _df_rf = ds_whole.sel(sujet=sujet, chan=chan, cond=cond, band=band, pre_post='post', phase_cycle=cycle_phase, feature='oc_ratio')['respfeatures']
                                    Y_post = _df_rf.values[~np.isnan(_df_rf.values)]

                                    _df_rf = ds_whole.sel(sujet=sujet, chan=chan, cond=cond, band=band, pre_post='pre', phase_cycle=cycle_phase, feature=feature)['respfeatures']
                                    Y_pre = _df_rf.values[~np.isnan(_df_rf.values)]

                                    _df_rf = ds_whole.sel(sujet=sujet, chan=chan, cond=cond, band=band, pre_post='post', phase_cycle=cycle_phase, feature=feature)['respfeatures']
                                    Y_post = _df_rf.values[~np.isnan(_df_rf.values)]

    #### compute OLS
    percentile_thresh_list = [2.5, 97.5]
    n_perm_spearman = 500
    time_phase_list = ['pre', 'post']

    reg_allsujet_data = []

    for sujet in sujet_list:

        print(sujet)

        if os.path.exists(os.path.join(path_precompute, 'TF', 'session', 'df_reg', f"df_reg_{sujet}.xlsx")):

            print(f"{sujet} ALREADY COMPUTED")
            os.chdir(os.path.join(path_precompute, 'TF', 'session', 'df_reg'))
            df_reg_patient = pd.read_excel(f"df_reg_{sujet}.xlsx")
        
        else:

            # chan_i, chan = 0, chanlist_sujet[0]
            def generate_reg_patient(sujet, chan_i, chan):
            # for chan_i, chan in enumerate(chanlist_sujet):

                chanlist_sujet = get_chanlist(sujet)[0]
                print_advancement(chan_i, len(chanlist_sujet), [25,50,75])

                reg_chan = []

                for cond in conditions:

                    for band in freq_band_dict:

                        for phase_protocol in time_phase_list:

                            for cycle_phase in phase_cycle_list:

                                for feature in rf_metric_allphase_cycle:

                                    _ds_Pxx = ds_whole.sel(sujet=sujet, chan=chan, cond=cond, band=band, pre_post='post', phase_cycle=cycle_phase, feature=feature)['Pxx']
                                    X = _ds_Pxx.values[~np.isnan(_ds_Pxx.values)]

                                    _df_rf = ds_whole.sel(sujet=sujet, chan=chan, cond=cond, band=band, pre_post=phase_protocol, phase_cycle=cycle_phase, feature=feature)['respfeatures']
                                    Y = _df_rf.values[~np.isnan(_df_rf.values)]

                                    _rho, _signi = spearman_permutation_test(X, Y, n_perm_spearman, percentile_thresh_list)

                                    X = sm.add_constant(X)

                                    _mdl = sm.OLS(Y, X).fit()

                                    with warnings.catch_warnings():
                                        warnings.simplefilter("ignore", category=UserWarning)
                                        warnings.simplefilter("ignore", category=ValueWarning)

                                        table = pd.read_html(
                                            StringIO(_mdl.summary().tables[1].as_html()),
                                            header=0,
                                            index_col=0
                                        )[0]

                                    coeff     = table["coef"].values[-1]
                                    pval      = table["P>|t|"].values[-1]

                                    _ROI = modify_loca_name([_ds_Pxx['ROI'].data.tolist()])[0][0]

                                    _df = pd.DataFrame({'sujet' : [sujet], 'chan' : [chan], 'ROI' : [_ROI], 'cond' : [cond], 'band' : [band], 'phase_protocol' : [phase_protocol], 'rf_metric' : [feature],
                                                        'phase_cycle' : [cycle_phase], 'rho' : [_rho], 'rho_signi' : [_signi*1], 'OLS_a' : coeff, 'OLS_p' : pval})

                                    reg_chan.append(_df)

                return pd.concat(reg_chan)

            chanlist_sujet = get_chanlist(sujet)[0]
            reg_patient_chan = joblib.Parallel(n_jobs = n_core, prefer = 'processes')(joblib.delayed(generate_reg_patient)(sujet, chan_i, chan) for chan_i, chan in enumerate(chanlist_sujet))

            df_reg_patient = pd.concat(reg_patient_chan)

            os.chdir(os.path.join(path_precompute, 'TF', 'session', 'df_reg'))
            df_reg_patient.to_excel(f"df_reg_{sujet}.xlsx")

        #### inspect
        if debug:

            _df = generate_reg_patient(sujet, chan_i, chan)
            _df.query(f"cond=='{cond}' and band=='{band}' and phase_protocol=='pre' and phase_cycle=='{cycle_phase}' and rf_metric=='oc_ratio'")
            _df.query(f"cond=='{cond}' and band=='{band}' and phase_protocol=='post' and phase_cycle=='{cycle_phase}' and rf_metric=='oc_ratio'")

            _df.query(f"cond=='{cond}' and band=='{band}' and phase_protocol=='pre' and phase_cycle=='{cycle_phase}' and rf_metric=='duration'")
            _df.query(f"cond=='{cond}' and band=='{band}' and phase_protocol=='post' and phase_cycle=='{cycle_phase}' and rf_metric=='duration'")


        #### add df
        reg_allsujet_data.append(df_reg_patient)

    reg_allsujet_data = pd.concat(reg_allsujet_data)

    os.chdir(os.path.join(path_precompute, 'TF', 'session', 'df_reg'))
    reg_allsujet_data.to_excel(f"df_reg.xlsx")

    reg_allsujet_data_R = reg_allsujet_data.copy()
    reg_allsujet_data_R[['resp', 'state']] = reg_allsujet_data['cond'].str.split('_', expand=True).rename(columns={0: 'resp', 1: 'state'})
    reg_allsujet_data_R = reg_allsujet_data_R.drop(columns=['cond'])

    os.chdir(os.path.join(path_precompute, 'TF', 'session', 'df_R'))
    reg_allsujet_data_R.to_excel(f"df_reg_R.xlsx")









################################
######## EXECUTE ########
################################


if __name__ == '__main__':

    export_Pxx_df()
    get_reg_allsujet()


                        