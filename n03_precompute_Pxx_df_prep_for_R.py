

import joblib


from n00_config_params import *
from n00bis_config_analysis_functions import *
from n01_manip_data import *

import joblib
import statsmodels.api as sm
from io import StringIO
from statsmodels.tools.sm_exceptions import ValueWarning
import warnings





################################
######## EXPORT DF ########
################################

def export_df():
        
    #### full cycle
    for band in freq_band_dict:

        print(band)
        
        #### load
        os.chdir(os.path.join(path_precompute, 'TF', 'session'))

        df_band_Pxx_allsujet = []

        for sujet in sujet_list:

            df_band_Pxx_allsujet.append(xr.open_dataarray(f"{sujet}_xr_Pxx.nc").sel(pre_post='post', band=band).median('phase').to_dataframe().reset_index().dropna())

        df_band_Pxx_allsujet = pd.concat(df_band_Pxx_allsujet)
        ROI_corr = modify_loca_name(df_band_Pxx_allsujet['ROI'].values)[:,0]
        df_band_Pxx_allsujet['ROI'] = ROI_corr

        df_band_Pxx_allsujet[['resp', 'state']] = df_band_Pxx_allsujet['cond'].str.split('_', expand=True).rename(columns={0: 'resp', 1: 'state'})
        df_band_Pxx_allsujet = df_band_Pxx_allsujet.drop(columns=['cond'])

        #### export
        os.chdir(os.path.join(path_precompute, 'TF', 'session', 'df_R'))
        df_band_Pxx_allsujet.to_excel(f"df_R_{band}.xlsx")

    #### phase
    for band in freq_band_dict:

        #### load
        os.chdir(os.path.join(path_precompute, 'TF', 'session'))

        df_band_Pxx_allsujet = []

        for sujet in sujet_list:

            df_band_Pxx_allsujet.append(xr.open_dataarray(f"{sujet}_xr_Pxx.nc").sel(pre_post='post', band=band).to_dataframe().reset_index().dropna())

        df_band_Pxx_allsujet = pd.concat(df_band_Pxx_allsujet)
        ROI_corr = modify_loca_name(df_band_Pxx_allsujet['ROI'].values)[:,0]
        df_band_Pxx_allsujet['ROI'] = ROI_corr

        df_band_Pxx_allsujet[['resp', 'state']] = df_band_Pxx_allsujet['cond'].str.split('_', expand=True).rename(columns={0: 'resp', 1: 'state'})
        df_band_Pxx_allsujet = df_band_Pxx_allsujet.drop(columns=['cond'])

        #### export
        os.chdir(os.path.join(path_precompute, 'TF', 'session', 'df_R'))
        df_band_Pxx_allsujet.to_excel(f"df_R_{band}_phase.xlsx")





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

    if os.path.exists(os.path.join(path_precompute, 'TF', 'session', 'df_reg', f"df_reg_allsujet.xlsx")):
        
        print('ALREADY COMPUTED')
    
    #### load Pxx
    das = []
    for sujet in sujet_list:
        fp = os.path.join(path_precompute, "TF", "session", f"{sujet}_xr_Pxx.nc")
        da = xr.load_dataarray(fp)
        da = da.expand_dims(sujet=[sujet])

        das.append(da)

    da_Pxx_allsujet_phase_cycle = xr.concat(das, dim="sujet")
    da_Pxx_allsujet_whole_cycle = da_Pxx_allsujet_phase_cycle.median('phase')
    
    #### load resp features
    rf_metrics_tot_sel = ['cycle_duration', 'inspi_duration', 'expi_duration', 'cycle_freq', 'inspi_volume',
       'expi_volume', 'total_amplitude', 'inspi_amplitude', 'expi_amplitude',
       'total_volume']
    
    resp_vars = {}
    oc_respfeatures = []

    for sujet in sujet_list:

        # load respiration features
        os.chdir(os.path.join(path_results, 'respi', 'oc_ratio'))
        _oc_respfeatures = pd.read_excel(f'{sujet}_df_oc.xlsx')
        oc_respfeatures.append(_oc_respfeatures)

        os.chdir(os.path.join(path_precompute, 'RESP', 'respfeatures')) 
        resp_features_cleaned = pd.read_excel(f'{sujet}_respfeatures_cleaned.xlsx')
        df_resp_cycle_cleaned = pd.read_excel(f"{sujet}_cycles_info_cleaned.xlsx")

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

        _cond_da = {}
        
        for cond_i, cond in enumerate(conditions):

            resp_pre  = resp_features_cleaned.iloc[cond_sel_dict_pre[cond]].drop(columns=['Unnamed: 0'])[rf_metrics_tot_sel]
            resp_post = resp_features_cleaned.iloc[cond_sel_dict[cond]].drop(columns=['Unnamed: 0'])[rf_metrics_tot_sel]

            # convert to numpy
            pre_vals  = resp_pre.to_numpy()
            post_vals = resp_post.to_numpy()
            feat_names = list(resp_pre.columns)

            n_cycle = min(pre_vals.shape[0], post_vals.shape[0])

            ds_resp = xr.Dataset(
                data_vars={
                    "respfeatures": (("phase_group", "cycle", "feature"),
                                    np.stack([pre_vals[:n_cycle], post_vals[:n_cycle]], axis=0))
                },
                coords={
                    "phase_group": ["pre", "post"],
                    "cycle": np.arange(n_cycle),
                    "feature": feat_names,
                    "sujet": sujet,
                    "cond": cond,
                }
            ).expand_dims(("sujet", "cond"))

            # collect
            _cond_da[cond] = ds_resp

        resp_vars[sujet] = xr.concat(list(_cond_da.values()), dim="cond").sortby("cond")

    df_oc_respfeatures = pd.concat(oc_respfeatures)
    ds_resp_all = xr.concat(list(resp_vars.values()), dim="sujet").sortby("sujet")

    time_phase_list = ['pre', 'post']
    
    rf_metrics_whole_cycle = ['cycle_duration', 'cycle_freq', 'total_amplitude', 'total_volume']
    rf_metrics_phase_cycle = ['inspi_duration', 'expi_duration', 'inspi_volume', 'expi_volume', 'inspi_amplitude', 'expi_amplitude']
    rf_oc_metrics = ['oc_ratio', 'oc_val']

    ds_respfeatures_whole_cycle = ds_resp_all.sel(feature=rf_metrics_whole_cycle)
    ds_respfeatures_phase_cycle = ds_resp_all.sel(feature=rf_metrics_phase_cycle)

    ds_respfeatures_whole_cycle = ds_respfeatures_whole_cycle.rename({"phase_group": "pre_post"})
    ds_respfeatures_phase_cycle = ds_respfeatures_phase_cycle.rename({"phase_group": "pre_post"})

    ds_all_whole_cycle = xr.merge([da_Pxx_allsujet_whole_cycle, ds_respfeatures_whole_cycle])
    ds_all_phase_cycle = xr.merge([da_Pxx_allsujet_phase_cycle, ds_respfeatures_phase_cycle])

    #### compute OLS whole
    percentile_thresh_list = [2.5, 97.5]
    n_perm_spearman = 500
    time_phase_list = ['post']

    reg_allsujet_data = []

    for sujet in sujet_list:

        print(sujet)

        # chan_i, chan = 0, chanlist_sujet[0]
        def generate_reg_patient(sujet, chan_i, chan):
        # for chan_i, chan in enumerate(chanlist_sujet):

            chanlist_sujet = get_chanlist(sujet)[0]
            print_advancement(chan_i, len(chanlist_sujet), [25,50,75])

            reg_chan = []

            for cond in conditions:

                for band in freq_band_dict:

                    for phase_protocol in time_phase_list:

                        #### classic respfeatures
                        for feature in rf_metrics_whole_cycle:

                            _ds = ds_all_whole_cycle.sel(sujet=sujet, chan=chan, cond=cond, band=band, pre_post=phase_protocol, feature=feature)
                            X, Y = _ds['Pxx'].values[~np.isnan(_ds['Pxx'].values)], _ds['respfeatures'].values[~np.isnan(_ds['respfeatures'].values)]

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


                            _df = pd.DataFrame({'sujet' : [sujet], 'chan' : [chan], 'cond' : [cond], 'band' : [band], 'phase_protocol' : [phase_protocol], 'rf_metric' : [feature],
                                                'rho' : [_rho], 'rho_signi' : [_signi*1], 'OLS_a' : coeff, 'OLS_p' : pval})

                            reg_chan.append(_df)

                        #### occlusion
                        for feature in rf_oc_metrics:

                            _ds = ds_all_whole_cycle.sel(sujet=sujet, chan=chan, cond=cond, band=band, pre_post=phase_protocol)
                            X = _ds['Pxx'].values[~np.isnan(_ds['Pxx'].values)]
                            Y = df_oc_respfeatures.query(f"sujet == '{sujet}' and cond == '{cond}'")[feature].values

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


                            _df = pd.DataFrame({'sujet' : [sujet], 'chan' : [chan], 'cond' : [cond], 'band' : [band], 'phase_protocol' : [phase_protocol], 'rf_metric' : [feature],
                                                'rho' : [_rho], 'rho_signi' : [_signi*1], 'OLS_a' : coeff, 'OLS_p' : pval})

                            reg_chan.append(_df)

            return pd.concat(reg_chan)

        chanlist_sujet = get_chanlist(sujet)[0]
        reg_patient_chan = joblib.Parallel(n_jobs = n_core, prefer = 'processes')(joblib.delayed(generate_reg_patient)(sujet, chan_i, chan) for chan_i, chan in enumerate(chanlist_sujet))

        reg_allsujet_data.append(pd.concat(reg_patient_chan))

    reg_allsujet_data = pd.concat(reg_allsujet_data)

    os.chdir(os.path.join(path_precompute, 'TF', 'session', 'df_reg'))
    reg_allsujet_data.to_excel(f"df_reg_allsujet.xlsx")

    reg_allsujet_data_R = reg_allsujet_data.copy()
    reg_allsujet_data_R[['resp', 'state']] = reg_allsujet_data['cond'].str.split('_', expand=True).rename(columns={0: 'resp', 1: 'state'})
    reg_allsujet_data_R = reg_allsujet_data_R.drop(columns=['cond'])
    reg_allsujet_data_R = reg_allsujet_data_R.query(f"phase_protocol == 'post'").drop(columns=['phase_protocol'])

    os.chdir(os.path.join(path_precompute, 'TF', 'session', 'df_R'))
    reg_allsujet_data_R.to_excel(f"df_reg_allsujet_whole_cycle.xlsx")


    return reg_allsujet_data








################################
######## EXECUTE ########
################################


if __name__ == '__main__':

    export_df()
    get_reg_allsujet()


                        