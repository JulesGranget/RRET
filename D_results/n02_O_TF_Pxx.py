



from A_config.n01_O_params import *
from A_config.n02_O_analysis_functions import *
from A_config.n03_X_manip_data import *

from nilearn import plotting
import joblib
import statsmodels.api as sm
from io import StringIO
from statsmodels.tools.sm_exceptions import ValueWarning
import warnings





########################################
######## ANAT ALLPATIENT ########
########################################


def export_allpatient_anat():

    os.chdir(os.path.join(path_results, 'anatomy'))

    if not os.path.exists("anatomy_allpatient.xlsx"):

        anat_allpatient = []

        for sujet in sujet_list:

            chanlist, localist = get_chanlist(sujet)
            localist_corr = modify_loca_name(localist)
            _df = pd.DataFrame([{'sujet' : sujet, 'chan_name' : _chan_name, 'loca' : localist_corr[i][0], 'type' : localist_corr[i][1], 'side' : localist_corr[i][2]} for i, _chan_name in enumerate(chanlist)])
            anat_allpatient.append(_df)

        anat_allpatient = pd.concat(anat_allpatient)
        
        all_df.to_excel("anatomy_allpatient.xlsx")







################################
######## TF ABSOLUTE ########
################################


def tf_absolute_patientwise():

    #### load
    sujet = sujet_list[0]

    os.chdir(os.path.join(path_precompute, 'TF', 'session'))
    tf_allcond_dict = {}
    time_phase_list = ['pre', 'post']
    
    for time_phase in time_phase_list:
        
        tf_allcond_dict[time_phase] = {}

        for cond in conditions:

            tf_allcond_dict[time_phase][cond] = np.load(f'{sujet}_{cond}_tf_allchan_stretch_{time_phase}.npy')   

    chanlist, localist_raw = get_chanlist(sujet)
    localist = modify_loca_name(localist_raw)
    localist_unique = np.unique([_loca[0] for _loca in localist if _loca[0].find('UNSORTED') == -1])
    
    #### plot
    os.chdir(os.path.join(path_results, 'TF', 'patient_wise'))
    time_vec = np.arange(stretch_point_TF)
    ticks_freq = [2, 8, 12, 30, 50, 60, 150]
    time_phase_list = ['pre', 'post']
    percentile_plot = [1,99]

    for loca_sel in localist_unique:

        #### vlim
        _tf_loca_dict = {}
        _tf_loca_vlim = []
        _chan_sel_vec = localist[:,0] == loca_sel
        
        for time_phase_sel_i, time_phase_sel in enumerate(time_phase_list):

            _tf_loca_dict[time_phase_sel] = {}
        
            for cond_i, cond in enumerate(conditions):

                _tf = np.median(np.median(tf_allcond_dict[time_phase_sel][cond][_chan_sel_vec], axis=1), axis=0) #median cycle then chan

                if debug:

                    plt.pcolormesh(_tf)
                    plt.show()

                _tf_loca_vlim.append(_tf)
                _tf_loca_dict[time_phase_sel][cond] = _tf

        _vlim = np.abs([np.percentile(np.stack(_tf_loca_vlim, axis=0), percentile_plot[0]), np.percentile(np.stack(_tf_loca_vlim, axis=0), percentile_plot[1])]).max()

        fig, axs = plt.subplots(ncols=4, nrows=2, figsize=(16,8))

        for time_phase_sel_i, time_phase_sel in enumerate(time_phase_list):
        
            for cond_i, cond in enumerate(conditions):
            
                ax = axs[time_phase_sel_i, cond_i]
                pcm = ax.pcolormesh(time_vec, frex, _tf_loca_dict[time_phase_sel][cond], vmin=-_vlim, vmax=_vlim, cmap='seismic')
                ax.vlines(x=int(stretch_point_TF/2), ymin=frex[0], ymax=frex[-1], color="k", linewidth=2)
                ax.set_yscale('log')
                
                ax.set_yticks(ticks_freq)
                ax.set_yticklabels([str(t) for t in ticks_freq])

                if cond_i == 0:
                    ax.set_ylabel(time_phase_sel)
                if time_phase_sel_i == 1:
                    ax.set_xlabel('Phase')
                if time_phase_sel_i == 0:
                    ax.set_title(f"{cond}")

        cbar = fig.colorbar(pcm, ax=axs, orientation="vertical", fraction=0.02, pad=0.04)
        plt.suptitle(f"{sujet} {loca_sel} ({_chan_sel_vec.sum()})")

        # plt.show()
        fig.savefig(f"{sujet}_{loca_sel}.jpg")

    plt.close('all')





def tf_absolute_allpatient():

    #### load
    os.chdir(os.path.join(path_results, 'anatomy'))
    df_loca_allpatient = pd.read_excel('anatomy_allpatient.xlsx')
    localist_unique = np.unique([_loca for _loca in df_loca_allpatient['loca'].values if _loca.find('UNSORTED') == -1])

    time_phase_list = ['pre', 'post']
    ticks_freq = [2, 8, 12, 30, 50, 60, 150]
    percentile_plot = [1,99]
    time_vec = np.arange(stretch_point_TF)

    #loca_sel_i, loca_sel = 2, localist_unique[2]
    # for loca_sel_i, loca_sel in enumerate(localist_unique):
    def plot_allsujet_loca(loca_sel_i, loca_sel):

        print_advancement(loca_sel_i, len(localist_unique), [25,50,75])

        #### load sujet
        os.chdir(os.path.join(path_precompute, 'TF', 'session'))
        loca_sel_vec = df_loca_allpatient['loca'].values == loca_sel
        _patient_list_unique = np.unique(df_loca_allpatient['sujet'].values[loca_sel_vec])
        _tf_allsujet = np.zeros((_patient_list_unique.size, len(time_phase_list), len(conditions), frex.size, stretch_point_TF))
        _ncycle_tot = np.zeros((len(time_phase_list), len(conditions)))

        for _sujet_i, _sujet in enumerate(_patient_list_unique):

            _sujet_chanlist, _sujet_localist = get_chanlist(_sujet)
            _sujet_localist = modify_loca_name(_sujet_localist)
            _sujet_chan_sel_vec = _sujet_localist[:,0] == loca_sel
            
            for time_phase_i, time_phase in enumerate(time_phase_list):
                
                for cond_i, cond in enumerate(conditions):

                    _tf = np.load(f'{_sujet}_{cond}_tf_allchan_stretch_{time_phase}.npy')[_sujet_chan_sel_vec]
                    _tf_allsujet[_sujet_i, time_phase_i, cond_i, :, :] = np.median(np.median(_tf, axis=1), axis=0)
                    _ncycle_tot[time_phase_i, cond_i] += _tf.shape[1]

        _tf_allsujet_median = np.median(_tf_allsujet, axis=0)

        #### plot
        _vlim = np.abs([np.percentile(np.stack(_tf_allsujet_median, axis=0), percentile_plot[0]), np.percentile(np.stack(_tf_allsujet_median, axis=0), percentile_plot[1])]).max()

        fig, axs = plt.subplots(ncols=4, nrows=2, figsize=(16,8))

        for time_phase_sel_i, time_phase_sel in enumerate(time_phase_list):
        
            for cond_i, cond in enumerate(conditions):
            
                ax = axs[time_phase_sel_i, cond_i]
                pcm = ax.pcolormesh(time_vec, frex, _tf_allsujet_median[time_phase_sel_i, cond_i], vmin=-_vlim, vmax=_vlim, cmap='seismic')
                ax.vlines(x=int(stretch_point_TF/2), ymin=frex[0], ymax=frex[-1], color="k", linewidth=2)
                ax.set_yscale('log')
                
                ax.set_yticks(ticks_freq)
                ax.set_yticklabels([str(t) for t in ticks_freq])

                if cond_i == 0:
                    ax.set_ylabel(time_phase_sel)
                if time_phase_sel_i == 1:
                    ax.set_xlabel('Phase')
                if time_phase_sel_i == 0:
                    ax.set_title(f"{cond} c{int(_ncycle_tot[time_phase_i, cond_i])}")

        cbar = fig.colorbar(pcm, ax=axs, orientation="vertical", fraction=0.02, pad=0.04)
        plt.suptitle(f"{loca_sel} s({len(_patient_list_unique)})")

        # plt.show()

        #### save
        os.chdir(os.path.join(path_results, 'TF', 'allsujet'))
        fig.savefig(f"allsujet_{loca_sel}.jpg")

    joblib.Parallel(n_jobs = n_core, prefer = 'processes')(joblib.delayed(plot_allsujet_loca)(loca_sel_i, loca_sel) for loca_sel_i, loca_sel in enumerate(localist_unique))
    
    plt.close('all')





################################
######## PXX ABSOLUTE ########
################################

def Pxx_absolute_patientwise():
    
    sujet = sujet_list[0]

    os.chdir(os.path.join(path_precompute, 'TF', 'session'))
    xr_Pxx = xr.load_dataarray(f'{sujet}_xr_Pxx.nc')    

    df_coords_allsujet = []
    for sujet in sujet_list:
        _df_coords = get_coords(sujet)
        _df_coords['sujet'] = [sujet] * _df_coords.shape[0]
        df_coords_allsujet.append(_df_coords)
    df_coords_allsujet = pd.concat(df_coords_allsujet)

    #### load Pxx
    das = []
    for sujet in sujet_list:
        fp = os.path.join(path_precompute, "TF", "session", f"{sujet}_xr_Pxx.nc")
        da = xr.load_dataarray(fp)
        da = da.expand_dims(sujet=[sujet])

        das.append(da)

    da_Pxx_allsujet_phase_cycle = xr.concat(das, dim="sujet")
    da_Pxx_allsujet_whole_cycle = da_Pxx_allsujet_phase_cycle.median('phase')

    #### load respfeatures
    time_phase_list = ['pre', 'post']

    resp_vars = {} 
    rf_metrics_tot_sel = ['cycle_duration', 'inspi_duration', 'expi_duration', 'cycle_freq', 'inspi_volume',
       'expi_volume', 'total_amplitude', 'inspi_amplitude', 'expi_amplitude',
       'total_volume']
    rf_metrics_whole_cycle = ['cycle_duration', 'cycle_freq', 'total_amplitude', 'total_volume']
    rf_metrics_phase_cycle = ['inspi_duration', 'expi_duration', 'inspi_volume', 'expi_volume', 'inspi_amplitude', 'expi_amplitude']


    for sujet in sujet_list:

        # load respiration features
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

    ds_resp_all = xr.concat(list(resp_vars.values()), dim="sujet").sortby("sujet")

    ds_respfeatures_whole_cycle = ds_resp_all.sel(feature=rf_metrics_whole_cycle)
    ds_respfeatures_phase_cycle = ds_resp_all.sel(feature=rf_metrics_phase_cycle)

    da_Pxx_allsujet_whole_cycle = da_Pxx_allsujet_whole_cycle.to_dataset(name="Pxx")

    ds_respfeatures_whole_cycle = ds_respfeatures_whole_cycle.rename({"phase_group": "pre_post"})
    ds_respfeatures_whole_cycle = ds_respfeatures_whole_cycle.reindex(cycle=ds_Pxx.coords["cycle"].values)

    ds_all_whole_cycle = xr.merge([da_Pxx_allsujet_whole_cycle, ds_respfeatures_whole_cycle])

    #### plot
    sujet = sujet_list[0]
    chanlist_sujet = get_chanlist(sujet)[0]
    chan = chanlist_sujet[0]

    ds_plot = ds_all_whole_cycle.sel(sujet=sujet, chan=chan)
    df_plot = ds_plot.to_dataframe().reset_index().dropna()
    df_plot = df_plot.pivot_table(index=[_col for _col in df_plot.columns if _col not in ['feature', 'respfeatures']], columns="feature", values="respfeatures").reset_index()

    sns.lmplot(data=df_plot, x='total_amplitude', y='Pxx', hue='cond', col='band', row='pre_post', sharex=False, sharey=False)
    plt.show()


    #### REG STATS WHOLE CYCLE
    percentile_thresh_list = [2.5, 97.5]
    n_perm_spearman = 500

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
        
    reg_allsujet_data = []

    for sujet in sujet_list:

        print(sujet)

        def generate_reg_patient(sujet, chan_i, chan):
        # for chan_i, chan in enumerate(chanlist_sujet):

            chanlist_sujet = get_chanlist(sujet)[0]
            print_advancement(chan_i, len(chanlist_sujet), [25,50,75])

            reg_chan = []

            for cond in conditions:

                for band in freq_band_dict:

                    for phase_protocol in time_phase_list:

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

            return pd.concat(reg_chan)

        chanlist_sujet = get_chanlist(sujet)[0]
        reg_patient_chan = joblib.Parallel(n_jobs = n_core, prefer = 'processes')(joblib.delayed(generate_reg_patient)(sujet, chan_i, chan) for chan_i, chan in enumerate(chanlist_sujet))

        reg_allsujet_data.append(pd.concat(reg_patient_chan))

    reg_allsujet_data = pd.concat(reg_allsujet_data)

    reg_allsujet_data = reg_allsujet_data.merge(
        df_coords_allsujet,
        on=["sujet", "chan"],
        how="left"
        )




    #### plot sujet wise all contacts
    metric = 'total_amplitude'
    band = 'gamma'
    sujet = sujet_list[0]

    cmap = plt.cm.seismic


    for metric in rf_metrics:

        for band in freq_band_dict:

            global_max = reg_allsujet_data.query(f"sujet == '{sujet}' and rf_metric == '{metric}' and band == '{band}' and coords_info_present == True")['OLS_a'].abs().max()
            # size_max = 20.0 
            size_markers = 10
            size_markers_signi_coeff = 2

            norm = matplotlib.colors.TwoSlopeNorm(vmin=-global_max, vcenter=0.0, vmax=global_max)

            # from nilearn import plotting

            fig, axes = plt.subplots(
                nrows=len(time_phase_list), ncols=len(conditions),
                figsize=(10*len(conditions), 5*len(time_phase_list)),
                constrained_layout=True
            )

            for r, phase in enumerate(time_phase_list):
                for c, cond in enumerate(conditions):

                    ax = axes[r, c]

                    _df = reg_allsujet_data.query(f"sujet == '{sujet}' and rf_metric == '{metric}' and band == '{band}' and cond == '{cond}' and phase_protocol == '{phase}' and coords_info_present == True")
                    coords_all = _df[['coords_x', 'coords_y', 'coords_z']].values

                    # colors = np.where(_df["rho_signi"].values == 1, "crimson", "royalblue")
                    colors_all = cmap(norm(_df["OLS_a"].values))

                    signi_circle = (_df["rho_signi"] == 1).to_numpy()

                    coords_signi = coords_all[signi_circle]
                    coords_no_signi = coords_all[~signi_circle]

                    colors_signi = colors_all[signi_circle]
                    colors_no_signi = colors_all[~signi_circle]

                    # sizes = _df["OLS_a"].abs().to_numpy()
                    # sizes_scaled = (sizes / global_max) * size_max
                    # sizes_scaled = np.clip(sizes_scaled, 2, None)

                    disp = nilearn.plotting.plot_glass_brain(
                        None,
                        display_mode="xz",
                        alpha=0.1,
                        axes=ax
                    )
                    # disp.add_markers(_df[['coords_x', 'coords_y', 'coords_z']].values, marker_color=list(colors), marker_size=list(sizes_scaled))

                    disp.add_markers(
                        coords_no_signi,
                        marker_color=list(colors_no_signi),
                        marker_size=size_markers
                    )

                    # circled markers: black outer + colored inner
                    disp.add_markers(
                        coords_signi,
                        marker_color="green",
                        marker_size=size_markers*size_markers_signi_coeff  # outline thickness
                    )
                    disp.add_markers(
                        coords_signi,
                        marker_color=list(colors_signi),
                        marker_size=size_markers
                    )

                    ax.set_title(f"{cond} : {phase}")

            plt.suptitle(f"{band}\n{metric}")
            # plt.tight_layout()
            plt.show()

            os.chdir(os.path.join(path_results, 'respi', 'sujet_wise_OLS'))
            fig.savefig(f"{sujet}_{band}_{metric}_OLS.png")


    #### plot sujet wise significant contacts only
    metric = 'total_amplitude'
    band = 'gamma'
    sujet = sujet_list[0]

    cmap = plt.cm.seismic


    for metric in rf_metrics:

        for band in freq_band_dict:

            global_max = reg_allsujet_data.query(f"sujet == '{sujet}' and rf_metric == '{metric}' and band == '{band}' and coords_info_present == True")['OLS_a'].abs().max()
            # size_max = 20.0 
            size_markers = 10
            size_markers_signi_coeff = 2

            norm = matplotlib.colors.TwoSlopeNorm(vmin=-global_max, vcenter=0.0, vmax=global_max)

            # from nilearn import plotting

            fig, axes = plt.subplots(
                nrows=len(time_phase_list), ncols=len(conditions),
                figsize=(10*len(conditions), 5*len(time_phase_list)),
                constrained_layout=True
            )

            for r, phase in enumerate(time_phase_list):
                for c, cond in enumerate(conditions):

                    ax = axes[r, c]

                    _df = reg_allsujet_data.query(f"sujet == '{sujet}' and rf_metric == '{metric}' and band == '{band}' and cond == '{cond}' and phase_protocol == '{phase}' and coords_info_present == True")
                    coords_all = _df[['coords_x', 'coords_y', 'coords_z']].values

                    colors_all = cmap(norm(_df["OLS_a"].values))

                    signi_circle = (_df["rho_signi"] == 1).to_numpy()

                    coords_signi = coords_all[signi_circle]
                    colors_signi = colors_all[signi_circle]

                    disp = nilearn.plotting.plot_glass_brain(
                        None,
                        display_mode="xz",
                        alpha=0.1,
                        axes=ax
                    )

                    disp.add_markers(
                        coords_signi,
                        marker_color=colors_signi,
                        marker_size=size_markers
                    )
                    

                    ax.set_title(f"{cond} : {phase}")

            plt.suptitle(f"{band}\n{metric}")
            # plt.tight_layout()
            plt.show()

            os.chdir(os.path.join(path_results, 'respi', 'sujet_wise_OLS'))
            fig.savefig(f"{sujet}_{band}_{metric}_OLS.png")



    #### plot allsujet significant contacts only
    metric = 'total_amplitude'
    band = 'gamma'

    cmap = plt.cm.seismic


    for metric in rf_metrics:

        for band in freq_band_dict:

            global_max = reg_allsujet_data.query(f"rf_metric == '{metric}' and band == '{band}' and coords_info_present == True and phase_protocol == 'post'")['OLS_a'].abs().max()
            # size_max = 20.0 
            size_markers = 10
            size_markers_signi_coeff = 2

            norm = matplotlib.colors.TwoSlopeNorm(vmin=-global_max, vcenter=0.0, vmax=global_max)

            # from nilearn import plotting

            fig, axes = plt.subplots(
                nrows=len(time_phase_list), ncols=len(conditions),
                figsize=(10*len(conditions), 5*len(time_phase_list)),
                constrained_layout=True
            )

            for r, phase in enumerate(time_phase_list):
                for c, cond in enumerate(conditions):

                    ax = axes[r, c]

                    _df = reg_allsujet_data.query(f"rf_metric == '{metric}' and band == '{band}' and cond == '{cond}' and coords_info_present == True")
                    coords_all = _df[['coords_x', 'coords_y', 'coords_z']].values

                    colors_all = cmap(norm(_df["OLS_a"].values))

                    signi_circle = (_df["rho_signi"] == 1).to_numpy()

                    coords_signi = coords_all[signi_circle]
                    colors_signi = colors_all[signi_circle]

                    disp = nilearn.plotting.plot_glass_brain(
                        None,
                        display_mode="xz",
                        alpha=0.1,
                        axes=ax
                    )

                    disp.add_markers(
                        coords_signi,
                        marker_color=colors_signi,
                        marker_size=size_markers
                    )
                    

                    ax.set_title(f"{cond} : {phase}")

            plt.suptitle(f"{band}\n{metric}")
            # plt.tight_layout()
            plt.show()

            os.chdir(os.path.join(path_results, 'respi', 'sujet_wise_OLS'))
            fig.savefig(f"{sujet}_{band}_{metric}_OLS.png")


    global_max = reg_allsujet_data.query(f"rf_metric == '{metric}' and band == '{band}' and coords_info_present == True and phase_protocol == 'post'")['OLS_a'].abs().max()
    # size_max = 20.0 
    size_markers = 10
    size_markers_signi_coeff = 2

    norm = matplotlib.colors.TwoSlopeNorm(vmin=-global_max, vcenter=0.0, vmax=global_max)

    # from nilearn import plotting

    fig, axes = plt.subplots(
        ncols=len(conditions),
        figsize=(10*len(conditions), 5*len(time_phase_list)),
        constrained_layout=True
    )

    for c, cond in enumerate(conditions):

        ax = axes[c]

        _df = reg_allsujet_data.query(f"rf_metric == '{metric}' and band == '{band}' and cond == '{cond}' and coords_info_present == True")
        coords_all = _df[['coords_x', 'coords_y', 'coords_z']].values

        colors_all = cmap(norm(_df["OLS_a"].values))

        signi_circle = (_df["rho_signi"] == 1).to_numpy()

        coords_signi = coords_all[signi_circle]
        colors_signi = colors_all[signi_circle]

        disp = nilearn.plotting.plot_glass_brain(
            None,
            display_mode="xz",
            alpha=0.1,
            axes=ax
        )

        disp.add_markers(
            coords_signi,
            marker_color=colors_signi,
            marker_size=size_markers
        )
        

        ax.set_title(f"{cond} : {phase}")

    plt.suptitle(f"{band}\n{metric}")
    # plt.tight_layout()
    plt.show()



################################
######## TF RELATIVE ########
################################

def tf_relative_patientwise():

    #### load
    sujet = sujet_list[0]

    os.chdir(os.path.join(path_precompute, 'TF', 'session'))
    tf_allcond_dict = {}
            
    for cond in conditions:

        tf_allcond_dict[cond] = np.load(f'{sujet}_{cond}_tf_allchan_stretch_post.npy')   

    chanlist, localist_raw = get_chanlist(sujet)
    localist = modify_loca_name(localist_raw)
    localist_unique = np.unique([_loca[0] for _loca in localist if _loca[0].find('UNSORTED') == -1])

    baseline_cond = 'rsp_ctrl'
    test_conditions = [_cond for _cond in conditions if _cond != baseline_cond]
    
    #### plot
    os.chdir(os.path.join(path_results, 'TF', 'patient_wise', 'RELATIVE'))
    time_vec = np.arange(stretch_point_TF)
    ticks_freq = [2, 8, 12, 30, 50, 60, 150]
    percentile_plot = [1,99]

    for loca_sel in localist_unique:

        #### vlim
        _tf_loca_dict = {}
        _tf_loca_vlim = []
        _chan_sel_vec = localist[:,0] == loca_sel

        _tf_baseline = tf_allcond_dict[baseline_cond][_chan_sel_vec]
                
        for cond_i, cond in enumerate(test_conditions):

            _tf_diff = np.median(tf_allcond_dict[cond][_chan_sel_vec], axis=1) - np.median(_tf_baseline, axis=1)
            _tf = np.median(_tf_diff, axis=0)

            if debug:

                plt.pcolormesh(_tf)
                plt.show()

            _tf_loca_vlim.append(_tf)
            _tf_loca_dict[cond] = _tf

        _vlim = np.abs([np.percentile(np.stack(_tf_loca_vlim, axis=0), percentile_plot[0]), np.percentile(np.stack(_tf_loca_vlim, axis=0), percentile_plot[1])]).max()

        fig, axs = plt.subplots(ncols=len(test_conditions), figsize=(16,8))
        
        for cond_i, cond in enumerate(test_conditions):
        
            ax = axs[cond_i]
            pcm = ax.pcolormesh(time_vec, frex, _tf_loca_dict[cond], vmin=-_vlim, vmax=_vlim, cmap='seismic')
            ax.vlines(x=int(stretch_point_TF/2), ymin=frex[0], ymax=frex[-1], color="k", linewidth=2)
            ax.set_yscale('log')
            
            ax.set_yticks(ticks_freq)
            ax.set_yticklabels([str(t) for t in ticks_freq])

            ax.set_title(f"{cond}")

        cbar = fig.colorbar(pcm, ax=axs, orientation="vertical", fraction=0.02, pad=0.04)
        plt.suptitle(f"{sujet} {loca_sel} ({_chan_sel_vec.sum()})")

        # plt.show()
        fig.savefig(f"{sujet}_{loca_sel}.jpg")

    plt.close('all')





def tf_relative_allpatient():

    #### load
    os.chdir(os.path.join(path_results, 'anatomy'))
    df_loca_allpatient = pd.read_excel('anatomy_allpatient.xlsx')
    localist_unique = np.unique([_loca for _loca in df_loca_allpatient['loca'].values if _loca.find('UNSORTED') == -1])

    ticks_freq = [2, 8, 12, 30, 50, 60, 150]
    percentile_plot = [1,99]
    time_vec = np.arange(stretch_point_TF)

    baseline_cond = 'rsp_ctrl'
    test_conditions = [_cond for _cond in conditions if _cond != baseline_cond]

    #loca_sel_i, loca_sel = 2, localist_unique[2]
    # for loca_sel_i, loca_sel in enumerate(localist_unique):
    def plot_allsujet_loca(loca_sel_i, loca_sel):

        print_advancement(loca_sel_i, len(localist_unique), [25,50,75])

        #### load sujet
        os.chdir(os.path.join(path_precompute, 'TF', 'session'))
        loca_sel_vec = df_loca_allpatient['loca'].values == loca_sel
        _patient_list_unique = np.unique(df_loca_allpatient['sujet'].values[loca_sel_vec])
        _tf_allsujet = np.zeros((_patient_list_unique.size, len(test_conditions), frex.size, stretch_point_TF))

        for _sujet_i, _sujet in enumerate(_patient_list_unique):

            _sujet_chanlist, _sujet_localist = get_chanlist(_sujet)
            _sujet_chan_sel_vec = _sujet_localist == loca_sel

            _tf_baseline = np.median(np.load(f'{_sujet}_{baseline_cond}_tf_allchan_stretch_post.npy')[_sujet_chan_sel_vec], axis=1)
                            
            for cond_i, cond in enumerate(test_conditions):

                _tf_cond = np.median(np.load(f'{_sujet}_{cond}_tf_allchan_stretch_post.npy')[_sujet_chan_sel_vec], axis=1)
                _tf_diff = np.median(_tf_cond - _tf_baseline, axis=0)
                _tf_allsujet[_sujet_i, cond_i, :, :] = _tf_diff

        _tf_allsujet_median = np.median(_tf_allsujet, axis=0)

        #### plot
        _vlim = np.abs([np.percentile(np.stack(_tf_allsujet_median, axis=0), percentile_plot[0]), np.percentile(np.stack(_tf_allsujet_median, axis=0), percentile_plot[1])]).max()

        fig, axs = plt.subplots(ncols=len(test_conditions), figsize=(16,5))
        
        for cond_i, cond in enumerate(test_conditions):
        
            ax = axs[cond_i]
            pcm = ax.pcolormesh(time_vec, frex, _tf_allsujet_median[cond_i], vmin=-_vlim, vmax=_vlim, cmap='seismic')
            ax.vlines(x=int(stretch_point_TF/2), ymin=frex[0], ymax=frex[-1], color="k", linewidth=2)
            ax.set_yscale('log')
            
            ax.set_yticks(ticks_freq)
            ax.set_yticklabels([str(t) for t in ticks_freq])

            ax.set_title(f"{cond}")

        cbar = fig.colorbar(pcm, ax=axs, orientation="vertical", fraction=0.02, pad=0.04)
        plt.suptitle(f"{loca_sel} s({len(_patient_list_unique)})")

        # plt.show()

        #### save
        os.chdir(os.path.join(path_results, 'TF', 'allsujet', 'RELATIVE'))
        fig.savefig(f"allsujet_{loca_sel}.jpg")

    joblib.Parallel(n_jobs = n_core, prefer = 'processes')(joblib.delayed(plot_allsujet_loca)(loca_sel_i, loca_sel) for loca_sel_i, loca_sel in enumerate(localist_unique))
    
    plt.close('all')








################################
######## PXX RELATIVE ########
################################

def Pxx_relative_patientwise_allpatient():


    baseline_cond = 'rsp_ctrl'
    test_conditions = [_cond for _cond in conditions if _cond != baseline_cond]

    df_coords_allsujet = get_df_loca_allsujet()

    #### load Pxx
    das = []
    for sujet in sujet_list:
        fp = os.path.join(path_precompute, "TF", "session", f"{sujet}_xr_Pxx.nc")
        da = xr.load_dataarray(fp)
        da = da.expand_dims(sujet=[sujet])

        das.append(da)

    da_Pxx_allsujet_phase_cycle = xr.concat(das, dim="sujet")
    da_Pxx_allsujet_whole_cycle = da_Pxx_allsujet_phase_cycle.median('phase')

    #### load respfeatures
    resp_vars = {} 
    rf_metrics_tot_sel = ['cycle_duration', 'inspi_duration', 'expi_duration', 'cycle_freq', 'inspi_volume',
       'expi_volume', 'total_amplitude', 'inspi_amplitude', 'expi_amplitude',
       'total_volume']
    rf_metrics_whole_cycle = ['cycle_duration', 'cycle_freq', 'total_amplitude', 'total_volume']
    rf_metrics_phase_cycle = ['inspi_duration', 'expi_duration', 'inspi_volume', 'expi_volume', 'inspi_amplitude', 'expi_amplitude']


    for sujet in sujet_list:

        # load respiration features
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

    ds_resp_all = xr.concat(list(resp_vars.values()), dim="sujet").sortby("sujet")

    ds_respfeatures_whole_cycle = ds_resp_all.sel(feature=rf_metrics_whole_cycle)
    ds_respfeatures_phase_cycle = ds_resp_all.sel(feature=rf_metrics_phase_cycle)

    da_Pxx_allsujet_whole_cycle = da_Pxx_allsujet_whole_cycle.to_dataset(name="Pxx")

    ds_respfeatures_whole_cycle = ds_respfeatures_whole_cycle.rename({"phase_group": "pre_post"})
    ds_respfeatures_whole_cycle = ds_respfeatures_whole_cycle.sel(pre_post='post')

    ds_all_whole_cycle = xr.merge([da_Pxx_allsujet_whole_cycle, ds_respfeatures_whole_cycle])

    #### plot
    sujet = sujet_list[0]
    chanlist_sujet = get_chanlist(sujet)[0]
    chan = chanlist_sujet[0]

    ds_plot = ds_all_whole_cycle.sel(sujet=sujet, chan=chan)
    df_plot = ds_plot.to_dataframe().reset_index().dropna()
    df_plot = df_plot.pivot_table(index=[_col for _col in df_plot.columns if _col not in ['feature', 'respfeatures']], columns="feature", values="respfeatures").reset_index()

    sns.lmplot(data=df_plot, x='total_amplitude', y='Pxx', hue='cond', col='band', row='pre_post', sharex=False, sharey=False)
    plt.show()

    
    #### REG STATS WHOLE CYCLE
    os.chdir(os.path.join(path_precompute, 'TF', 'session', 'df_reg'))
    reg_allsujet_data = pd.read_excel(f"df_reg_allsujet.xlsx")
    reg_allsujet_data = reg_allsujet_data.merge(df_coords_allsujet, on=["sujet", "chan"], how="left")
    
    reg_allsujet_data_diff = []
    reg_allsujet_baseline = reg_allsujet_data.query(f"phase_protocol == 'post' and cond == 'rsp_ctrl'")

    for cond in test_conditions:

        _df = reg_allsujet_data.query(f"phase_protocol == 'post' and cond == '{cond}'")
        _df_diff = _df.copy()
        _df_diff['OLS_a'] = _df['OLS_a'].values - reg_allsujet_baseline['OLS_a'].values
        
        reg_allsujet_data_diff.append(_df_diff)

    reg_allsujet_data_diff = pd.concat(reg_allsujet_data_diff)

    #### plot sujet wise all contacts
    metric = 'total_amplitude'
    band = 'gamma'
    sujet = sujet_list[0]

    cmap = plt.cm.seismic
    rf_metrics = rf_metrics_whole_cycle
    size_markers = 10

    for sujet in sujet_list:

        for band in freq_band_dict:

            nrows = len(rf_metrics)
            ncols = len(test_conditions)

            fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10 * ncols, 5 * nrows), constrained_layout=True)

            for r, metric in enumerate(rf_metrics):

                _df_all = reg_allsujet_data_diff.query(f"sujet == '{sujet}' and rf_metric == '{metric}' and band == '{band}' and coords_info_present == True")
                global_max = _df_all["OLS_a"].abs().max()

                norm = matplotlib.colors.TwoSlopeNorm(vmin=-global_max, vcenter=0.0, vmax=global_max)

                for c, cond in enumerate(test_conditions):

                    ax = axes[r, c]

                    _df = reg_allsujet_data_diff.query(f"sujet == '{sujet}' and rf_metric == '{metric}' and band == '{band}' and cond == '{cond}' and coords_info_present == True")

                    coords_all = _df[['coords_x', 'coords_y', 'coords_z']].values
                    colors_all = cmap(norm(_df["OLS_a"].values))

                    signi_circle = (_df["rho_signi"] == 1).to_numpy()
                    coords_no_signi = coords_all[~signi_circle]
                    colors_no_signi = colors_all[~signi_circle]

                    disp = plotting.plot_glass_brain(
                        None,
                        display_mode="xz",
                        alpha=0.1,
                        axes=ax
                    )

                    disp.add_markers(
                        coords_no_signi,
                        marker_color=list(colors_no_signi),
                        marker_size=size_markers
                    )

                    if r == 0:
                        ax.set_title(f"{cond}")

                    if c == 0:
                        ax.text(
                            -0.05, 0.5, metric,
                            transform=ax.transAxes,
                            rotation=90,
                            va="center",
                            ha="right",
                            fontsize=12
                        )

            plt.suptitle(f"{sujet} — {band}", fontsize=16)

            # plt.show()

            outdir = os.path.join(path_results, 'Pxx', 'patient_wise', 'RELATIVE')
            fig.savefig(os.path.join(outdir, f"allplot_DIFF_{sujet}_{band}_OLS.png"), dpi=200)
            plt.close(fig)





    #### plot allsujet contacts only
    band = 'gamma'

    cmap = plt.cm.seismic
    rf_metrics = rf_metrics_whole_cycle
    size_markers = 10

    for band in freq_band_dict:

        nrows = len(rf_metrics)
        ncols = len(test_conditions)

        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10 * ncols, 5 * nrows), constrained_layout=True)

        for r, metric in enumerate(rf_metrics):

            _df_all = reg_allsujet_data_diff.query(f"rf_metric == '{metric}' and band == '{band}' and coords_info_present == True")
            global_max = _df_all["OLS_a"].abs().max()

            norm = matplotlib.colors.TwoSlopeNorm(vmin=-global_max, vcenter=0.0, vmax=global_max)

            for c, cond in enumerate(test_conditions):

                ax = axes[r, c]

                _df = reg_allsujet_data_diff.query(f"rf_metric == '{metric}' and band == '{band}' and cond == '{cond}' and coords_info_present == True")

                coords_all = _df[['coords_x', 'coords_y', 'coords_z']].values
                colors_all = cmap(norm(_df["OLS_a"].values))

                signi_circle = (_df["rho_signi"] == 1).to_numpy()
                coords_no_signi = coords_all[~signi_circle]
                colors_no_signi = colors_all[~signi_circle]

                disp = plotting.plot_glass_brain(
                    None,
                    display_mode="xz",
                    alpha=0.1,
                    axes=ax
                )

                disp.add_markers(
                    coords_no_signi,
                    marker_color=list(colors_no_signi),
                    marker_size=size_markers
                )

                if r == 0:
                    ax.set_title(f"{cond}")

                if c == 0:
                    ax.text(
                        -0.05, 0.5, metric,
                        transform=ax.transAxes,
                        rotation=90,
                        va="center",
                        ha="right",
                        fontsize=12
                    )

        plt.suptitle(f"{band}", fontsize=16)

        # plt.show()

        outdir = os.path.join(path_results, 'Pxx', 'allsujet', 'RELATIVE')
        fig.savefig(os.path.join(outdir, f"allsujet_DIFF_{band}_OLS.png"), dpi=200)
        plt.close(fig)


    
    #### plot allsujet allROI
    band = 'gamma'

    rf_metrics = rf_metrics_whole_cycle

    for band in freq_band_dict:

        _df_all = reg_allsujet_data_diff.query(f"rf_metric == '{metric}' and band == '{band}' and coords_info_present == True")


        ROI_signi_sel = reg_allsujet_data_diff.query(f"band == '{band}' and rf_metric == '{metric}' and loca == 'Amygdala'")["ROI"].values
        _df_plot = df_R_Pxx.query(f"band == '{band}' and term != '(Intercept)' and ROI in {ROI_signi_sel.tolist()}").copy()

        _df_plot["sig"] = _df_plot["pvalue"].apply(p_to_stars)

        fig = px.bar(_df_plot, x="ROI", y="estimate", color="term", barmode="group", text="sig", title=band)

        fig.update_layout(xaxis_tickangle=-45, yaxis_title="estimate", xaxis_title="ROI", legend_title="term")
        fig.update_traces(textposition="outside", cliponaxis=False)

        # fig.show()

        outdir = os.path.join(path_results, "LMM", "fig")
        fig.write_html(os.path.join(outdir, f"{band}_signiROI_barplot_LMM.html"),
                   include_plotlyjs="cdn")
        





        nrows = len(rf_metrics)
        ncols = len(test_conditions)

        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10 * ncols, 5 * nrows), constrained_layout=True)

        for r, metric in enumerate(rf_metrics):

            _df_all = reg_allsujet_data_diff.query(f"rf_metric == '{metric}' and band == '{band}' and coords_info_present == True")
            global_max = _df_all["OLS_a"].abs().max()

            norm = matplotlib.colors.TwoSlopeNorm(vmin=-global_max, vcenter=0.0, vmax=global_max)

            for c, cond in enumerate(test_conditions):

                ax = axes[r, c]

                _df = reg_allsujet_data_diff.query(f"rf_metric == '{metric}' and band == '{band}' and cond == '{cond}' and coords_info_present == True")

                coords_all = _df[['coords_x', 'coords_y', 'coords_z']].values
                colors_all = cmap(norm(_df["OLS_a"].values))

                signi_circle = (_df["rho_signi"] == 1).to_numpy()
                coords_no_signi = coords_all[~signi_circle]
                colors_no_signi = colors_all[~signi_circle]

                disp = plotting.plot_glass_brain(
                    None,
                    display_mode="xz",
                    alpha=0.1,
                    axes=ax
                )

                disp.add_markers(
                    coords_no_signi,
                    marker_color=list(colors_no_signi),
                    marker_size=size_markers
                )

                if r == 0:
                    ax.set_title(f"{cond}")

                if c == 0:
                    ax.text(
                        -0.05, 0.5, metric,
                        transform=ax.transAxes,
                        rotation=90,
                        va="center",
                        ha="right",
                        fontsize=12
                    )

        plt.suptitle(f"{band}", fontsize=16)

        # plt.show()

        outdir = os.path.join(path_results, 'Pxx', 'allsujet', 'RELATIVE')
        fig.savefig(os.path.join(outdir, f"allsujet_DIFF_{band}_OLS.png"), dpi=200)
        plt.close(fig)



























################################
######## EXECUTE ########
################################


if __name__ == '__main__':


    export_allpatient_anat()
    
    
    tf_absolute_patientwise()
    tf_absolute_allpatient()
    Pxx_absolute_patientwise()


    tf_relative_patientwise()
    tf_relative_allpatient()

    Pxx_relative_patientwise_allpatient()

    
    






















    


    


















