



from n00_config_params import *
from n00bis_config_analysis_functions import *
from n01_manip_data import *

import nilearn







################################
######## EXECUTE ########
################################


if __name__ == '__main__':


    for cond in conditions:
    
        for sujet in sujet_list_allcond[cond]:

            raise

    os.chdir(os.path.join(path_precompute, 'TF'))
    df_Pxx = pd.read_excel(f'{sujet}_{cond}_df_Pxx.xlsx')    

    os.chdir(os.path.join(path_precompute, 'TF'))
    tf_cond = np.load(f'{sujet}_{cond}_tf_conv.npy')   

    time_vec_sec = np.arange(tf_cond[0,0].size)/srate - (tf_cond[0,0].size/srate)/2

    plot_time = [-4, 4]
    plot_vec_sel = (time_vec_sec >= plot_time[0]) & (time_vec_sec <= plot_time[1])
    plot_vec_sec = time_vec_sec[plot_vec_sel]

    chan_list, loca_list = get_chanlist(sujet)

    data_allcond, resp_allcond, chanlist, localist = get_data_sujet(sujet)
    data, resp = data_allcond[cond], resp_allcond[cond]
    del data_allcond, resp_allcond
    resp_med = np.median(resp, axis=0)[plot_vec_sel]

    for chan_i in range(10):

        fig, ax1 = plt.subplots()

        pcm = ax1.pcolormesh(plot_vec_sec, frex, tf_cond[chan_i][:, plot_vec_sel])
        ax1.set_yscale('log')

        ticks = [2, 8, 12, 30, 50, 60, 150]
        ax1.set_yticks(ticks)
        ax1.set_yticklabels([str(t) for t in ticks])

        ax1.set_xlabel('Time (s)')
        ax1.set_ylabel('Frequency (Hz)')
        ax1.set_title(f"{sujet} {cond}\n{chan_list[chan_i]} {loca_list[chan_i]}")

        # fig.colorbar(pcm, ax=ax1, label='Power')

        ax2 = ax1.twinx()          
        ax2.plot(plot_vec_sec, resp_med, color='r')  
        ax2.set_ylabel('Resp (a.u.)')     

        fig.tight_layout()

        plt.show()







    #### all df

    df_Pxx_allsujet = []

    for cond in conditions:
    
        for sujet in sujet_list_allcond[cond]:

            print(f"{sujet} {cond}")

            print('load')
            os.chdir(os.path.join(path_precompute, 'TF'))
            _df_Pxx = pd.read_excel(f'{sujet}_{cond}_df_Pxx.xlsx')    

            os.chdir(os.path.join(path_precompute, 'RESP')) 
            respfeatures_pre = pd.read_excel(f'{sujet}_{cond}_df_respfeatures_pre.xlsx')
            respfeatures_post = pd.read_excel(f'{sujet}_{cond}_df_respfeatures_post.xlsx')

            print('aggregate')
            _df_respfeature_all = []

            for i, r in _df_Pxx.iterrows():

                _cycle, _phase = r['cycle'], r['phase']
                
                if _phase in ['pre_inspi', 'pre_expi']:
                    _df_respfeature_all.append(respfeatures_pre.iloc[_cycle])
                elif _phase in ['inspi', 'expi']:
                    _df_respfeature_all.append(respfeatures_post.iloc[_cycle])

            print('concat')
            _df_allconcatenated = pd.concat([_df_Pxx, pd.DataFrame(_df_respfeature_all).set_index(_df_Pxx.index)], axis=1)
            _df_allconcatenated = _df_allconcatenated.drop(columns='Unnamed: 0')

            df_Pxx_allsujet.append(_df_allconcatenated)

    df_Pxx_allsujet = pd.concat(df_Pxx_allsujet, ignore_index=True, sort=False)

    rf_metrics = ['cycle_duration', 'inspi_duration', 'expi_duration', 'cycle_freq', 'inspi_volume',
       'expi_volume', 'total_amplitude', 'inspi_amplitude', 'expi_amplitude',
       'total_volume']

    df_plot = df_Pxx_allsujet.query(f"sujet == '{sujet_list[0]}' and chan == 'LDh1'")
    sns.lmplot(data=df_plot, x='total_amplitude', y='Pxx', hue='phase', col='band')

    #### channel wise
    df_plot = df_Pxx_allsujet.query(f"sujet == '{sujet_list[0]}' and chan == 'LDh1'")
    for _metric in rf_metrics:
        sns.lmplot(data=df_plot, x=f'{_metric}', y='Pxx', hue='cond', col='band', row='phase', sharex=False, sharey=False)
        plt.show()

    #### FOR ONE SUJET
        #### generate df_reg
    _sujet = sujet_list[0]

    df_sujet = df_Pxx_allsujet.query(f"sujet == '{_sujet}'")
    df_reg_list = []

    for _cond in df_sujet['cond'].unique():

        print(_cond)

        for _chan_i, _chan in enumerate(df_sujet['chan'].unique()):

            print_advancement(_chan_i, df_sujet['chan'].unique().shape[0], steps=[25,50,75])

            for _band in df_sujet['band'].unique():

                for _phase in df_sujet['phase'].unique():

                    _df = df_sujet.query(f"cond == '{_cond}' and chan == '{_chan}' and band == '{_band}' and phase == '{_phase}'")

                    _ROI = _df['ROI'].values[0]

                    for _metric in rf_metrics:

                        import statsmodels.formula.api as smf

                        model = smf.ols(f"Pxx ~ {_metric}", data=_df).fit()
                        _coef = model.params[_metric]       # pandas Series
                        _pval = model.pvalues[_metric]

                        df_reg_list.append(pd.DataFrame({'sujet' : [_sujet], 'cond' : [_cond],  'chan' : [_chan],  'ROI' : [_ROI],  'band' : [_band],  'phase' : [_phase], 
                                                    'metric' : [_metric], 'coef' : [_coef], 'pval' : [_pval]}))
                                                
    df_reg = pd.concat(df_reg_list)

        #### phase / cond
    rf_metrics = ['cycle_duration', 'inspi_duration', 'expi_duration', 'cycle_freq', 'inspi_volume',
       'expi_volume', 'total_amplitude', 'inspi_amplitude', 'expi_amplitude',
       'total_volume']
    
    metric = 'total_amplitude'
    band = 'gamma'

    for metric in rf_metrics:

        for band in df_sujet['band'].unique():

            df_plot = df_reg.query(f"metric == '{metric}' and band == '{band}'")
            conds = df_plot['cond'].unique()
            phases = df_plot['phase'].unique()

            global_max = df_plot['coef'].abs().max()
            size_max   = 20.0 

            fig, axes = plt.subplots(
                nrows=len(phases), ncols=len(conds),
                figsize=(10*len(conds), 3*len(phases)),
                constrained_layout=True
            )

            if len(phases) == 1 and len(conds) == 1:
                axes = np.array([[axes]])
            elif len(phases) == 1:
                axes = axes[np.newaxis, :]
            elif len(conds) == 1:
                axes = axes[:, np.newaxis]

            for r, phase in enumerate(phases):
                for c, cond in enumerate(conds):

                    ax = axes[r, c]

                    _df = df_plot.query("cond == @cond and phase == @phase")
                    if _df.empty:
                        ax.axis('off')
                        ax.set_title(f"{cond}\n{phase}")
                        continue

                    colors = np.where(_df["pval"].to_numpy() < 0.05, "crimson", "royalblue")

                    sizes = _df["coef"].abs().to_numpy()
                    sizes_scaled = (sizes / global_max) * size_max
                    sizes_scaled = np.clip(sizes_scaled, 2, None)

                    disp = nilearn.plotting.plot_glass_brain(
                        None,
                        display_mode="lyrz",
                        alpha=0.1,
                        axes=ax
                    )
                    disp.add_markers(_sujet_coords, marker_color=list(colors), marker_size=list(sizes_scaled))

                    ax.set_title(f"{cond} : {phase}")

            plt.suptitle(f"{band}\n{metric}")
            plt.tight_layout()
            # plt.show()

            os.chdir(os.path.join(path_results, 'respi', 'sujet_wise_OLS'))
            fig.savefig(f"{sujet}_{band}_{metric}_OLS.png")

















    #### 2D
    df_base = df_reg.query("cond == 'oc_ctrl' and metric == 'cycle_duration'")
    bands  = df_base['band'].unique()
    phases = df_base['phase'].unique()

    global_max = df_base['coef'].abs().max()
    size_max   = 20.0 


    fig, axes = plt.subplots(
        nrows=len(phases), ncols=len(bands),
        figsize=(4*len(bands), 3*len(phases)),
        constrained_layout=True
    )

    if len(phases) == 1 and len(bands) == 1:
        axes = np.array([[axes]])
    elif len(phases) == 1:
        axes = axes[np.newaxis, :]
    elif len(bands) == 1:
        axes = axes[:, np.newaxis]

    for r, phase in enumerate(phases):
        for c, band in enumerate(bands):

            ax = axes[r, c]

            df_plot = df_base.query("band == @band and phase == @phase")
            if df_plot.empty:
                ax.axis('off')
                ax.set_title(f"{band}\n{phase}")
                continue

            colors = np.where(df_plot["pval"].to_numpy() < 0.05, "crimson", "royalblue")

            sizes = df_plot["coef"].abs().to_numpy()
            sizes_scaled = (sizes / global_max) * size_max
            sizes_scaled = np.clip(sizes_scaled, 2, None)

            disp = plotting.plot_glass_brain(
                None,
                display_mode="lyrz",
                alpha=0.1,
                axes=ax
            )
            disp.add_markers(_sujet_coords, marker_color=list(colors), marker_size=list(sizes_scaled))

            ax.set_title(f"{band}\n{phase}")

    plt.show()















    loca_list_allsujet = df_Pxx_allsujet['ROI'].values

    loca_list_corr = modify_loca_name(loca_list_allsujet)
    df_Pxx_allsujet['ROI'] = loca_list_corr

    df_Pxx_allsujet = df_Pxx_allsujet.drop(columns=['Unnamed: 0'])
    df_Pxx_allsujet = df_Pxx_allsujet.reset_index(drop=True)

    # loca_list_filtered = np.array([_loca for _loca in df_Pxx_allsujet['ROI'].unique() if _loca.find('Unknown') == -1 and _loca.find('REF') == -1 and _loca.find('Ref') == -1])

    loca_list_filtered_allsujet = np.array(['Amygdala', 'White-Matter', 'insula-ant', 'parsopercularis',
       'parahippocampal', 'inferiorparietal', 'inferiortemporal',
       'superiortemporal', 'fusiform', 'parstriangularis',
       'rostralmiddlefrontal', 'superiorparietal',
       'rostralanteriorcingulate', 'lateralorbitofrontal',
       'middletemporal', 'Hippocampus', 'temporalpole',
       'insula-pos', 'lateraloccipital', 'lingual', 'precuneus',
       'isthmuscingulate', 'cuneus', 'pericalcarine', 'supramarginal',
       'transversetemporal', 'medialorbitofrontal', 'superiorfrontal',
       'entorhinal', 'posteriorcingulate', 'insula', 'postcentral', 'Thalamus', 'caudalmiddlefrontal',
       'paracentral', 'precentral', 'Brain-Stem', 'Putamen',
       'caudalanteriorcingulate', 'Accumbens-area', 'Caudate',
       'superiortempora'])
    
    #_loca = loca_list_filtered_allsujet[0]
    for _loca in loca_list_filtered_allsujet:
        
        df_plot = df_Pxx_allsujet.query(f"ROI == '{_loca}'")
        df_plot = df_plot.dropna()
        sns.catplot(data=df_plot, kind='box', x='phase', y='Pxx', hue='cond', col='band', showfliers=False)
        plt.suptitle(_loca)
        plt.show()

    #### EXPORT TF

    plot_time = [-4, 4]

    for cond in conditions:

        for _loca in loca_list_filtered_allsujet[:10]:

            sujet_removed = 0

            for sujet_i, sujet in enumerate(sujet_list_allcond[cond]):

                print(f"{sujet} {cond}")

                os.chdir(os.path.join(path_precompute, 'TF'))
                _tf = np.load(f'{sujet}_{cond}_tf_conv.npy')    

                _chan_list, _loca_list = get_chanlist(sujet)
                _loca_list_filtered = modify_loca_name(_loca_list)
                _loca_sel = _loca_list_filtered == _loca

                _half_time = _tf.shape[-1] * 1/srate - (_tf.shape[-1] * 1/srate)/2
                time_vec_sec = np.linspace(-_half_time, _half_time, _tf.shape[-1])
                plot_vec_sel = (time_vec_sec >= plot_time[0]) & (time_vec_sec <= plot_time[1])
                plot_vec_sec = time_vec_sec[plot_vec_sel]

                _tf_filtered_loca = np.median(_tf[_loca_sel], axis=0)
                _tf_filtered = _tf_filtered_loca[:,plot_vec_sel].reshape(1, nfrex, -1)

                if np.isnan(_tf_filtered).sum() != 0:
                    sujet_removed += 1
                    continue

                if sujet_i == 0:
                    tf_allsujet = _tf_filtered
                else:
                    tf_allsujet = np.concatenate([tf_allsujet, _tf_filtered], axis=0)

            tf_med = np.median(tf_allsujet, axis=0)
            vlim = np.abs([tf_med.min(), tf_med.max()]).max()

            fig, ax1 = plt.subplots()

            pcm = ax1.pcolormesh(plot_vec_sec, frex, tf_med, vmin=-vlim, vmax=vlim, cmap='seismic')
            ax1.set_yscale('log')

            ticks = [2, 8, 12, 30, 50, 60, 150]
            ax1.set_yticks(ticks)
            ax1.set_yticklabels([str(t) for t in ticks])

            ax1.set_xlabel('Time (s)')
            ax1.set_ylabel('Frequency (Hz)')
            ax1.set_title(f"{_loca} \n {cond} \n nsujet:{tf_allsujet.shape[0]} / rmv sujet:{sujet_removed}")

            fig.colorbar(pcm, ax=ax1, label='Power')

            fig.tight_layout()

            # plt.show()

            os.chdir(os.path.join(path_results, 'TF', 'allsujet'))
            fig.savefig(f"{_loca}_{cond}.png")

            plt.close('all')

    #### EXPORT TF DIFF

    plot_time = [-4, 4]

    for cond in conditions:

        for _loca in loca_list_filtered_allsujet[:10]:

            sujet_removed = 0

            for sujet_i, sujet in enumerate(sujet_list_allcond[cond]):

                print(f"{sujet} {cond}")

                os.chdir(os.path.join(path_precompute, 'TF'))
                _tf = np.load(f'{sujet}_{cond}_tf_conv.npy')    

                _chan_list, _loca_list = get_chanlist(sujet)
                _loca_list_filtered = modify_loca_name(_loca_list)
                _loca_sel = _loca_list_filtered == _loca

                _half_time = _tf.shape[-1] * 1/srate - (_tf.shape[-1] * 1/srate)/2
                time_vec_sec = np.linspace(-_half_time, _half_time, _tf.shape[-1])
                plot_vec_sel = (time_vec_sec >= plot_time[0]) & (time_vec_sec <= plot_time[1])
                plot_vec_sec = time_vec_sec[plot_vec_sel]

                _tf_filtered_loca = np.median(_tf[_loca_sel], axis=0)
                _tf_filtered = _tf_filtered_loca[:,plot_vec_sel].reshape(1, nfrex, -1)

                if np.isnan(_tf_filtered).sum() != 0:
                    sujet_removed += 1
                    continue

                if sujet_i == 0:
                    tf_allsujet = _tf_filtered
                else:
                    tf_allsujet = np.concatenate([tf_allsujet, _tf_filtered], axis=0)

            tf_med = np.median(tf_allsujet, axis=0)
            vlim = np.abs([tf_med.min(), tf_med.max()]).max()

            fig, ax1 = plt.subplots()

            pcm = ax1.pcolormesh(plot_vec_sec, frex, tf_med, vmin=-vlim, vmax=vlim, cmap='seismic')
            ax1.set_yscale('log')

            ticks = [2, 8, 12, 30, 50, 60, 150]
            ax1.set_yticks(ticks)
            ax1.set_yticklabels([str(t) for t in ticks])

            ax1.set_xlabel('Time (s)')
            ax1.set_ylabel('Frequency (Hz)')
            ax1.set_title(f"{_loca} \n {cond} \n nsujet:{tf_allsujet.shape[0]} / rmv sujet:{sujet_removed}")

            fig.colorbar(pcm, ax=ax1, label='Power')

            fig.tight_layout()

            # plt.show()

            os.chdir(os.path.join(path_results, 'TF', 'allsujet'))
            fig.savefig(f"{_loca}_{cond}.png")

            plt.close('all')








    






