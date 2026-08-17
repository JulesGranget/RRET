



from statsmodels.stats.multitest import multipletests
import itertools
from A_config.n01_O_params import *
from A_config.n02_O_analysis_functions import *
from A_config.n03_X_manip_data import *






########################################
######## EXPORT OC STATS ########
########################################

def export_challenge_stats():

    #### params
    rf_metric_sel = ['cycle_duration', 'inspi_duration', 'expi_duration', 'cycle_freq', 'cycle_ratio',
       'inspi_volume', 'expi_volume', 'total_amplitude', 'inspi_amplitude', 'expi_amplitude', 'total_volume']

    #### get data
    df_chl_allsujet = []
    
    for sujet in sujet_list:

        path_respfeatures = os.path.join(path_precompute, 'RESP', 'respfeatures')
        filename = f"{sujet}_respfeatures_cleaned_label.xlsx"
        _df_rf = pd.read_excel(os.path.join(path_respfeatures, filename))
        _oc_r = rscore(_df_rf[rf_metric_sel].values, median=None, mad=None, axis=0)
        _df_rf_r = pd.DataFrame(_oc_r, columns=rf_metric_sel)
        _df_rf_r = pd.concat([_df_rf[['sujet', 'cond']], _df_rf_r], axis=1)
        df_chl_allsujet.append(_df_rf_r)

    df_chl_allsujet = pd.concat(df_chl_allsujet)
    df_median_chl = df_chl_allsujet.groupby(['sujet', 'cond']).median().reset_index()[['sujet', 'cond'] + rf_metric_sel]

    #### median
    df_stats_chl = []

    for metric_rf in rf_metric_sel:
        
        df_pivot = df_median_chl.pivot(index='sujet', columns='cond', values=metric_rf)
        statfried, p_fried = scipy.stats.friedmanchisquare(df_pivot['oc_chl'], df_pivot['oc_ctrl'], df_pivot['rsp_chl'], df_pivot['rsp_ctrl'])

        results = []

        for c1, c2 in itertools.combinations(conditions, 2):
            
            stat, p = scipy.stats.wilcoxon(df_pivot[c1], df_pivot[c2])
        
            results.append({"cond1": c1, "cond2": c2, "stat": stat, "p": p})

        posthoc = pd.DataFrame(results)
        posthoc["p_corrected"] = multipletests(posthoc["p"], method="bonferroni")[1]
        posthoc['friedman_p'] = [p_fried] * posthoc.shape[0]
        posthoc['friedman_statval'] = [statfried] * posthoc.shape[0]
        posthoc['rf_metric'] = [metric_rf] * posthoc.shape[0]

        df_stats_chl.append(posthoc)

    df_stats_chl = pd.concat(df_stats_chl)

    #### plot
    for metric_rf in rf_metric_sel:

        fig_median_chl, ax = plt.subplots(figsize=(10,8))

        sns.swarmplot(
            data=df_median_chl,
            x='cond',
            y=metric_rf,
            hue='sujet',
            size=10,
            ax=ax
        )

        # plt.show()

        path_savefig = os.path.join(path_results, 'respi', 'rf_analysis')
        filename_savefig = f"{metric_rf}_rscore.png"
        fig_median_chl.savefig(os.path.join(path_savefig, filename_savefig))

    filename_savedf = f"df_stat_rf_metric.xlsx"
    df_stats_chl.to_excel(os.path.join(path_savefig, filename_savedf))







################################
######## EXPORT PLOT ########
################################

#sujet = sujet_list[0]
def export_plot_oc_resp():

    sujet = sujet_list[0]

    print(sujet)

    #### extract
    os.chdir(os.path.join(path_precompute, 'RESP', 'session'))

    resp_stretch = {}

    for cond in conditions:
        resp_stretch[cond] = np.load(f"{sujet}_{cond}_stretch_resp_post.npy")
    
    oc_ratio = {}
    oc_val = {}

    for cond in conditions:

        oc_ratio[cond] = []
        oc_val[cond] = []

        for cycle_i in range(resp_stretch[cond].shape[0]): 

            sig = resp_stretch[cond][cycle_i]
            sig_diff = np.diff(sig)
            start_inspi_dec = np.argmin(sig_diff[:int(stretch_point_TF/4)])
            start_inspi_asc = np.argmax(sig_diff[int(stretch_point_TF/4):int(stretch_point_TF/2)]) + int(stretch_point_TF/4)

            oc_top_i = np.where(sig_diff[start_inspi_dec:] > 0)[0][0] + start_inspi_dec
            oc_top_val = sig[oc_top_i]
            median_inspi_trough = np.median(sig[start_inspi_dec:start_inspi_asc])

            _oc_ratio = (oc_top_val - median_inspi_trough) / median_inspi_trough
            oc_ratio[cond].append(_oc_ratio)
            oc_val[cond].append(oc_top_val)

            # if _oc_ratio > 10 and cond == 'rsp_chl':
            #     raise

            if debug:

                time_vec = np.arange(sig.size)[:-1]
                plt.plot(time_vec, sig[:-1])
                plt.vlines([start_inspi_asc, start_inspi_dec], ymin=sig.min(), ymax=sig.max(), colors='g')
                plt.hlines([median_inspi_trough], xmin=0, xmax=stretch_point_TF, colors='b')
                plt.scatter([oc_top_i], sig[oc_top_i], color='r')
                plt.plot(time_vec,sig_diff)
                plt.show()

    if debug:         

        cond = 'oc_ctrl'
        cond = 'oc_chl'

        for cycle_i in range(resp_stretch[cond].shape[0]): 

            sig = resp_stretch[cond][cycle_i][:-1]
            sig_diff = np.diff(resp_stretch[cond][cycle_i])
            start_inspi_dec = np.argmin(sig_diff[:int(stretch_point_TF/4)])
            start_inspi_asc = np.argmax(sig_diff[int(stretch_point_TF/4):int(stretch_point_TF/2)]) + int(stretch_point_TF/4)

            oc_top_i = np.where(sig_diff[start_inspi_dec:] > 0)[0][0] + start_inspi_dec

            median_inspi_trough = np.median(sig[start_inspi_dec:start_inspi_asc])

            time_vec = np.arange(sig.size)
            plt.plot(time_vec, sig)
            plt.vlines([start_inspi_asc, start_inspi_dec], ymin=sig.min(), ymax=sig.max(), colors='g')
            plt.hlines([median_inspi_trough], xmin=0, xmax=stretch_point_TF, colors='b')
            plt.scatter([oc_top_i], sig[oc_top_i], color='r')
            plt.plot(time_vec,sig_diff)
            plt.title(f"{sujet} {cond} {cycle_i}")
            plt.show()

    #### plot RB_OC example
    cycle_i_RB = 0
    cond = 'oc_ctrl'

    sig = resp_stretch[cond][cycle_i_RB][:-1]
    sig_diff = np.diff(resp_stretch[cond][cycle_i_RB])
    start_inspi_dec = np.argmin(sig_diff[:int(stretch_point_TF/4)])
    start_inspi_asc = np.argmax(sig_diff[int(stretch_point_TF/4):int(stretch_point_TF/2)]) + int(stretch_point_TF/4)

    oc_top_i = np.where(sig_diff[start_inspi_dec:] > 0)[0][0] + start_inspi_dec

    median_inspi_trough = np.median(sig[start_inspi_dec:start_inspi_asc])

    time_vec = np.arange(sig.size)

    fig_oc_example, ax = plt.subplots()
    ax.plot(time_vec, sig)
    ax.hlines([median_inspi_trough], xmin=0, xmax=stretch_point_TF, colors='b')
    # ax.hlines([0], xmin=0, xmax=stretch_point_TF, colors='r')
    ax.vlines([stretch_point_TF/2], ymin=sig.min()+0.1*sig.min(), ymax=sig.max()+0.1*sig.max(), colors='k')

    ax.tick_params(axis="x", labelsize=12)
    ax.tick_params(axis="y", labelsize=12)
    ax.scatter([oc_top_i], sig[oc_top_i], color='r')
    ax.set_title(f"{sujet} {cond} {cycle_i_RB}")
    ax.set_ylabel(f"r-zscore", fontsize=14)
    ax.set_xlabel(f"Phase", fontsize=14)

    plt.ylim(sig.min()+0.1*sig.min(),sig.max()+0.1*sig.max())

    # plt.show()

    path_savefig = os.path.join(path_results, 'respi', 'oc_ratio', 'plot')
    filename_savefig = f"OC_example_OCCTRL.png"
    fig_oc_example.savefig(os.path.join(path_savefig, filename_savefig))

    filename_savefig = f"fig02c_OC_example_OCCTRL.svg"
    fig_oc_example.savefig(os.path.join(path_paper_figure_export, filename_savefig))

    #### plot OC_CHL
    cycle_i_OC = 25
    cond = 'oc_chl'

    sig = resp_stretch[cond][cycle_i_OC][:-1]
    sig_diff = np.diff(resp_stretch[cond][cycle_i_OC])
    start_inspi_dec = np.argmin(sig_diff[:int(stretch_point_TF/4)])
    start_inspi_asc = np.argmax(sig_diff[int(stretch_point_TF/4):int(stretch_point_TF/2)]) + int(stretch_point_TF/4)

    oc_top_i = np.where(sig_diff[start_inspi_dec:] > 0)[0][0] + start_inspi_dec

    median_inspi_trough = np.median(sig[start_inspi_dec:start_inspi_asc])

    time_vec = np.arange(sig.size)

    fig_oc_example, ax = plt.subplots()
    ax.plot(time_vec, sig)
    ax.hlines([median_inspi_trough], xmin=0, xmax=stretch_point_TF, colors='b')
    # ax.hlines([0], xmin=0, xmax=stretch_point_TF, colors='r')
    ax.vlines([stretch_point_TF/2], ymin=sig.min()+0.1*sig.min(), ymax=sig.max()+0.1*sig.max(), colors='k')

    ax.tick_params(axis="x", labelsize=12)
    ax.tick_params(axis="y", labelsize=12)
    ax.scatter([oc_top_i], sig[oc_top_i], color='r')
    ax.set_title(f"{sujet} {cond} {cycle_i_OC}")
    ax.set_ylabel(f"r-zscore", fontsize=14)
    ax.set_xlabel(f"Phase", fontsize=14)

    plt.ylim(sig.min()+0.1*sig.min(),sig.max()+0.1*sig.max())

    # plt.show()

    path_savefig = os.path.join(path_results, 'respi', 'oc_ratio', 'plot')
    filename_savefig = f"OC_example_OCCHL.png"
    fig_oc_example.savefig(os.path.join(path_savefig, filename_savefig))

    filename_savefig = f"fig02c_OC_example_OCCHL.svg"
    fig_oc_example.savefig(os.path.join(path_paper_figure_export, filename_savefig))

    




def export_Pxx_diff_with_patientwise_info():

    #### load
    df_Pxx_allband_raw = []
    
    os.chdir(os.path.join(path_precompute, 'TF', 'session', 'df_R'))
    
    for band in freq_band_dict:
        df_Pxx_allband_raw.append(pd.read_excel(f"df_R_{band}.xlsx").drop(columns=['Unnamed: 0']).query(f"phase_cycle != 'whole' and pre_post == 'post'"))

    df_Pxx_allband_raw = pd.concat(df_Pxx_allband_raw)

    #### reduce
    df_Pxx_allband_raw['cond'] = df_Pxx_allband_raw['resp'] + '_' + df_Pxx_allband_raw['state']
    
    df_Pxx_allband = df_Pxx_allband_raw.groupby(['chan', 'phase_cycle', 'band', 'pre_post', 'ROI', 'sujet', 'resp', 'state', 'cond']).median().reset_index().drop(columns=['cycle', 'pre_post', 'resp', 'state'])
    df_Pxx_allband = df_Pxx_allband.drop(columns=['chan']).groupby(['phase_cycle', 'band', 'ROI', 'sujet', 'cond']).median().reset_index()

    #### plot
    os.chdir(os.path.join(path_results, 'Pxx', 'allsujet', 'allpatient_patientlinked'))

    example_plot = [['theta', 'Hippocampus'], ['beta', 'lateralorbitofrontal'], ['gamma', 'Amygdala']]

    #band = 'theta'
    for band in freq_band_dict:

        #ROI_sel = 'Hippocampus'
        for ROI_sel in ROI_short_list:

            df_plot = df_Pxx_allband.query(f"ROI == '{ROI_sel}' and band == '{band}' and cond in ['rsp_ctrl', 'oc_ctrl']")

            df_plot['cond'] = df_plot['cond'].replace({'rsp_ctrl' : 'RB', 'oc_ctrl' : 'O'})

            phase_order = ['inspi', 'expi']
            cond_order = ['RB', 'O']

            df_plot['phase_cycle'] = pd.Categorical(df_plot['phase_cycle'], categories=phase_order, ordered=True)
            df_plot['cond'] = pd.Categorical(df_plot['cond'], categories=cond_order, ordered=True)
            df_plot = df_plot.sort_values(['sujet', 'phase_cycle', 'cond'])

            fig, ax = plt.subplots(figsize=(7, 6))

            palette = {'RB': 'tab:blue', 'O': 'tab:orange'}

            sns.barplot(data=df_plot, x='phase_cycle', y='Pxx', hue='cond', order=phase_order, hue_order=cond_order,
                palette=palette, errorbar=None, alpha=0.4, dodge=True, ax=ax)

            phase_x = {phase: i for i, phase in enumerate(phase_order)}

            cond_offset = {'RB': -0.2, 'O': 0.2}

            for (sujet, phase), df_sub in df_plot.groupby(['sujet', 'phase_cycle'],observed=True):

                df_sub = df_sub.sort_values('cond')

                if df_sub['cond'].nunique() != len(cond_order):
                    continue

                x = [phase_x[phase] + cond_offset[cond] for cond in df_sub['cond']]

                y = df_sub['Pxx'].to_numpy()

                ax.plot(x, y, color='gray', alpha=0.45, linewidth=1, zorder=2)

            for cond in cond_order:

                df_cond = df_plot[df_plot['cond'] == cond]

                x = [phase_x[phase] + cond_offset[cond] for phase in df_cond['phase_cycle']]

                ax.scatter(x, df_cond['Pxx'], color=palette[cond], edgecolor='black', linewidth=0.4, s=45, zorder=3)

            ax.set_xlabel("Respiratory phase", fontsize=18)
            ax.set_ylabel("Pxx", fontsize=18)

            ax.tick_params(axis="x", labelsize=18)
            ax.tick_params(axis="y", labelsize=18)
            ax.legend(title=None)

            plt.suptitle(f"{ROI_sel} {band}")
            plt.tight_layout()

            # plt.show()

            fig.savefig(f"{ROI_sel}_{band}_2cond_patientlinked.png")

            _current_pair = [band, ROI_sel]

            if np.array([True for token in example_plot if _current_pair == token]).any():

                fig.savefig(os.path.join(path_paper_figure_export, f"fig04_{ROI_sel}_{band}_2cond_patientlinked.svg"))



def figure_oc_ratio_and_power():

    #### config
    path_load_tf = os.path.join(path_precompute, 'TF', 'session')
    path_load_resp_cycle = os.path.join(path_precompute, 'RESP', 'session')
    path_load_ocratio = os.path.join(path_precompute, 'RESP', 'respfeatures')

    #### get oc ratio extremum
    df_oc_ratio = pd.concat([pd.read_excel(os.path.join(path_load_ocratio, f"{sujet}_df_oc.xlsx")) for sujet in sujet_list]).drop(columns='Unnamed: 0')
    df_oc_ratio = df_oc_ratio.query(f"cond == 'oc_ctrl'")
    diff_oc_ratio = np.array([df_oc_ratio.query(f"sujet == '{sujet}'")['oc_ratio'].max() - df_oc_ratio.query(f"sujet == '{sujet}'")['oc_ratio'].min() for sujet in sujet_list])
    sujet_sel = sujet_list[np.argmax(diff_oc_ratio)]
    resp_cycle_i_list = [df_oc_ratio.query(f"sujet == '{sujet_sel}'")['oc_ratio'].argmin(), df_oc_ratio.query(f"sujet == '{sujet_sel}'")['oc_ratio'].argmax()]

    #### load
    cond_sel = 'oc_ctrl'
    tf = np.load(os.path.join(path_load_tf, f"{sujet_sel}_{cond_sel}_tf_allchan_stretch_post.npy"))
    resp_cycle = np.load(os.path.join(path_load_resp_cycle, f"{sujet_sel}_{cond_sel}_stretch_resp_post.npy"))
    df_loca = get_df_loca_allsujet().query(f"sujet == '{sujet_sel}'").drop(columns='Unnamed: 0').reset_index()

    #### explore cycles
    max_cycle = resp_cycle[resp_cycle_i_list[-1]]
    
    for cycle_i in range(resp_cycle.shape[0]):

        plt.plot(scipy.stats.zscore(max_cycle))
        plt.plot(scipy.stats.zscore(resp_cycle[cycle_i]))
        plt.title(cycle_i)
        plt.show()

    resp_cycle_i_list = [13, resp_cycle_i_list[-1]] #7,13,20 good cycle 

    #### chunk
    band_sel = 'theta'
    freq_sel_vec = (freq_band_dict[band_sel][0] < frex) & (frex < freq_band_dict[band_sel][1])
    tf_band = np.median(tf[:,:,freq_sel_vec], axis=2)

    resp_cycle_chunk = resp_cycle[resp_cycle_i_list]
    tf_band_chunk = scipy.signal.savgol_filter(tf_band[:,resp_cycle_i_list], window_length=20, polyorder=2, axis=-1)

    #### explore
    ROI_sujet = [_loca for _loca in df_loca['loca'].unique() if _loca in ROI_short_list]
    chan_i_sel = df_loca.query(f"loca in {ROI_sujet}").index

    for chan_i in chan_i_sel:

        loca = df_loca.iloc[chan_i]['loca']
    
        fig, axs = plt.subplots(ncols=2, figsize=(15,8))

        resp_cycle_min, resp_cycle_max = scipy.stats.zscore(resp_cycle_chunk, axis=1).min(), scipy.stats.zscore(resp_cycle_chunk, axis=1).max()
        tf_band_min, tf_band_max = scipy.stats.zscore(tf_band_chunk[chan_i], axis=1).min(), scipy.stats.zscore(tf_band_chunk[chan_i], axis=1).max()
        vmin, vmax = np.array([resp_cycle_min, tf_band_min]).min()*1.1, np.array([resp_cycle_max, tf_band_max]).max()*1.1

        data_plot = {'resp' : resp_cycle_chunk, band_sel : tf_band_chunk[chan_i]}

        for c, plot_type in enumerate(['resp', band_sel]):

            ax = axs[c]

            ax.plot(scipy.stats.zscore(data_plot[plot_type][0]))
            ax.plot(scipy.stats.zscore(data_plot[plot_type][1]))

            ax.set_title(f"chan{chan_i}:{loca}, {plot_type}")
            ax.set_ylim(vmin,vmax)
            ax.vlines(int(stretch_point_TF/2), ymin=vmin, ymax=vmax, color='k')
            ax.set_ylabel('zscore', fontsize=16)
            ax.set_xlabel('Phase', fontsize=16)
            ax.tick_params(axis="x", labelsize=16)
            ax.tick_params(axis="y", labelsize=16)

        plt.show()

    #### save
    chan_i_to_save = 98

    loca = df_loca.iloc[chan_i_to_save]['loca']
        
    fig, axs = plt.subplots(ncols=2, figsize=(15,6))

    resp_cycle_min, resp_cycle_max = scipy.stats.zscore(resp_cycle_chunk, axis=1).min(), scipy.stats.zscore(resp_cycle_chunk, axis=1).max()
    tf_band_min, tf_band_max = scipy.stats.zscore(tf_band_chunk[chan_i_to_save], axis=1).min(), scipy.stats.zscore(tf_band_chunk[chan_i_to_save], axis=1).max()
    vmin, vmax = np.array([resp_cycle_min, tf_band_min]).min()*1.1, np.array([resp_cycle_max, tf_band_max]).max()*1.1

    data_plot = {'resp' : resp_cycle_chunk, band_sel : tf_band_chunk[chan_i_to_save]}

    for c, plot_type in enumerate(['resp', band_sel]):

        ax = axs[c]

        ax.plot(scipy.stats.zscore(data_plot[plot_type][0]))
        ax.plot(scipy.stats.zscore(data_plot[plot_type][1]))

        ax.set_title(f"{sujet_sel}, chan{chan_i_to_save}:{loca}, {plot_type}")
        ax.set_ylim(vmin,vmax)
        ax.vlines(int(stretch_point_TF/2), ymin=vmin, ymax=vmax, color='k')
        ax.set_ylabel('zscore', fontsize=16)
        ax.set_xlabel('Phase', fontsize=16)
        ax.tick_params(axis="x", labelsize=16)
        ax.tick_params(axis="y", labelsize=16)

    # plt.show()

    fig.savefig(os.path.join(path_paper_figure_export, f"fig05a_occlusion_example.svg"))





def figure_oc_ratio_to_power_reg():

    #### config
    path_load_df = os.path.join(path_precompute, 'TF', 'session', 'df_R', 'df_reg_ALLROI_ALLDATA_R.xlsx')
    path_load_ocratio = os.path.join(path_precompute, 'RESP', 'respfeatures')

    #### get oc ratio extremum
    cond_sel = 'oc_ctrl'
    df_oc_ratio = pd.concat([pd.read_excel(os.path.join(path_load_ocratio, f"{sujet}_df_oc.xlsx")) for sujet in sujet_list]).drop(columns='Unnamed: 0')
    df_oc_ratio = df_oc_ratio.query(f"cond == '{cond_sel}'")
    diff_oc_ratio = np.array([df_oc_ratio.query(f"sujet == '{sujet}'")['oc_ratio'].max() - df_oc_ratio.query(f"sujet == '{sujet}'")['oc_ratio'].min() for sujet in sujet_list])
    sujet_sel = sujet_list[np.argmax(diff_oc_ratio)]

    #### load
    df_reg = pd.read_excel(path_load_df)
    df_reg_sujet = df_reg.query(f"sujet == '{sujet_sel}'")
    plot_ROI_list = [_loca for _loca in df_reg_sujet['ROI'].unique() if _loca in ROI_short_list]
    df_reg_sujet = df_reg_sujet.query(f"rf_metric == 'oc_ratio' and cond == '{cond_sel}' and ROI in {plot_ROI_list}")

    #### plot explore
    band_sel = 'theta'
    palette = {'inspi' : 'blue', 'expi' : 'red'}

    for ROI_sel in plot_ROI_list:

        df_plot = df_reg_sujet.query(f"sujet == '{sujet_sel}' and band == '{band_sel}' and ROI == '{ROI_sel}'")

        g = sns.lmplot(df_plot, x='rf_metric_val', y='Pxx', hue='phase_cycle', palette=palette)

        for ax in g.axes.flat:
            ax.tick_params(axis="x", labelsize=14)
            ax.tick_params(axis="y", labelsize=14)

            g.set_xlabels("oc_ratio", fontsize=16)
            g.set_ylabels("Pxx", fontsize=16)

        plt.suptitle(f"{sujet_sel}, {band_sel}, {ROI_sel}")

        plt.show()

    plt.close('all')

    #### plot save
    band_sel = 'theta'
    ROI_sel = 'medialorbitofrontal'
    palette = {'inspi' : 'blue', 'expi' : 'red'}

    df_plot = df_reg_sujet.query(f"sujet == '{sujet_sel}' and band == '{band_sel}' and ROI == '{ROI_sel}'")

    g = sns.lmplot(df_plot, x='rf_metric_val', y='Pxx', hue='phase_cycle', palette=palette, ci=False)

    for ax in g.axes.flat:
        ax.tick_params(axis="x", labelsize=14)
        ax.tick_params(axis="y", labelsize=14)

        g.set_xlabels("oc_ratio", fontsize=16)
        g.set_ylabels("Pxx", fontsize=16)

    plt.suptitle(f"{sujet_sel}, {band_sel}, {ROI_sel}")

    # plt.show()

    plt.savefig(os.path.join(path_paper_figure_export, f"fig05a_patient_example_ocratio_reg.svg"))

    plt.close('all')

    #### plot allsujet explore
    band_sel = 'beta'
    palette = {'inspi' : 'blue', 'expi' : 'red'}

    df_reg_allsujet = df_reg.query(f"cond in 'oc_ctrl' and rf_metric == 'oc_ratio'")

    plot_ROI_list = df_reg_allsujet.groupby(['ROI', 'sujet']).count().reset_index().groupby(['ROI']).count().reset_index()[['ROI', 'sujet']].query(f"sujet > 3")['ROI'].values

    for ROI_sel in plot_ROI_list:

        fig_allsujet, axs = plt.subplots(ncols=len(phase_cycle_list), figsize=(15, 7))
        
        for c, phase_cycle_sel in enumerate(phase_cycle_list):
    
            ax = axs[c]
    
            df_plot = df_reg_allsujet.query(f"band == '{band_sel}' and ROI == '{ROI_sel}' and phase_cycle == '{phase_cycle_sel}'")
    
            for sujet in df_plot["sujet"].unique():
    
                df_sub = df_plot[df_plot["sujet"] == sujet]
    
                sns.regplot(data=df_sub, x="rf_metric_val", y="Pxx", ci=None, ax=ax)
    
            ax.tick_params(axis="x", labelsize=14)
            ax.tick_params(axis="y", labelsize=14)
    
            ax.set_xlabel("oc_ratio", fontsize=16)
            ax.set_ylabel("Pxx", fontsize=16)
    
            ax.set_title(phase_cycle_sel, fontsize=16)
    
        plt.suptitle(f"{band_sel}, {ROI_sel}", fontsize=18)
    
        plt.tight_layout()
        plt.show()

    for band_sel in ['theta', 'beta', 'gamma']:

        for ROI_sel in plot_ROI_list:
                        
            df_plot = df_reg_allsujet.query(f"band == '{band_sel}' and ROI == '{ROI_sel}'")

            g = sns.lmplot(data=df_plot, x="rf_metric_val", y="Pxx", hue='phase_cycle', palette=palette, ci=None)
            
            for ax in g.axes.flat:
                ax.tick_params(axis="x", labelsize=14)
                ax.tick_params(axis="y", labelsize=14)

                g.set_xlabels("oc_ratio", fontsize=16)
                g.set_ylabels("Pxx", fontsize=16)

            plt.suptitle(f"{band_sel}, {ROI_sel}", fontsize=18)

            plt.tight_layout()
            plt.show()
        
    plt.close('all')


    #### plot allsujet save
    band_sel = 'beta'
    ROI_sel = 'precentral'
    palette = {'inspi' : 'blue', 'expi' : 'red'}

    df_plot = df_reg_allsujet.query(f"band == '{band_sel}' and ROI == '{ROI_sel}'")
    
    g = sns.lmplot(data=df_plot, x="rf_metric_val", y="Pxx", hue='phase_cycle', palette=palette, ci=None, legend=None)
    
    for ax in g.axes.flat:
        ax.tick_params(axis="x", labelsize=14)
        ax.tick_params(axis="y", labelsize=14)

        g.set_xlabels("oc_ratio", fontsize=16)
        g.set_ylabels("Pxx", fontsize=16)

    plt.suptitle(f"{band_sel}, {ROI_sel}", fontsize=18)

    plt.tight_layout()

    # plt.show()

    plt.savefig(os.path.join(path_paper_figure_export, f"fig05a_allpatient_example_ocratio_reg.svg"))

    plt.close('all')










########################################
######## EXPORT OC STATS ########
########################################

#sujet = sujet_list[0]
def export_plot_oc_stats():

    #### get data
    df_oc_allsujet = []
    
    for sujet in sujet_list:

        path_respfeatures = os.path.join(path_precompute, 'RESP', 'respfeatures')
        filename = f"{sujet}_df_oc.xlsx"
        _df_oc = pd.read_excel(os.path.join(path_respfeatures, filename))
        df_oc_allsujet.append(_df_oc)

    df_oc_allsujet = pd.concat(df_oc_allsujet).drop(columns=['Unnamed: 0'])
    df_median_oc = df_oc_allsujet.groupby(['sujet', 'cond']).median().reset_index()
    df_median_oc_patientlinked = df_oc_allsujet.groupby(['sujet', 'cond']).median('oc_ratio').query(f"cond in ['oc_ctrl', 'oc_chl']").reset_index().drop(columns=['cycle_i', 'oc_val'])

    #### sujetwise
    fig_allsujet_oc, ax = plt.subplots(figsize=(10,8))
    sns.boxplot(df_oc_allsujet, x='sujet', y='oc_ratio', hue='cond', showfliers=False, ax=ax)
    # plt.show()

    path_savefig = os.path.join(path_results, 'respi', 'oc_ratio', 'plot')
    filename_savefig = f"OC_sujetwise.png"
    fig_allsujet_oc.savefig(os.path.join(path_savefig, filename_savefig))

    #### median
    df_wide = df_median_oc.pivot(index='sujet', columns='cond', values='oc_ratio')

    # run test
    statfried, p = scipy.stats.friedmanchisquare(df_wide['oc_chl'], df_wide['oc_ctrl'], df_wide['rsp_chl'], df_wide['rsp_ctrl'])

    results = []

    for c1, c2 in itertools.combinations(conditions, 2):
        
        stat, p = scipy.stats.wilcoxon(df_wide[c1], df_wide[c2])
    
        results.append({"cond1": c1, "cond2": c2, "stat": stat, "p": p})

    posthoc = pd.DataFrame(results)
    posthoc["p_corrected"] = multipletests(posthoc["p"], method="bonferroni")[1]
    posthoc['friedman_p'] = [p] * posthoc.shape[0]
    posthoc['friedman_statval'] = [statfried] * posthoc.shape[0]

    # fig_median_oc, ax = plt.subplots(figsize=(10,8))
    # sns.swarmplot(df_median_oc, x='cond', y='oc_ratio', hue='sujet', ax=ax, size=10)
    # plt.show()

    fig_median_oc, ax = plt.subplots(figsize=(10,8))
    sns.swarmplot(df_median_oc.query(f"cond in ['oc_ctrl', 'oc_chl']"), x='cond', y='oc_ratio', hue='sujet', ax=ax, size=10, order=['oc_ctrl', 'oc_chl'])
    # plt.show()

    filename_savedf = f"df_stat_oc.xlsx"
    posthoc.to_excel(os.path.join(path_savefig, filename_savedf))

    #### patient linked
    df_plot = df_median_oc_patientlinked.copy()
    df_plot['cond'] = df_plot['cond'].replace({'oc_chl' : str('O+Ch'), 'oc_ctrl' : str('O')})

    fig_patientlinked, ax = plt.subplots(figsize=(6, 8))

    g = sns.barplot(
        data=df_plot,
        x="cond",
        y="oc_ratio",
        order=[str("O"), str("O+Ch")],
        estimator=np.median,
        errorbar=None,
        palette={str("O"): "steelblue", str("O+Ch"): "orange"},
        alpha=0.4,
        ax=ax,
        legend=False,
    )

    sns.pointplot(
        data=df_plot,
        x="cond",
        y="oc_ratio",
        hue="sujet",
        order=[str("O"), str("O+Ch")],
        estimator=np.median,
        errorbar=None,
        dodge=False,
        markers="o",
        markersize=8,
        linestyles="-",
        linewidth=2,
        alpha=0.7,
        palette=["gray"] * df_plot["sujet"].nunique(),
        ax=ax,
        legend=False,
    )

    # axis-label fontsize
    ax.xaxis.label.set_size(20)
    ax.yaxis.label.set_size(20)

    # tick-label fontsize
    ax.tick_params(axis="x", labelsize=20)
    ax.tick_params(axis="y", labelsize=20)

    # title
    ax.set_title("Occlusion ratio", fontsize=20)

    plt.tight_layout()

    # plt.show()

    path_savefig = os.path.join(path_results, 'respi', 'oc_ratio', 'plot')
    filename_savefig = f"psycho_patientlinked.png"
    fig_patientlinked.savefig(os.path.join(path_savefig, filename_savefig))

    filename_savefig = f"fig02d_OCratio_median.svg"
    fig_patientlinked.savefig(os.path.join(path_paper_figure_export, filename_savefig))

    plt.close('all')



################################
######## EXECUTE ########
################################


if __name__ == '__main__':

    export_challenge_stats()
    export_plot_oc_resp()
    export_plot_oc_stats()
    export_Pxx_diff_with_patientwise_info()

    figure_oc_ratio_and_power()
    figure_oc_ratio_to_power_reg()

                        