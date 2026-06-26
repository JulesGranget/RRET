



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

    #### load
    os.chdir(os.path.join(path_precompute, 'RESP', 'respfeatures')) 
    respfeatures = pd.read_excel(f'{sujet}_respfeatures_cleaned_label.xlsx')

    #### extract
    cond = 'oc_ctrl'

    os.chdir(os.path.join(path_precompute, 'RESP', 'session'))
    resp_stretch = np.load(f"{sujet}_{cond}_stretch_resp_post.npy")
    _respfeature = respfeatures.query(f"cond == '{cond}'")
    
    oc_ratio = []
    oc_val = []

    for cycle_i in range(resp_stretch.shape[0]): 

        sig = resp_stretch[cycle_i]
        sig_diff = np.diff(sig)
        start_inspi_dec = np.argmin(sig_diff[:int(stretch_point_TF/4)])
        start_inspi_asc = np.argmax(sig_diff[int(stretch_point_TF/4):int(stretch_point_TF/2)]) + int(stretch_point_TF/4)

        oc_top_i = np.where(sig_diff[start_inspi_dec:] > 0)[0][0] + start_inspi_dec
        oc_top_val = sig[oc_top_i]
        median_inspi_trough = np.median(sig[start_inspi_dec:start_inspi_asc])

        _oc_ratio = (oc_top_val - median_inspi_trough) / median_inspi_trough
        oc_ratio.append(_oc_ratio)
        oc_val.append(oc_top_val)

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

        for cycle_i in range(resp_stretch.shape[0]): 

            sig = resp_stretch[cycle_i][:-1]
            sig_diff = np.diff(resp_stretch[cycle_i])
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

    #### plot oc example
    cycle_i =0

    sig = resp_stretch[cycle_i][:-1]
    sig_diff = np.diff(resp_stretch[cycle_i])
    start_inspi_dec = np.argmin(sig_diff[:int(stretch_point_TF/4)])
    start_inspi_asc = np.argmax(sig_diff[int(stretch_point_TF/4):int(stretch_point_TF/2)]) + int(stretch_point_TF/4)

    oc_top_i = np.where(sig_diff[start_inspi_dec:] > 0)[0][0] + start_inspi_dec

    median_inspi_trough = np.median(sig[start_inspi_dec:start_inspi_asc])

    time_vec = np.arange(sig.size)

    fig_oc_example, ax = plt.subplots()
    ax.plot(time_vec, sig)
    ax.hlines([median_inspi_trough], xmin=0, xmax=stretch_point_TF, colors='b')
    # ax.hlines([0], xmin=0, xmax=stretch_point_TF, colors='r')
    ax.vlines([stretch_point_TF/2], ymin=sig.min(), ymax=sig.max(), colors='r')
    ax.scatter([oc_top_i], sig[oc_top_i], color='r')
    ax.set_title(f"{sujet} {cond} {cycle_i}")
    # plt.show()

    path_savefig = os.path.join(path_results, 'respi', 'oc_ratio', 'plot')
    filename_savefig = f"OC_example.png"
    fig_oc_example.savefig(os.path.join(path_savefig, filename_savefig))

    filename_savefig = f"fig03d_OC_example.svg"
    fig_oc_example.savefig(os.path.join(path_paper_figure_export, filename_savefig))






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

    fig_median_oc, ax = plt.subplots(figsize=(10,8))
    sns.swarmplot(df_median_oc, x='cond', y='oc_ratio', hue='sujet', ax=ax, size=10)
    # plt.show()

    fig_median_oc, ax = plt.subplots(figsize=(10,8))

    sns.swarmplot(
        data=df_median_oc,
        x='cond',
        y='oc_ratio',
        hue='sujet',
        size=10,
        ax=ax
    )

    # plt.show()

    path_savefig = os.path.join(path_results, 'respi', 'oc_ratio', 'plot')
    filename_savefig = f"OC_median.png"
    fig_median_oc.savefig(os.path.join(path_savefig, filename_savefig))

    filename_savefig = f"fig03e_OC_median.svg"
    fig_median_oc.savefig(os.path.join(path_paper_figure_export, filename_savefig))

    filename_savedf = f"df_stat_oc.xlsx"
    posthoc.to_excel(os.path.join(path_savefig, filename_savedf))




################################
######## EXECUTE ########
################################


if __name__ == '__main__':

    export_challenge_stats()
    export_plot_oc_resp()
    export_plot_oc_stats()

                        