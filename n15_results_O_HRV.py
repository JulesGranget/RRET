



from n00_config_O_params import *
from n00bis_config_O_analysis_functions import *
from n00ter_X_manip_data import *

import plotly.express as px
import plotly.graph_objects as go








########################################
######## EXPORT ANAT PLOT ########
########################################



def export_fig_ecg():

    #### import data
    path_import_df_ecg = os.path.join(path_precompute, 'ECG', f"df_ecg.xlsx")
    df_ecg = pd.read_excel(path_import_df_ecg)

    path_import_xr_ecg = os.path.join(path_precompute, 'ECG', f"ecg_xr.nc")
    xr_ecg = xr.load_dataarray(path_import_xr_ecg)

    #### config
    HRV_list_short = ['HRV_MCV', 'HRV_Mad', 'HRV_Median', 'HRV_RMSSD']
    export_path_fig = os.path.join(path_results, 'ECG')

    #### sig plot    
    for phase in ["start", "end"]:

        fig = go.Figure()

        # time vector
        time = xr_ecg.coords["time"].values

        # RB median across patients
        rb_median = xr_ecg.sel(cond="RB", phase=phase).median(dim="patient").values

        fig.add_trace(go.Scatter(x=time, y=rb_median, mode="lines", name="RB median", line=dict(width=4, dash="dash")))

        # CHL vector for each patient
        for patient in xr_ecg.coords["patient"].values:

            chl_vec = xr_ecg.sel(patient=patient, cond="CHL", phase=phase).values

            fig.add_trace(go.Scatter(x=time, y=chl_vec, mode="lines", name=f"CHL {patient}"))

        fig.update_layout(template="simple_white", title=f"ECG - {phase}", xaxis_title="Time (s)", yaxis_title="ECG", width=900, height=500)
        # fig.show()

        fig.write_html(os.path.join(export_path_fig, f"iHR_time_seg_{phase}.html"))
        fig.write_image(os.path.join(export_path_fig, f"iHR_time_seg_{phase}.svg"))

    #### CHL
    df_plot = pd.melt(df_ecg, id_vars=['sujet', 'trial_num', 'cond', 'time_seg', 'HR_count', 'ecg_length', 'included', 'unpleasantness', 'anxiety'], value_vars=HRV_list_short)
    fig = sns.catplot(df_plot, kind='box', x='sujet', y='value', hue='cond', col='variable', row='time_seg', sharey=False, showfliers=False)
    # plt.show()

    fig.savefig(os.path.join(export_path_fig, f"CHL.png"))
    fig.savefig(os.path.join(export_path_fig, f"CHL.svg"))

    #### TIME DIFF
    df_plot = pd.melt(df_ecg, id_vars=['sujet', 'trial_num', 'cond', 'time_seg', 'HR_count', 'ecg_length', 'included', 'unpleasantness', 'anxiety'], value_vars=HRV_list_short)

    df_diff = df_plot.query(f"time_seg == 'start'").copy()
    df_diff[f'diff_time'] = df_plot.query(f"time_seg == 'end'")['value'].values - df_plot.query(f"time_seg == 'start'")['value'].values
    
    fig = sns.catplot(df_diff, kind='box', x='sujet', y='diff_time', hue='cond', col='variable', sharey=False, showfliers=False)
    for ax in fig.axes.flat:
        ax.axhline(0, color='k', linestyle='-', linewidth=2)
    plt.suptitle(f"end - start")
    # plt.show()

    fig.savefig(os.path.join(export_path_fig, f"TIME.png"))
    fig.savefig(os.path.join(export_path_fig, f"TIME.svg"))

    #### psycho
    psycho_metrics = ['unpleasantness', 'anxiety']
    df_psycho = pd.melt(df_diff, id_vars=[col for col in df_diff.columns if col not in psycho_metrics], value_vars=psycho_metrics, value_name='psycho_score', var_name='psycho_metric')
    
    for psycho_metric in psycho_metrics:
        fig = sns.catplot(df_psycho.query(f"psycho_metric == '{psycho_metric}'"), kind='swarm', x='psycho_score', y='diff_time', hue='sujet', col='variable', row='cond', sharey=False)
        plt.suptitle(f"{psycho_metric} \n end - start")
        # plt.show()

        fig.savefig(os.path.join(export_path_fig, f"PSYCHO_{psycho_metric}.png"))
        fig.savefig(os.path.join(export_path_fig, f"PSYCHO_{psycho_metric}.svg"))

    plt.close('all')

















    #### import df
    df_allpatient_plot = get_df_loca_allsujet_raw()

    #### export exclusion count
    unsorted_selected = []
    outside_our_ROI_selected = []

    for _, row_val in df_allpatient_plot.iterrows():

        if row_val['Select'] == 1:

            if row_val['loca'].find('UNSORTED') != -1:

                unsorted_selected.append(row_val['loca'])

            if row_val['loca'].find('UNSORTED') == -1 and row_val['loca'] not in ROI_short_list and row_val['loca'] not in ['WM', 'unknown']:

                outside_our_ROI_selected.append(row_val['loca'])

    df_export_count = pd.DataFrame({'total_contact' : [df_allpatient_plot.shape[0]], 'total_contact_removed' : [df_allpatient_plot.query(f"Select == 1").shape[0]], 
                    'WM' : [df_allpatient_plot.query(f"Select == 1 and loca == 'WM'").shape[0]],
                    'unknown' : [df_allpatient_plot.query(f"Select == 1 and loca == 'unknown'").shape[0]],
                    'UNSORTED' : np.array(unsorted_selected).shape[0],
                    'outside_our_ROI' : np.array(outside_our_ROI_selected).shape[0],
                   'total_contact_removed_filtered' : [df_allpatient_plot.query(f"Select == 1 and loca in {ROI_short_list}").shape[0]],
                   })
    
    filename_export = os.path.join(path_results, 'anatomy', f"df_exclusion_count.xlsx")
    df_export_count.to_excel(filename_export)

    #### export count subject
    df_loca_allsujet = get_df_loca_allsujet()
    
    df_sujet_count = df_loca_allsujet.query(f"loca in {ROI_short_list}")
    df_sujet_count = df_sujet_count.groupby("loca")["sujet"].nunique().reset_index(name="count")

    for sujet_n_thresh in np.arange(1,5):
        
        df_plot = df_sujet_count.query(f"count >= {sujet_n_thresh}")
        fig = px.bar(df_plot, x="loca", y="count", title=f"Sujet count thresh:{sujet_n_thresh}")
        # fig.show()

        os.chdir(os.path.join(path_results, 'anatomy'))
        fig.write_html(f"sujet_count_{sujet_n_thresh}thresh.html")

    #### export loca count
    df_loca_allsujet = get_df_loca_allsujet()
    
    df_loca_allsujet_filt = df_loca_allsujet.query(f"loca in {ROI_short_list}")
    df_count = df_loca_allsujet_filt.groupby(["sujet", "loca"]).count().reset_index()[['sujet', 'loca', 'chan']]

    df_ROI_scount = df_count.groupby('loca').nunique('sujet')['sujet'].reset_index(name='count')
    localist_thresh = df_ROI_scount.query(f"count >= {sujet_thresh}")['loca'].values.tolist()

    df_count = df_count.query(f"loca in {ROI_short_list}").query(f"loca in {localist_thresh}")

    df_plot = (
        df_count
        .pivot_table(index="loca", columns="sujet", values="chan", aggfunc="sum", fill_value=0)
    )

    df_plot["total"] = df_plot.sum(axis=1)
    df_plot = df_plot.sort_values("total", ascending=False).drop(columns="total")

    fig_anat, ax = plt.subplots(figsize=(8, 8))
    df_plot.plot(kind="bar", stacked=True, ax=ax)

    ax.set_xlabel("loca")
    ax.set_ylabel("Total chan count")
    ax.set_title("Total chan per loca, colored by sujet contribution")
    ax.legend(title="sujet", bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    # plt.show()

    path_savefig = os.path.join(path_results, 'anatomy')
    filename_savefig = f"allsujet_anat.png"
    fig_anat.savefig(os.path.join(path_savefig, filename_savefig))

    #### export df anat
    df_export = df_count.pivot_table(index="loca", columns="sujet", values="chan", aggfunc="sum", fill_value=0)
    df_export['total'] = df_export.values.sum(axis=1)

    output_df = os.path.join(path_results, "anatomy")
    filename = os.path.join(output_df, "df_anat_allsujet_count.xlsx")
    df_export.to_excel(filename)







################################
######## EXECUTE ########
################################


if __name__ == '__main__':

    export_fig_ecg()



                        