



from n00_config_O_params import *
from n00bis_config_O_analysis_functions import *
from n00ter_X_manip_data import *


import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go





################################
######## GENERATE ALLDF ########
################################





def p_to_stars(p):
    if p < 0.001:
        return '***'
    elif p < 0.01:
        return '**'
    elif p < 0.05:
        return '*'
    else:
        return ''













########################
######## OLSa ########
########################





def export_res_OLSa_reg_allcond():

    #### params
    sujet_thresh = 3
    pre_post_list = ['pre', 'post']
    rf_metric_short_list = ['cycle_freq', 'amplitude', 'oc_ratio', 'oc_val']

    #### load df
    os.chdir(os.path.join(path_precompute, 'TF', 'session', 'df_reg'))

    df_loca_allsujet = get_df_loca_allsujet()

    df_reg = pd.read_excel(f"df_reg_ALLROI.xlsx")
    
    label_name_nsujet = {}
    
    for ROI in ROI_short_list:

        _nsujet = df_loca_allsujet.query(f"loca == '{ROI}'")['sujet'].unique().size
        label_name_nsujet[ROI] = f"{ROI} s({_nsujet})"

    for ROI in ROI_short_list:
        df_reg['ROI'] = df_reg['ROI'].replace({ROI: label_name_nsujet[ROI]})

    df_ROI_scount = df_loca_allsujet.groupby('loca').nunique('sujet')['sujet'].reset_index(name='count')

    for ROI in ROI_short_list:
        df_ROI_scount['loca'] = df_ROI_scount['loca'].replace({ROI: label_name_nsujet[ROI]})

    localist_thresh = df_ROI_scount.query(f"count >= {sujet_thresh}")['loca'].values.tolist()

    ROI_plot_short_list = []

    for _ROI_thresh in localist_thresh:
        for _ROI_short in ROI_short_list:
            if _ROI_thresh.find(_ROI_short) != -1:
                ROI_plot_short_list.append(_ROI_thresh)

    df_reg_short = df_reg.query(f"ROI in {ROI_plot_short_list}")

    #### plot allcond
    for band in freq_band_dict:

        for phase_protocol in pre_post_list:

            for phase_cycle in phase_cycle_list:

                # --- Filter once ---
                df_plot = (
                    df_reg_short
                    .query(f"band == '{band}' and phase_protocol == '{phase_protocol}' and phase_cycle == '{phase_cycle}'")
                )

                rf_metrics = rf_metric_short_list
                conditions = list(df_plot["cond"].unique())

                fig = make_subplots(
                    rows=1,
                    cols=len(rf_metrics),
                    subplot_titles=rf_metrics,
                    shared_yaxes=False
                )

                for col_i, rf_metric in enumerate(rf_metrics, start=1):
                    df_sub = df_plot.query("rf_metric == @rf_metric")

                    # Keep a stable ROI order (you can change sorting if you want)
                    rois = list(df_sub["ROI"].dropna().unique())
                    x_base = np.arange(len(rois))

                    # bar geometry (controls grouping)
                    n_cond = len(conditions)
                    group_width = 0.8                      # total width reserved per ROI group
                    bar_w = group_width / n_cond           # each condition bar width
                    offsets = (np.arange(n_cond) - (n_cond - 1) / 2) * bar_w  # symmetric shifts

                    # dynamic vertical offset for star
                    y_max = df_sub["OLS_a"].abs().max()
                    y_off = (y_max * 0.05) if (y_max and np.isfinite(y_max)) else 0.05

                    for j, cond in enumerate(conditions):
                        df_c = df_sub.query("cond == @cond").set_index("ROI")

                        # y values aligned to the ROI list
                        y = np.array([df_c["OLS_a"].get(roi, np.nan) for roi in rois], dtype=float)
                        signi = np.array([df_c["rho_signi"].get(roi, 0) for roi in rois], dtype=float)

                        x = x_base + offsets[j]

                        # --- Bars (no outline) ---
                        fig.add_trace(
                            go.Bar(
                                x=x,
                                y=y,
                                width=bar_w * 0.95,
                                name=cond,
                                marker=dict(line=dict(width=0)),
                                showlegend=(col_i == 1)  # show legend only once
                            ),
                            row=1,
                            col=col_i
                        )

                        # --- Stars aligned to this condition's x shift ---
                        mask = (signi == 1) & np.isfinite(y)
                        if mask.any():
                            y_star = np.where(y[mask] >= 0, y[mask] + y_off, y[mask] - y_off)

                            fig.add_trace(
                                go.Scatter(
                                    x=x[mask],
                                    y=y_star,
                                    mode="text",
                                    text=["★"] * int(mask.sum()),
                                    textposition="middle center",
                                    textfont=dict(size=16),
                                    showlegend=False
                                ),
                                row=1,
                                col=col_i
                            )

                    # zero line
                    fig.add_hline(y=0, line_dash="dash", line_color="red", row=1, col=col_i)

                    # Set ROI labels on the x-axis of THIS subplot
                    axis_id = "" if col_i == 1 else str(col_i)
                    fig.update_xaxes(
                        tickmode="array",
                        tickvals=x_base,
                        ticktext=rois,
                        tickangle=45,
                        row=1,
                        col=col_i
                    )

                fig.update_layout(
                    barmode="overlay",   # we're manually positioning bars, so overlay is correct
                    template="simple_white",
                    height=520,
                    width=420 * len(rf_metrics),
                    title=f"{phase_protocol} {phase_cycle} {band}",
                    bargap=0.15
                )

                # fig.show()

                filepath = os.path.join(path_results, 'Pxx', 'reg_with_RF', 'allcond', f"REG_RF_{band}_{phase_cycle}_{phase_protocol}.html")
                fig.write_html(filepath)





########################
######## LMM CTRL ########
########################





def export_reg_CTRL_example():

    #### load
    filename = os.path.join(path_precompute, 'TF', 'session', 'df_R', f"df_reg_ALLROI_ALLDATA_R.xlsx")
    df_reg_allROI_alldata_R = pd.read_excel(filename)
    band_list = ["theta", "beta", "gamma"]
    sujet_thresh = 3

    df_loca_allsujet = get_df_loca_allsujet()

    label_name_nsujet = {}
    
    for ROI in ROI_short_list:

        _nsujet = df_loca_allsujet.query(f"loca == '{ROI}'")['sujet'].unique().size
        _ncontact = df_loca_allsujet.query(f"loca == '{ROI}'").groupby("sujet").count()['chan'].sum()        

        label_name_nsujet[ROI] = f"{ROI} s({_nsujet}) c({_ncontact})"

    for ROI in ROI_short_list:
        df_reg_allROI_alldata_R['ROI'] = df_reg_allROI_alldata_R['ROI'].replace({ROI: label_name_nsujet[ROI]})

    df_ROI_scount = df_loca_allsujet.groupby('loca').nunique('sujet')['sujet'].reset_index(name='count')

    for ROI in ROI_short_list:
        df_ROI_scount['loca'] = df_ROI_scount['loca'].replace({ROI: label_name_nsujet[ROI]})

    localist_thresh = df_ROI_scount.query(f"count >= {sujet_thresh}")['loca'].values.tolist()

    ROI_plot_short_list = []

    for _ROI_thresh in localist_thresh:
        for _ROI_short in ROI_short_list:
            if _ROI_thresh.find(_ROI_short) != -1:
                ROI_plot_short_list.append(_ROI_thresh)

    #### plot

    ROI = ROI_plot_short_list[0]
    band = "theta"
    
    for ROI in ROI_plot_short_list:

        for band in band_list:
    
            df_plot = df_reg_allROI_alldata_R.query(f"ROI == '{ROI}' and band == '{band}' and cond == 'rsp_ctrl'")

            sns.lmplot(data=df_plot, x="rf_metric_val", y="Pxx", hue="sujet", col="phase_cycle", row="rf_metric", sharex=False, sharey=False, ci=None,  # remove confidence band
                scatter_kws={"alpha": 0.2},  # make points transparent
                line_kws={"linewidth": 3}  # optional: thicker regression line
            )

            plt.subplots_adjust(top=0.92)
            plt.suptitle(f"{ROI} {band}")

            # plt.show()

            filename = os.path.join(path_results, "Pxx", "reg_with_RF", "ctrl", "examples", f"lmplot_{ROI}_{band}.png")
            plt.savefig(filename)

            plt.close('all')



def export_res_LMM_reg_CTRL():

    #### generate df
    #### params
    sujet_thresh = 3
    pre_post_list = ['pre', 'post']
    rf_metric_short_list = ['amplitude', 'oc_ratio', 'oc_val']
    band_list_short = ['theta', 'beta', 'gamma']

    #### load df
    folder_import_dir = os.path.join(path_results, 'LMM', 'df')

    df_loca_allsujet = get_df_loca_allsujet()

    df_LMM_CTRL = []

    for phase_cycle_sel in phase_cycle_list:
        
        for rf_metric_sel in rf_metric_short_list:
            
            for band in band_list_short:
            
                for ROI in ROI_short_list:
    
                    _filename = os.path.join(folder_import_dir, f"CTRL_RES_{band}_post_{phase_cycle_sel}_{rf_metric_sel}_{ROI}_LMM.xlsx")
                    df_LMM_CTRL.append(pd.read_excel(_filename))

    df_LMM_CTRL = pd.concat(df_LMM_CTRL)
    
    label_name_nsujet = {}
    
    for ROI in ROI_short_list:

        _nsujet = df_loca_allsujet.query(f"loca == '{ROI}'")['sujet'].unique().size
        _ncontact = df_loca_allsujet.query(f"loca == '{ROI}'").groupby("sujet").count()['chan'].sum()        

        label_name_nsujet[ROI] = f"{ROI} s({_nsujet}) c({_ncontact})"

    for ROI in ROI_short_list:
        df_LMM_CTRL['ROI'] = df_LMM_CTRL['ROI'].replace({ROI: label_name_nsujet[ROI]})

    df_ROI_scount = df_loca_allsujet.groupby('loca').nunique('sujet')['sujet'].reset_index(name='count')

    for ROI in ROI_short_list:
        df_ROI_scount['loca'] = df_ROI_scount['loca'].replace({ROI: label_name_nsujet[ROI]})

    localist_thresh = df_ROI_scount.query(f"count >= {sujet_thresh}")['loca'].values.tolist()

    ROI_plot_short_list = []

    for _ROI_thresh in localist_thresh:
        for _ROI_short in ROI_short_list:
            if _ROI_thresh.find(_ROI_short) != -1:
                ROI_plot_short_list.append(_ROI_thresh)

    df_LMM_CTRL_short = df_LMM_CTRL.query(f"ROI in {ROI_plot_short_list}")


    #### plot ctrl
    phase_cycle_list = ["inspi", "expi"]
    phase_protocol = "post"
    val_signi = 0.05

    for band in freq_band_dict:

        df_plot = df_LMM_CTRL_short.query("band == @band")

        fig = make_subplots(
            rows=1,
            cols=len(rf_metric_short_list),
            subplot_titles=rf_metric_short_list,
            shared_yaxes=False
        )

        for col_i, rf_metric in enumerate(rf_metric_short_list, start=1):
            df_sub = df_plot.query("rf_metric == @rf_metric")

            # stable ROI order
            rois = list(df_sub["ROI"].dropna().unique())
            x_base = np.arange(len(rois))

            # group geometry (2 hue levels = inspi/expi)
            n_phase = len(phase_cycle_list)
            group_width = 0.8
            bar_w = group_width / n_phase
            offsets = (np.arange(n_phase) - (n_phase - 1) / 2) * bar_w

            # y offset for stars (per subplot)
            y_max = df_sub["estimate"].abs().max()
            y_off = (y_max * 0.05) if (y_max and np.isfinite(y_max)) else 0.05

            for j, phase_cycle in enumerate(phase_cycle_list):
                # IMPORTANT FIX:
                # ensure ONE value per ROI (otherwise y becomes 2D and breaks masking)
                df_c = (
                    df_sub.query("phase_cycle == @phase_cycle")
                        .groupby("ROI", as_index=True)[["estimate", "p.value"]]
                        .mean()
                )

                y = np.array([df_c["estimate"].get(roi, np.nan) for roi in rois], dtype=float)
                p = np.array([df_c["p.value"].get(roi, np.nan) for roi in rois], dtype=float)
                signi = (p < val_signi)

                x = x_base + offsets[j]

                # Bars (no outline)
                fig.add_trace(
                    go.Bar(
                        x=x,
                        y=y,
                        width=bar_w * 0.95,
                        name=phase_cycle,
                        marker=dict(line=dict(width=0)),
                        showlegend=(col_i == 1)
                    ),
                    row=1,
                    col=col_i
                )

                # Stars aligned to the correct hue bar
                mask = signi & np.isfinite(y)
                if mask.any():
                    y_star = np.where(y[mask] >= 0, y[mask] + y_off, y[mask] - y_off)

                    fig.add_trace(
                        go.Scatter(
                            x=x[mask],
                            y=y_star,
                            mode="text",
                            text=["★"] * int(mask.sum()),
                            textposition="middle center",
                            textfont=dict(size=16),
                            showlegend=False
                        ),
                        row=1,
                        col=col_i
                    )

            # zero reference line
            fig.add_hline(y=0, line_dash="dash", line_color="red", row=1, col=col_i)

            # ROI tick labels centered on ROI groups (x_base)
            fig.update_xaxes(
                tickmode="array",
                tickvals=x_base,
                ticktext=rois,
                tickangle=45,
                row=1,
                col=col_i
            )

        fig.update_layout(
            barmode="overlay",  # we manually position bars on x, so overlay is correct
            template="simple_white",
            height=520,
            width=420 * len(rf_metric_short_list),
            title=f"{phase_protocol} {band}",
            bargap=0.15
        )

        fig.show()


        filepath = os.path.join(path_results, 'Pxx', 'reg_with_RF', 'ctrl', f"REG_RF_{band}_{phase_protocol}.html")
        fig.write_html(filepath)










########################
######## LMM OC ########
########################





def export_reg_OC_example():

    #### load
    filename = os.path.join(path_precompute, 'TF', 'session', 'df_R', f"df_reg_ALLROI_ALLDATA_R.xlsx")
    df_reg_allROI_alldata_R = pd.read_excel(filename)
    band_list = ["theta", "beta", "gamma"]
    cond_list_sel = ['oc_ctrl', 'oc_chl']
    sujet_thresh = 3

    df_loca_allsujet = get_df_loca_allsujet()

    label_name_nsujet = {}
    
    for ROI in ROI_short_list:

        _nsujet = df_loca_allsujet.query(f"loca == '{ROI}'")['sujet'].unique().size
        _ncontact = df_loca_allsujet.query(f"loca == '{ROI}'").groupby("sujet").count()['chan'].sum()        

        label_name_nsujet[ROI] = f"{ROI} s({_nsujet}) c({_ncontact})"

    for ROI in ROI_short_list:
        df_reg_allROI_alldata_R['ROI'] = df_reg_allROI_alldata_R['ROI'].replace({ROI: label_name_nsujet[ROI]})

    df_ROI_scount = df_loca_allsujet.groupby('loca').nunique('sujet')['sujet'].reset_index(name='count')

    for ROI in ROI_short_list:
        df_ROI_scount['loca'] = df_ROI_scount['loca'].replace({ROI: label_name_nsujet[ROI]})

    localist_thresh = df_ROI_scount.query(f"count >= {sujet_thresh}")['loca'].values.tolist()

    ROI_plot_short_list = []

    for _ROI_thresh in localist_thresh:
        for _ROI_short in ROI_short_list:
            if _ROI_thresh.find(_ROI_short) != -1:
                ROI_plot_short_list.append(_ROI_thresh)

    df_reg_allROI_alldata_R = df_reg_allROI_alldata_R.query(f"cond in {cond_list_sel}")

    #### plot
    ROI = ROI_plot_short_list[0]
    band = "theta"
    
    for ROI in ROI_plot_short_list:

        for band in band_list:

            for rf_metric in df_reg_allROI_alldata_R['rf_metric'].unique():
    
                df_plot = df_reg_allROI_alldata_R.query(f"ROI == '{ROI}' and band == '{band}' and rf_metric == '{rf_metric}'")

                sns.lmplot(data=df_plot, x="rf_metric_val", y="Pxx", hue="sujet", col="phase_cycle", row="cond", sharex=False, sharey=False, ci=None,  # remove confidence band
                    scatter_kws={"alpha": 0.2},  # make points transparent
                    line_kws={"linewidth": 3}  # optional: thicker regression line
                )

                plt.subplots_adjust(top=0.92)
                plt.suptitle(f"{ROI} {band} {rf_metric}")

                # plt.show()

                filename = f"lmplot_OC_{ROI}_{band}_{rf_metric}.png"
                filepath = os.path.join(path_results, "Pxx", "reg_with_RF", "oc", "examples", filename)
                plt.savefig(filepath)

                plt.close('all')







def export_res_LMM_reg_OC_ALLCOND():

    #### generate df
    #### params
    sujet_thresh = 3
    rf_metric_short_list = ['amplitude', 'oc_ratio', 'oc_val']
    band_list_short = ['theta', 'beta', 'gamma']
    phase_cycle_list = ["inspi", "expi"]
    phase_protocol = "post"
    val_signi = 0.05
    estimate_list = ["rf_metric_val", "condoc_chl", "rf_metric_val:condoc_chl"]

    #### load df
    folder_import_dir = os.path.join(path_results, 'LMM', 'df')

    df_loca_allsujet = get_df_loca_allsujet()

    df_LMM_OC = []

    for phase_cycle_sel in phase_cycle_list:
        
        for rf_metric_sel in rf_metric_short_list:
            
            for band in band_list_short:
            
                for ROI in ROI_short_list:
    
                    _filename = os.path.join(folder_import_dir, f"OC_ALLCOND_RES_{band}_post_{phase_cycle_sel}_{rf_metric_sel}_{ROI}_LMM.xlsx")
                    df_LMM_OC.append(pd.read_excel(_filename))

    df_LMM_OC = pd.concat(df_LMM_OC)
    
    label_name_nsujet = {}
    
    for ROI in ROI_short_list:

        _nsujet = df_loca_allsujet.query(f"loca == '{ROI}'")['sujet'].unique().size
        _ncontact = df_loca_allsujet.query(f"loca == '{ROI}'").groupby("sujet").count()['chan'].sum()        

        label_name_nsujet[ROI] = f"{ROI} s({_nsujet}) c({_ncontact})"

    for ROI in ROI_short_list:
        df_LMM_OC['ROI'] = df_LMM_OC['ROI'].replace({ROI: label_name_nsujet[ROI]})

    df_ROI_scount = df_loca_allsujet.groupby('loca').nunique('sujet')['sujet'].reset_index(name='count')

    for ROI in ROI_short_list:
        df_ROI_scount['loca'] = df_ROI_scount['loca'].replace({ROI: label_name_nsujet[ROI]})

    localist_thresh = df_ROI_scount.query(f"count >= {sujet_thresh}")['loca'].values.tolist()

    ROI_plot_short_list = []

    for _ROI_thresh in localist_thresh:
        for _ROI_short in ROI_short_list:
            if _ROI_thresh.find(_ROI_short) != -1:
                ROI_plot_short_list.append(_ROI_thresh)

    df_LMM_OC_short = df_LMM_OC.query(f"ROI in {ROI_plot_short_list}")


    #### plot 
    color_map = {
        "inspi": "#1f77b4",
        "expi":  "#d62728" 
    }

    ROI_order = ['Amygdala s(8) c(31)', 'Hippocampus s(7) c(57)', 'insula-ant s(5) c(21)', 'insula-pos s(4) c(18)', 
                 'lateralorbitofrontal s(5) c(21)', 'medialorbitofrontal s(5) c(7)', 'postcentral s(3) c(17)', 'precentral s(4) c(12)']

    for band in band_list_short:

        df_plot = df_LMM_OC_short.query("band == @band")

        fig = make_subplots(
            rows=len(estimate_list),
            cols=len(rf_metric_short_list),
            subplot_titles=rf_metric_short_list,
            shared_yaxes=False
        )

        for col_i, rf_metric in enumerate(rf_metric_short_list, start=1):

            y_max = df_plot.query("rf_metric == @rf_metric and term != '(Intercept)'")["estimate"].abs().max()

            for row_i, estimate in enumerate(estimate_list, start=1):

                df_sub = df_plot.query("rf_metric == @rf_metric and term == @estimate")
                
                for j, phase_cycle in enumerate(phase_cycle_list):

                    df_c = df_sub.query("phase_cycle == @phase_cycle")

                    signi_text = [p_to_stars(_p) for _p in df_c['p.value']]

                    fig.add_bar(
                        x=df_c['ROI'],
                        y=df_c['estimate'],
                        name=phase_cycle,
                        text=signi_text,
                        textposition="outside",
                        marker_color=color_map[phase_cycle],
                        row=row_i,
                        col=col_i
                    )

                    margin = y_max * 0.50   # 15% extra space

                    fig.update_yaxes(
                        range=[-y_max - margin, y_max + margin],
                        row=row_i,
                        col=col_i
                    )

                    fig.update_xaxes(
                        categoryorder="array",
                        categoryarray=ROI_order,
                        row=row_i,
                        col=col_i
                    )

                # ROI tick labels centered on ROI groups (x_base)
                if row_i == 3:
                    fig.update_xaxes(
                        tickmode="array",
                        # tickvals=x,
                        ticktext=ROI_order,
                        tickangle=45,
                        row=row_i,
                        col=col_i
                    )

                if col_i == 1:
                    fig.update_yaxes(
                        title_text=estimate,
                        row=row_i,
                        col=col_i
                    )

        fig.update_layout(
            barmode="group",  # we manually position bars on x, so overlay is correct
            template="simple_white",
            height=320 * len(estimate_list),
            width=420 * len(rf_metric_short_list),
            title=f"{phase_protocol} {band}",
            bargap=0.15
        )

        # fig.show()

        filepath = os.path.join(path_results, 'Pxx', 'reg_with_RF', 'oc', f"ALLCOND_REG_RF_{band}_{phase_protocol}.html")
        fig.write_html(filepath)

    #### export df values
    df_export = df_LMM_OC_short.query(f"ROI in {ROI_plot_short_list} and term != '(Intercept)' and rf_metric == 'oc_ratio'")
    
    filepath = os.path.join(path_results, 'Pxx', 'reg_with_RF', 'oc')
    filename = os.path.join(filepath, "ALLCOND_df_REG_oc_ratio.xlsx")
    df_export.to_excel(filename)







def export_res_LMM_reg_OC_UNIQUECOND():

    #### generate df
    #### params
    sujet_thresh = 3
    rf_metric_short_list = ['amplitude', 'oc_ratio', 'oc_val']
    band_list_short = ['theta', 'beta', 'gamma']
    phase_cycle_list = ["inspi", "expi"]
    cond_list_oc = ['oc_ctrl', 'oc_chl']
    phase_protocol = "post"
    val_signi = 0.05

    #### load df
    folder_import_dir = os.path.join(path_results, 'LMM', 'df')

    df_loca_allsujet = get_df_loca_allsujet()

    df_LMM_OC = []

    for cond in cond_list_oc:

        for phase_cycle_sel in phase_cycle_list:
            
            for rf_metric_sel in rf_metric_short_list:
                
                for band in band_list_short:
                
                    for ROI in ROI_short_list:
        
                        _filename = os.path.join(folder_import_dir, f"OC_UNIQUECOND_RES_{cond}_{band}_post_{phase_cycle_sel}_{rf_metric_sel}_{ROI}_LMM.xlsx")
                        df_LMM_OC.append(pd.read_excel(_filename))

    df_LMM_OC = pd.concat(df_LMM_OC)

    label_name_nsujet = {}
    
    for ROI in ROI_short_list:

        _nsujet = df_loca_allsujet.query(f"loca == '{ROI}'")['sujet'].unique().size
        _ncontact = df_loca_allsujet.query(f"loca == '{ROI}'").groupby("sujet").count()['chan'].sum()        

        label_name_nsujet[ROI] = f"{ROI} s({_nsujet}) c({_ncontact})"

    for ROI in ROI_short_list:
        df_LMM_OC['ROI'] = df_LMM_OC['ROI'].replace({ROI: label_name_nsujet[ROI]})

    df_ROI_scount = df_loca_allsujet.groupby('loca').nunique('sujet')['sujet'].reset_index(name='count')

    for ROI in ROI_short_list:
        df_ROI_scount['loca'] = df_ROI_scount['loca'].replace({ROI: label_name_nsujet[ROI]})

    localist_thresh = df_ROI_scount.query(f"count >= {sujet_thresh}")['loca'].values.tolist()

    ROI_plot_short_list = []

    for _ROI_thresh in localist_thresh:
        for _ROI_short in ROI_short_list:
            if _ROI_thresh.find(_ROI_short) != -1:
                ROI_plot_short_list.append(_ROI_thresh)

    df_LMM_OC_short = df_LMM_OC.query(f"ROI in {ROI_plot_short_list}")

    #### plot
    color_map = {
        "inspi": "#1f77b4",
        "expi":  "#d62728" 
    }

    ROI_order = ['Amygdala s(8) c(31)', 'Hippocampus s(7) c(57)', 'insula-ant s(5) c(21)', 'insula-pos s(4) c(18)', 
                 'lateralorbitofrontal s(5) c(21)', 'medialorbitofrontal s(5) c(7)', 'postcentral s(3) c(17)', 'precentral s(4) c(12)']

    for band in band_list_short:

        df_plot = df_LMM_OC_short.query("band == @band")

        fig = make_subplots(
            rows=len(cond_list_oc),
            cols=len(rf_metric_short_list),
            subplot_titles=rf_metric_short_list*2,
            shared_yaxes=False
        )

        for col_i, rf_metric in enumerate(rf_metric_short_list, start=1):

            y_max = df_plot.query("rf_metric == @rf_metric and term != '(Intercept)'")["estimate"].abs().max()

            for row_i, cond in enumerate(cond_list_oc, start=1):
                
                df_sub = df_plot.query("rf_metric == @rf_metric and cond == @cond and term == 'rf_metric_val'")

                for j, phase_cycle in enumerate(phase_cycle_list):

                    df_c = df_sub.query("phase_cycle == @phase_cycle")

                    signi_text = [p_to_stars(_p) for _p in df_c['p.value']]

                    fig.add_bar(
                        x=df_c['ROI'],
                        y=df_c['estimate'],
                        name=phase_cycle,
                        text=signi_text,
                        textposition="outside",
                        marker_color=color_map[phase_cycle],
                        row=row_i,
                        col=col_i
                    )

                    margin = y_max * 0.50   # 15% extra space

                    fig.update_yaxes(
                        range=[-y_max - margin, y_max + margin],
                        row=row_i,
                        col=col_i
                    )

                    fig.update_xaxes(
                        categoryorder="array",
                        categoryarray=ROI_order,
                        row=row_i,
                        col=col_i
                    )

                # zero reference line
                # fig.add_hline(y=0, line_dash="dash", line_color="red", row=row_i, col=col_i)

                # if col_i == 1:  
                #     fig.update_yaxes(
                #         title_text=cond,
                #         row=row_i,
                #         col=col_i
                #     )

                # margin = y_max * 0.50   # 15% extra space

                # fig.update_yaxes(
                #     range=[-y_max - margin, y_max + margin],
                #     row=row_i,
                #     col=col_i
                # )


        fig.update_layout(
            barmode="group",  # we manually position bars on x, so overlay is correct
            template="simple_white",
            xaxis_tickangle=-45,
            height=420 * len(cond_list_oc),
            width=420 * len(rf_metric_short_list),
            title=f"{phase_protocol} {band}",
            bargap=0.15
        )

        # fig.show()

        filepath = os.path.join(path_results, 'Pxx', 'reg_with_RF', 'oc', f"UNIQUECOND_REG_RF_{band}_{phase_protocol}.html")
        fig.write_html(filepath)

    #### export df values
    df_export = df_LMM_OC_short.query(f"ROI in {ROI_plot_short_list} and term != '(Intercept)' and rf_metric == 'oc_ratio'")
    
    filepath = os.path.join(path_results, 'Pxx', 'reg_with_RF', 'oc')
    filename = os.path.join(filepath, "UNIQUECOND_df_REG_oc_ratio.xlsx")
    df_export.to_excel(filename)






################################
######## EXECUTE ########
################################


if __name__ == '__main__':


    export_res_OLSa_reg_allcond()
    export_res_LMM_reg_CTRL()

    export_res_LMM_reg_OC_ALLCOND()
    export_res_LMM_reg_OC_UNIQUECOND()

    export_reg_OC_example()






