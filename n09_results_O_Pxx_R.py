



from n00_config_O_params import *
from n00bis_config_O_analysis_functions import *
from n00ter_X_manip_data import *


import plotly.express as px
from plotly.subplots import make_subplots









################################
######## GENERATE ALLDF ########
################################

def generate_all_df_Pxx_wholecycle():

    os.chdir(os.path.join(path_results, 'LMM', 'df'))

    df_loca_allsujet = get_df_loca_allsujet()
    ROI_list = [_loca for _loca in df_loca_allsujet['loca'].unique() if _loca.find('UNSORTED') == -1]

    df_R_Pxx = []

    for phase_cycle in phase_cycle_list:

        for band in freq_band_dict:

            for ROI in ROI_list:

                _df = pd.read_excel(f"{band}_Pxx_lmm_{ROI}_{phase_cycle}_res.xlsx")
                _df_add = pd.concat([pd.DataFrame({'phase_cycle' : [phase_cycle]*_df.shape[0], 'band' : [band]*_df.shape[0], 'ROI' : [ROI]*_df.shape[0]}), _df[['term', 'estimate', 'p.value']]], axis=1)
                
                df_R_Pxx.append(_df_add)

    df_R_Pxx = pd.concat(df_R_Pxx)
    df_R_Pxx = df_R_Pxx.rename(columns={'p.value' : 'pvalue'})

    return df_R_Pxx







def generate_all_df_OLSa():

    os.chdir(os.path.join(path_results, 'LMM', 'df'))

    df_loca_allsujet = get_df_loca_allsujet()
    ROI_list = [_loca for _loca in df_loca_allsujet['loca'].unique() if _loca.find('UNSORTED') == -1]

    df_R_OLSa = []

    pre_post_list = ['pre', 'post']

    for phase_cycle in phase_cycle_list:

        for pre_post in pre_post_list:

            for rf_metric in rf_metric_list:

                for band in freq_band_dict:

                    for ROI in ROI_list:

                        _df = pd.read_excel(f"{band}_{rf_metric}_OLS_a_lmm_{ROI}_{phase_cycle}_{pre_post}_res.xlsx")
                        _df_add = pd.concat([pd.DataFrame({'phase_cycle' : [phase_cycle]*_df.shape[0], 'pre_post' : [pre_post]*_df.shape[0], 'band' : [band]*_df.shape[0],
                                                           'rf_metric' : [rf_metric]*_df.shape[0], 'ROI' : [ROI]*_df.shape[0]}), _df[['term', 'estimate', 'p.value']]], axis=1)
                        
                        df_R_OLSa.append(_df_add)

    df_R_OLSa = pd.concat(df_R_OLSa)
    df_R_OLSa = df_R_OLSa.rename(columns={'p.value' : 'pvalue'})

    return df_R_OLSa





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
######## PXX ########
########################




def export_res_Pxx():

    #### params
    sujet_thresh = 3

    #### load df
    df_R_Pxx = generate_all_df_Pxx_wholecycle()

    df_loca_allsujet = get_df_loca_allsujet()
    ROI_list = [_loca for _loca in df_loca_allsujet['loca'].unique() if _loca.find('UNSORTED') == -1]

    label_name_nsujet = {}
    
    for ROI in ROI_list:

        _nsujet = df_loca_allsujet.query(f"loca == '{ROI}'")['sujet'].unique().size
        label_name_nsujet[ROI] = f"{ROI} s({_nsujet})"

    for ROI in ROI_list:
        df_R_Pxx['ROI'] = df_R_Pxx['ROI'].replace({ROI: label_name_nsujet[ROI]})

    df_ROI_scount = df_loca_allsujet.groupby('loca').nunique('sujet')['sujet'].reset_index(name='count')

    for ROI in ROI_list:
        df_ROI_scount['loca'] = df_ROI_scount['loca'].replace({ROI: label_name_nsujet[ROI]})

    localist_thresh = df_ROI_scount.query(f"count >= {sujet_thresh}")['loca'].values.tolist()

    ROI_plot_short_list = []

    for _ROI_thresh in localist_thresh:
        for _ROI_short in ROI_short_list:
            if _ROI_thresh.find(_ROI_short) != -1:
                ROI_plot_short_list.append(_ROI_thresh)

    #### plot
        #### allROI
    phase_cycle = 'whole'

    for band in freq_band_dict:
        
        _df_plot = df_R_Pxx.query(f"band == '{band}' and term != '(Intercept)' and phase_cycle == '{phase_cycle}'").copy()
        _df_plot["sig"] = _df_plot["pvalue"].apply(p_to_stars)

        fig = px.bar(_df_plot, x="ROI", y="estimate", color="term", barmode="group", text="sig", title=f"{band} {phase_cycle} allROI")

        fig.update_layout(xaxis_tickangle=-45, yaxis_title="estimate", xaxis_title="ROI", legend_title="term")
        fig.update_traces(textposition="outside", cliponaxis=False)

        # fig.show()

        outdir = os.path.join(path_results, "LMM", "fig", 'Pxx')
        fig.write_html(os.path.join(outdir, f"{phase_cycle}_{band}_allROI_barplot_LMM.html"),
                    include_plotlyjs="cdn")

        #### signiROI 
    phase_cycle = 'whole'
    
    for band in freq_band_dict:
        
        ROI_signi_sel = df_R_Pxx.query(f"band == '{band}' and term == 'respoc:statechl' and pvalue < 0.05 and phase_cycle == '{phase_cycle}'")["ROI"].values
        _df_plot = df_R_Pxx.query(f"band == '{band}' and term != '(Intercept)' and ROI in {ROI_signi_sel.tolist()}").copy()

        _df_plot["sig"] = _df_plot["pvalue"].apply(p_to_stars)

        fig = px.bar(_df_plot, x="ROI", y="estimate", color="term", barmode="group", text="sig", title=f"{band} {phase_cycle} signiROI")

        fig.update_layout(xaxis_tickangle=-45, yaxis_title="estimate", xaxis_title="ROI", legend_title="term")
        fig.update_traces(textposition="outside", cliponaxis=False)

        # fig.show()

        outdir = os.path.join(path_results, "LMM", "fig", 'Pxx')
        fig.write_html(os.path.join(outdir, f"{phase_cycle}_{band}_signiROI_barplot_LMM.html"),
                include_plotlyjs="cdn")

        #### threshROI
    phase_cycle = 'whole'

    for band in freq_band_dict:
        
        _df_plot = df_R_Pxx.query(f"ROI in {localist_thresh} and band == '{band}' and term != '(Intercept)' and phase_cycle == '{phase_cycle}'").copy()
        _df_plot["sig"] = _df_plot["pvalue"].apply(p_to_stars)

        fig = px.bar(_df_plot, x="ROI", y="estimate", color="term", barmode="group", text="sig", title=f"{band} {phase_cycle} threshROI")

        fig.update_layout(xaxis_tickangle=-45, yaxis_title="estimate", xaxis_title="ROI", legend_title="term")
        fig.update_traces(textposition="outside", cliponaxis=False)

        # fig.show()

        outdir = os.path.join(path_results, "LMM", "fig", 'Pxx')
        fig.write_html(os.path.join(outdir, f"{phase_cycle}_{band}_threshROI_barplot_LMM.html"),
                    include_plotlyjs="cdn")


        #### inspi/expi
    color_map = {"inspi": "steelblue",
                "expi": "darkorange",}

    LMM_param_list = ['respoc', 'statechl', 'respoc:statechl']

    for band in freq_band_dict:

        fig = make_subplots(
            rows=len(LMM_param_list),
            cols=1,
            shared_xaxes=True,
            subplot_titles=LMM_param_list
        )

        for r, LMM_param in enumerate(LMM_param_list, start=1):

            _df_plot = df_R_Pxx.query(f"band == '{band}' and term == '{LMM_param}' and phase_cycle != 'whole' and ROI in {ROI_plot_short_list}").copy()

            _df_plot["sig"] = _df_plot["pvalue"].apply(p_to_stars)

            for phase_cycle in _df_plot["phase_cycle"].unique():
                df_term = _df_plot[_df_plot["phase_cycle"] == phase_cycle]

                fig.add_bar(
                    x=df_term["ROI"],
                    y=df_term["estimate"],
                    name=phase_cycle,
                    text=df_term["sig"],
                    textposition="outside",
                    marker_color=color_map[phase_cycle],
                    row=r,
                    col=1
                )

        fig.update_layout(
            title=f"{band} all ROI",
            barmode="group",
            xaxis_tickangle=-45,
            yaxis_title="estimate",
            legend_title="phase_cycle",
            height=300 * len(phase_cycle_list)  # scale height
        )

        # fig.show()

        outdir = os.path.join(path_results, "LMM", "fig", "Pxx", "summary")
        fig.write_html(
            os.path.join(outdir, f"{band}_threshROI_LMM.html"),
            include_plotlyjs="cdn"
        )













########################
######## OLSa ########
########################





def export_res_OLSa():

    #### params
    sujet_thresh = 3
    pre_post_list = ['pre', 'post']
    ROI_short_list = ['Amygdala', 'insula-ant', 'lateralorbitofrontal', 'Hippocampus', 'insula-pos', 'medialorbitofrontal', 'insula', 'postcentral',
                      'precentral', 'Brain-Stem']
    rf_metric_short_list = ['cycle_freq', 'amplitude', 'oc_ratio', 'oc_val']

    #### load df
    df_R_OLSa = generate_all_df_OLSa()

    df_loca_allsujet = get_df_loca_allsujet()
    ROI_list = [_loca for _loca in df_loca_allsujet['loca'].unique() if _loca.find('UNSORTED') == -1]

    label_name_nsujet = {}
    
    for ROI in ROI_list:

        _nsujet = df_loca_allsujet.query(f"loca == '{ROI}'")['sujet'].unique().size
        label_name_nsujet[ROI] = f"{ROI} s({_nsujet})"

    for ROI in ROI_list:
        df_R_OLSa['ROI'] = df_R_OLSa['ROI'].replace({ROI: label_name_nsujet[ROI]})

    df_ROI_scount = df_loca_allsujet.groupby('loca').nunique('sujet')['sujet'].reset_index(name='count')

    for ROI in ROI_list:
        df_ROI_scount['loca'] = df_ROI_scount['loca'].replace({ROI: label_name_nsujet[ROI]})

    localist_thresh = df_ROI_scount.query(f"count >= {sujet_thresh}")['loca'].values.tolist()

    ROI_plot_short_list = []

    for _ROI_thresh in localist_thresh:
        for _ROI_short in ROI_short_list:
            if _ROI_thresh.find(_ROI_short) != -1:
                ROI_plot_short_list.append(_ROI_thresh)

    #### plot
        #### allROI
    phase_cycle = 'whole'

    for phase_cycle in phase_cycle_list:

        for pre_post in pre_post_list:

            for rf_metric in rf_metric_list:

                for band in freq_band_dict:
                    
                    _df_plot = df_R_OLSa.query(f"pre_post == '{pre_post}' and rf_metric == '{rf_metric}' and band == '{band}' and term != '(Intercept)' and phase_cycle == '{phase_cycle}'").copy()
                    _df_plot["sig"] = _df_plot["pvalue"].apply(p_to_stars)

                    fig = px.bar(_df_plot, x="ROI", y="estimate", color="term", barmode="group", text="sig", title=f"{rf_metric} {pre_post} {band} {phase_cycle} allROI")

                    fig.update_layout(xaxis_tickangle=-45, yaxis_title="estimate", xaxis_title="ROI", legend_title="term")
                    fig.update_traces(textposition="outside", cliponaxis=False)

                    # fig.show()

                    outdir = os.path.join(path_results, "LMM", "fig", 'OLSa', 'allfig')
                    fig.write_html(os.path.join(outdir, f"{rf_metric}_{pre_post}_{phase_cycle}_{band}_allROI_barplot_LMM.html"),
                                include_plotlyjs="cdn")

        #### signiROI 
    phase_cycle = 'whole'

    phase_cycle = 'whole'

    for phase_cycle in phase_cycle_list:

        for pre_post in pre_post_list:

            for rf_metric in rf_metric_list:

                for band in freq_band_dict:

                    ROI_signi_sel = df_R_OLSa.query(f"pre_post == '{pre_post}' and rf_metric == '{rf_metric}' and band == '{band}' and term == 'respoc:statechl' and pvalue < 0.05 and phase_cycle == '{phase_cycle}'")["ROI"].values
                    
                    _df_plot = df_R_OLSa.query(f"pre_post == '{pre_post}' and rf_metric == '{rf_metric}' and band == '{band}' and term != '(Intercept)' and phase_cycle == '{phase_cycle}' and ROI in {ROI_signi_sel.tolist()}").copy()
                    _df_plot["sig"] = _df_plot["pvalue"].apply(p_to_stars)

                    fig = px.bar(_df_plot, x="ROI", y="estimate", color="term", barmode="group", text="sig", title=f"{rf_metric} {pre_post} {band} {phase_cycle} allROI")

                    fig.update_layout(xaxis_tickangle=-45, yaxis_title="estimate", xaxis_title="ROI", legend_title="term")
                    fig.update_traces(textposition="outside", cliponaxis=False)

                    # fig.show()

                    outdir = os.path.join(path_results, "LMM", "fig", 'OLSa', 'allfig')
                    fig.write_html(os.path.join(outdir, f"{rf_metric}_{pre_post}_{phase_cycle}_{band}_signiROI_barplot_LMM.html"),
                include_plotlyjs="cdn")

        #### threshROI
    phase_cycle = 'whole'

    for phase_cycle in phase_cycle_list:

        for pre_post in pre_post_list:

            for rf_metric in rf_metric_list:

                for band in freq_band_dict:
                    
                    _df_plot = df_R_OLSa.query(f"ROI in {localist_thresh} and pre_post == '{pre_post}' and rf_metric == '{rf_metric}' and band == '{band}' and term != '(Intercept)' and phase_cycle == '{phase_cycle}'").copy()
                    _df_plot["sig"] = _df_plot["pvalue"].apply(p_to_stars)

                    fig = px.bar(_df_plot, x="ROI", y="estimate", color="term", barmode="group", text="sig", title=f"{rf_metric} {pre_post} {band} {phase_cycle} allROI")

                    fig.update_layout(xaxis_tickangle=-45, yaxis_title="estimate", xaxis_title="ROI", legend_title="term")
                    fig.update_traces(textposition="outside", cliponaxis=False)

                    # fig.show()

                    outdir = os.path.join(path_results, "LMM", "fig", 'OLSa', 'allfig')
                    fig.write_html(os.path.join(outdir, f"{rf_metric}_{pre_post}_{phase_cycle}_{band}_allROI_barplot_LMM.html"),
                                include_plotlyjs="cdn")

        #### inspi/expi
    color_map = {"inspi": "steelblue",
                "expi": "darkorange",}

    LMM_param_list = ['respoc', 'statechl', 'respoc:statechl']

    for rf_metric in rf_metric_list:

        for pre_post in pre_post_list:

            for band in freq_band_dict:

                fig = make_subplots(
                    rows=len(LMM_param_list),
                    cols=1,
                    shared_xaxes=True,
                    subplot_titles=LMM_param_list
                )

                for r, LMM_param in enumerate(LMM_param_list, start=1):

                    _df_plot = df_R_OLSa.query(f"pre_post == '{pre_post}' and rf_metric == '{rf_metric}' and band == '{band}' and term == '{LMM_param}' and phase_cycle != 'whole' and ROI in {localist_thresh}").copy()

                    _df_plot["sig"] = _df_plot["pvalue"].apply(p_to_stars)

                    for phase_cycle in _df_plot["phase_cycle"].unique():
                        df_term = _df_plot[_df_plot["phase_cycle"] == phase_cycle]

                        fig.add_bar(
                            x=df_term["ROI"],
                            y=df_term["estimate"],
                            name=phase_cycle,
                            text=df_term["sig"],
                            textposition="outside",
                            marker_color=color_map[phase_cycle],
                            row=r,
                            col=1
                        )

                fig.update_layout(
                    title=f"{pre_post} {rf_metric} {band}",
                    barmode="group",
                    xaxis_tickangle=-45,
                    yaxis_title="estimate",
                    legend_title="phase_cycle",
                    height=300 * len(phase_cycle_list)  # scale height
                )

                # fig.show()

                outdir = os.path.join(path_results, "LMM", "fig", "OLSa", "inspi_expi")
                fig.write_html(
                    os.path.join(outdir, f"{band}_threshROI_LMM.html"),
                    include_plotlyjs="cdn"
                )

    #### inspect

    if debug:


        df_R_OLSa.query(f"phase_cycle == 'whole' and pre_post == 'pre' and band == 'gamma' and term == 'respoc' and ROI in 'Amygdala s(8)'")
        df_R_OLSa.query(f"phase_cycle == 'whole' and pre_post == 'post' and band == 'gamma' and term == 'respoc' and ROI in 'Amygdala s(8)'")


    #### specific ROI
    color_map = {"inspi": "steelblue",
                "expi": "darkorange",}

    LMM_param_list = ['respoc', 'statechl', 'respoc:statechl']

    for pre_post in pre_post_list:

        for band in freq_band_dict:

            fig = make_subplots(
                rows=len(rf_metric_short_list),
                cols=len(LMM_param_list),
                shared_xaxes=True,
                subplot_titles=LMM_param_list
            )

            for r, rf_metric in enumerate(rf_metric_short_list, start=1):

                fig.add_annotation(
                                text=rf_metric,
                                xref="paper",
                                yref="y domain",
                                x=-1,              # move left of plots
                                y=0.5,
                                showarrow=False,
                                font=dict(size=14),
                                textangle=-90,
                                row=r,
                                col=1
                            )


                for c, LMM_param in enumerate(LMM_param_list, start=1):

                    _df_plot = df_R_OLSa.query(f"phase_cycle != 'whole' and pre_post == '{pre_post}' and rf_metric == '{rf_metric}' and band == '{band}' and term == '{LMM_param}' and ROI in {ROI_plot_short_list}").copy()

                    _df_plot["sig"] = _df_plot["pvalue"].apply(p_to_stars)

                    for phase_cycle in _df_plot["phase_cycle"].unique():

                        df_term = _df_plot.query(f"phase_cycle == '{phase_cycle}'")

                        fig.add_bar(
                            x=df_term["ROI"],
                            y=df_term["estimate"],
                            name=phase_cycle,
                            text=df_term["sig"],
                            textposition="outside",
                            marker_color=color_map[phase_cycle],
                            row=r,
                            col=c
                        )

            fig.update_layout(
                title=f"{pre_post} {band}",
                barmode="group",
                xaxis_tickangle=-45,
                height=300 * len(phase_cycle_list)  # scale height
            )

            # fig.show()

            outdir = os.path.join(path_results, "LMM", "fig", "OLSa", "summary")
            fig.write_html(
                os.path.join(outdir, f"{pre_post}_{band}_threshROI_LMM.html"),
                include_plotlyjs="cdn"
            )













################################
######## EXECUTE ########
################################


if __name__ == '__main__':



    export_res_Pxx()

    export_res_OLSa()
        


    























