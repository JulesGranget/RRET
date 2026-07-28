



from A_config.n01_O_params import *
from A_config.n02_O_analysis_functions import *
from A_config.n03_X_manip_data import *


import plotly.express as px
from plotly.subplots import make_subplots









################################
######## GENERATE ALLDF ########
################################

def generate_all_df_Pxx_wholecycle():

    phase_cycle_list = ['inspi', 'expi']

    os.chdir(os.path.join(path_results, 'LMM', 'df'))

    df_R_Pxx = []

    for phase_cycle in phase_cycle_list:

        for band in freq_band_dict:

            for ROI in ROI_short_list:

                _df = pd.read_excel(f"RES_{band}_{phase_cycle}_{ROI}_Pxx.xlsx")
                _df_add = pd.concat([pd.DataFrame({'phase_cycle' : [phase_cycle]*_df.shape[0], 'band' : [band]*_df.shape[0], 'ROI' : [ROI]*_df.shape[0]}), _df[['term', 'estimate', 'p.value']]], axis=1)
                
                df_R_Pxx.append(_df_add)

    df_R_Pxx = pd.concat(df_R_Pxx)
    df_R_Pxx = df_R_Pxx.rename(columns={'p.value' : 'pvalue'})

    return df_R_Pxx




def generate_all_df_Pxx_MECACO2():

    phase_cycle_list = ['inspi', 'expi']

    os.chdir(os.path.join(path_results, 'LMM', 'df'))

    df_R_Pxx = []

    for phase_cycle in phase_cycle_list:

        for band in freq_band_dict:

            for ROI in ROI_short_list_MECACO2:

                _df = pd.read_excel(f"MECACO2_{band}_{phase_cycle}_{ROI}_LMM.xlsx")
                _df_add = pd.concat([pd.DataFrame({'phase_cycle' : [phase_cycle]*_df.shape[0], 'band' : [band]*_df.shape[0], 'ROI' : [ROI]*_df.shape[0]}), _df[['term', 'estimate', 'p.value']]], axis=1)
                
                df_R_Pxx.append(_df_add)

    df_R_Pxx = pd.concat(df_R_Pxx)
    df_R_Pxx = df_R_Pxx.rename(columns={'p.value' : 'pvalue'})

    return df_R_Pxx


    


def generate_all_df_Pxx_MECACO2_REG():

    phase_cycle_list = ['inspi', 'expi']

    os.chdir(os.path.join(path_results, 'LMM', 'MECACO2', 'df'))

    df_R_Pxx = []

    for phase_cycle in phase_cycle_list:

        for band in freq_band_dict:

            for ROI in ROI_short_list_MECACO2:

                _df = pd.read_excel(f"RES_{phase_cycle}_{band}_{ROI}_Pxx_OC_LMM.xlsx")
                _df_add = pd.concat([pd.DataFrame({'phase_cycle' : [phase_cycle]*_df.shape[0], 'band' : [band]*_df.shape[0], 'ROI' : [ROI]*_df.shape[0]}), _df[['term', 'estimate', 'p.value']]], axis=1)
                
                df_R_Pxx.append(_df_add)

    df_R_Pxx = pd.concat(df_R_Pxx)
    df_R_Pxx = df_R_Pxx.rename(columns={'p.value' : 'pvalue'})

    return df_R_Pxx


    





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




def export_res_Pxx_general():

    phase_cycle_list = ['inspi', 'expi']

    #### load df
    df_R_Pxx = generate_all_df_Pxx_wholecycle()

    df_loca_allsujet = get_df_loca_allsujet()

    label_name_nsujet = {}
    
    for ROI in ROI_short_list:

        _nsujet = df_loca_allsujet.query(f"loca == '{ROI}'")['sujet'].unique().size
        label_name_nsujet[ROI] = f"{ROI} s({_nsujet})"

    for ROI in ROI_short_list:
        df_R_Pxx['ROI'] = df_R_Pxx['ROI'].replace({ROI: label_name_nsujet[ROI]})

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
        #### allROI
    # phase_cycle = 'whole'

    # for band in freq_band_dict:
        
    #     _df_plot = df_R_Pxx.query(f"band == '{band}' and term != '(Intercept)' and phase_cycle == '{phase_cycle}'").copy()
    #     _df_plot["sig"] = _df_plot["pvalue"].apply(p_to_stars)

    #     fig = px.bar(_df_plot, x="ROI", y="estimate", color="term", barmode="group", text="sig", title=f"{band} {phase_cycle} allROI")

    #     fig.update_layout(xaxis_tickangle=-45, yaxis_title="estimate", xaxis_title="ROI", legend_title="term")
    #     fig.update_traces(textposition="outside", cliponaxis=False)

    #     # fig.show()

    #     outdir = os.path.join(path_results, "LMM", "fig", 'Pxx')
    #     fig.write_html(os.path.join(outdir, f"{phase_cycle}_{band}_allROI_barplot_LMM.html"),
    #                 include_plotlyjs="cdn")

    #     #### signiROI 
    # phase_cycle = 'whole'
    
    # for band in freq_band_dict:
        
    #     ROI_signi_sel = df_R_Pxx.query(f"band == '{band}' and term == 'respoc:statechl' and pvalue < 0.05 and phase_cycle == '{phase_cycle}'")["ROI"].values
    #     _df_plot = df_R_Pxx.query(f"band == '{band}' and term != '(Intercept)' and ROI in {ROI_signi_sel.tolist()}").copy()

    #     _df_plot["sig"] = _df_plot["pvalue"].apply(p_to_stars)

    #     fig = px.bar(_df_plot, x="ROI", y="estimate", color="term", barmode="group", text="sig", title=f"{band} {phase_cycle} signiROI")

    #     fig.update_layout(xaxis_tickangle=-45, yaxis_title="estimate", xaxis_title="ROI", legend_title="term")
    #     fig.update_traces(textposition="outside", cliponaxis=False)

    #     # fig.show()

    #     outdir = os.path.join(path_results, "LMM", "fig", 'Pxx')
    #     fig.write_html(os.path.join(outdir, f"{phase_cycle}_{band}_signiROI_barplot_LMM.html"),
    #             include_plotlyjs="cdn")

    #     #### threshROI
    # phase_cycle = 'whole'

    # for band in freq_band_dict:
        
    #     _df_plot = df_R_Pxx.query(f"ROI in {localist_thresh} and band == '{band}' and term != '(Intercept)' and phase_cycle == '{phase_cycle}'").copy()
    #     _df_plot["sig"] = _df_plot["pvalue"].apply(p_to_stars)

    #     fig = px.bar(_df_plot, x="ROI", y="estimate", color="term", barmode="group", text="sig", title=f"{band} {phase_cycle} threshROI")

    #     fig.update_layout(xaxis_tickangle=-45, yaxis_title="estimate", xaxis_title="ROI", legend_title="term")
    #     fig.update_traces(textposition="outside", cliponaxis=False)

    #     # fig.show()

    #     outdir = os.path.join(path_results, "LMM", "fig", 'Pxx')
    #     fig.write_html(os.path.join(outdir, f"{phase_cycle}_{band}_threshROI_barplot_LMM.html"),
    #                 include_plotlyjs="cdn")


        #### inspi/expi
    color_map = {
        "inspi": "#1f77b4",
        "expi":  "#d62728" 
    }

    LMM_param_list = ['respoc', 'statechl', 'respoc:statechl']
    
    ROI_order = ['Amygdala s(8)', 'Hippocampus s(7)', 'insula-ant s(5)', 'insula-pos s(4)', 'lateralorbitofrontal s(5)', 'medialorbitofrontal s(5)', 'postcentral s(3)', 'precentral s(4)']

    #band = 'theta'
    for band in freq_band_dict:


        fig = make_subplots(
            rows=len(LMM_param_list),
            cols=1,
            shared_xaxes=True,
            subplot_titles=LMM_param_list
        )

        y_max = df_R_Pxx.query(f"band == '{band}' and term != '(Intercept)' and phase_cycle != 'whole' and ROI in {ROI_plot_short_list}").copy()["estimate"].abs().max()

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

                margin = y_max * 0.50   # 15% extra space

                fig.update_yaxes(
                    range=[-y_max - margin, y_max + margin],
                    row=r,
                    col=1
                )

                fig.update_xaxes(
                    categoryorder="array",
                    categoryarray=ROI_order,
                    row=r,
                    col=1
                )

        fig.update_layout(
            title=f"{band} all ROI",
            template="simple_white",
            barmode="group",
            xaxis_tickangle=-45,
            yaxis_title="estimate",
            legend_title="phase_cycle",
            height=320 * len(phase_cycle_list),  # scale height
            width=500,
        )

        # fig.show()

        outdir = os.path.join(path_results, "LMM", "fig", "Pxx", "summary")
        fig.write_html(
            os.path.join(outdir, f"{band}_threshROI_LMM.html"),
            include_plotlyjs="cdn"
        )

        fig.write_image(
            os.path.join(path_paper_figure_export, f"{band}_threshROI_LMM.svg")
        )

    #### export df data
    df_export = df_R_Pxx.query(f"phase_cycle != 'whole' and ROI in {ROI_plot_short_list} and term != '(Intercept)'")
    
    output_df = os.path.join(path_results, "LMM", "fig", "Pxx", "summary")
    filename = os.path.join(output_df, "df_Pxx_LMM.xlsx")
    df_export.to_excel(filename)









def export_res_Pxx_MECACO2():

    #### config
    sujet_thresh_MECACO2 = 1

    #### load df
    df_R_Pxx = generate_all_df_Pxx_MECACO2()

    df_loca_allsujet = get_df_loca_allsujet()
    df_loca_allsujet = df_loca_allsujet.query(f"sujet in {sujet_list_interaction_analysis}")

    label_name_nsujet = {}
    
    for ROI in ROI_short_list_MECACO2:

        _nsujet = df_loca_allsujet.query(f"loca == '{ROI}'")['sujet'].unique().size
        label_name_nsujet[ROI] = f"{ROI} s({_nsujet})"

    for ROI in ROI_short_list_MECACO2:
        df_R_Pxx['ROI'] = df_R_Pxx['ROI'].replace({ROI: label_name_nsujet[ROI]})

    df_ROI_scount = df_loca_allsujet.query(f"loca in {ROI_short_list_MECACO2}").groupby('loca').nunique('sujet')['sujet'].reset_index(name='count')

    for ROI in ROI_short_list_MECACO2:
        df_ROI_scount['loca'] = df_ROI_scount['loca'].replace({ROI: label_name_nsujet[ROI]})

    localist_thresh = df_ROI_scount.query(f"count >= {sujet_thresh_MECACO2}")['loca'].values.tolist()

    ROI_plot_short_list = []

    for _ROI_thresh in localist_thresh:
        for _ROI_short in ROI_short_list_MECACO2:
            if _ROI_thresh.find(_ROI_short) != -1:
                ROI_plot_short_list.append(_ROI_thresh)

        #### inspi/expi
    color_map = {
        "inspi": "#1f77b4",
        "expi":  "#d62728" 
    }

    LMM_param_list = ['oc_condoc', 'condBOTH', 'condCO2', 'condMECA',
       'oc_condoc:condBOTH', 'oc_condoc:condCO2', 'oc_condoc:condMECA']
    
    ROI_order = ['Amygdala s(3)', 'Hippocampus s(2)', 'insula-ant s(1)', 'lateralorbitofrontal s(2)', 'medialorbitofrontal s(1)', 'postcentral s(1)']

    #band = 'theta'
    for band in freq_band_dict:


        fig = make_subplots(
            rows=len(LMM_param_list),
            cols=1,
            shared_xaxes=True,
            subplot_titles=LMM_param_list
        )

        y_max = df_R_Pxx.query(f"band == '{band}' and term != '(Intercept)' and ROI in {ROI_plot_short_list}").copy()["estimate"].abs().max()

        for r, LMM_param in enumerate(LMM_param_list, start=1):

            _df_plot = df_R_Pxx.query(f"band == '{band}' and term == '{LMM_param}' and ROI in {ROI_plot_short_list}").copy()

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

                margin = y_max * 0.50   # 15% extra space

                fig.update_yaxes(
                    range=[-y_max - margin, y_max + margin],
                    row=r,
                    col=1
                )

                fig.update_xaxes(
                    categoryorder="array",
                    categoryarray=ROI_order,
                    row=r,
                    col=1
                )

        fig.update_layout(
            title=f"{band} all ROI",
            template="simple_white",
            barmode="group",
            xaxis_tickangle=-45,
            yaxis_title="estimate",
            legend_title="phase_cycle",
            height=320 * len(LMM_param_list),  # scale height
            width=500,
        )

        # fig.show()

        outdir = os.path.join(path_results, "LMM", "fig", "Pxx", "summary")
        fig.write_html(
            os.path.join(outdir, f"{band}_MECACO2_threshROI_LMM.html"),
            include_plotlyjs="cdn"
        )

        fig.write_image(
            os.path.join(path_paper_figure_export, f"{band}_MECACO2_threshROI_LMM.svg")
        )

    #### export df data
    df_export = df_R_Pxx.query(f"ROI in {ROI_plot_short_list} and term != '(Intercept)'")
    
    output_df = os.path.join(path_results, "LMM", "fig", "Pxx", "summary")
    filename = os.path.join(output_df, "df_Pxx_MECACO2_LMM.xlsx")
    df_export.to_excel(filename)






def export_res_Pxx_MECACO2_REG():

    #### config
    sujet_thresh_MECACO2 = 1

    #### load df
    df_R_Pxx = generate_all_df_Pxx_MECACO2_REG()

    df_loca_allsujet = get_df_loca_allsujet()
    df_loca_allsujet = df_loca_allsujet.query(f"sujet in {sujet_list_interaction_analysis}")

    label_name_nsujet = {}
    
    for ROI in ROI_short_list_MECACO2:

        _nsujet = df_loca_allsujet.query(f"loca == '{ROI}'")['sujet'].unique().size
        label_name_nsujet[ROI] = f"{ROI} s({_nsujet})"

    for ROI in ROI_short_list_MECACO2:
        df_R_Pxx['ROI'] = df_R_Pxx['ROI'].replace({ROI: label_name_nsujet[ROI]})

    df_ROI_scount = df_loca_allsujet.query(f"loca in {ROI_short_list_MECACO2}").groupby('loca').nunique('sujet')['sujet'].reset_index(name='count')

    for ROI in ROI_short_list_MECACO2:
        df_ROI_scount['loca'] = df_ROI_scount['loca'].replace({ROI: label_name_nsujet[ROI]})

    localist_thresh = df_ROI_scount.query(f"count >= {sujet_thresh_MECACO2}")['loca'].values.tolist()

    ROI_plot_short_list = []

    for _ROI_thresh in localist_thresh:
        for _ROI_short in ROI_short_list_MECACO2:
            if _ROI_thresh.find(_ROI_short) != -1:
                ROI_plot_short_list.append(_ROI_thresh)

        #### inspi/expi
    color_map = {
        "inspi": "#1f77b4",
        "expi":  "#d62728" 
    }

    LMM_param_list = ['condBOTH', 'condCO2', 'condMECA', 'oc_ratio',
       'condBOTH:oc_ratio', 'condCO2:oc_ratio', 'condMECA:oc_ratio']

    #### all ROI
    ROI_order = ['Amygdala s(3)', 'Hippocampus s(2)', 'insula-ant s(1)', 'lateralorbitofrontal s(2)', 'medialorbitofrontal s(1)', 'postcentral s(1)']

    #band = 'theta'
    for band in freq_band_dict:


        fig = make_subplots(
            rows=len(LMM_param_list),
            cols=1,
            shared_xaxes=True,
            subplot_titles=LMM_param_list
        )

        y_max = df_R_Pxx.query(f"band == '{band}' and term != '(Intercept)' and ROI in {ROI_plot_short_list}").copy()["estimate"].abs().max()

        for r, LMM_param in enumerate(LMM_param_list, start=1):

            _df_plot = df_R_Pxx.query(f"band == '{band}' and term == '{LMM_param}' and ROI in {ROI_plot_short_list}").copy()

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

                margin = y_max * 0.50   # 15% extra space

                fig.update_yaxes(
                    range=[-y_max - margin, y_max + margin],
                    row=r,
                    col=1
                )

                fig.update_xaxes(
                    categoryorder="array",
                    categoryarray=ROI_order,
                    row=r,
                    col=1
                )

        fig.update_layout(
            title=f"{band} all ROI",
            template="simple_white",
            barmode="group",
            xaxis_tickangle=-45,
            yaxis_title="estimate",
            legend_title="phase_cycle",
            height=320 * len(LMM_param_list),  # scale height
            width=500,
        )

        # fig.show()

        outdir = os.path.join(path_results, "LMM", "MECACO2", "fig")
        fig.write_html(
            os.path.join(outdir, f"{band}_REG_MECACO2_threshROI_LMM.html"),
            include_plotlyjs="cdn"
        )

        fig.write_image(
            os.path.join(path_paper_figure_export, f"{band}_REG_MECACO2_threshROI_LMM.svg")
        )

    #### only one ROI example reg
    outdir = os.path.join(path_results, "LMM", "MECACO2", "fig")

    filename = os.path.join(path_precompute, 'TF', 'session', 'df_R', f"df_reg_MECACO2_R.xlsx")
    df_reg_MECACO2 = pd.read_excel(filename)

    ROI_sel = 'Amygdala'

    df_plot_reg_one_ROI = df_reg_MECACO2.query(f"ROI == '{ROI_sel}'")

    band_sel = 'theta'
    sns.lmplot(df_plot_reg_one_ROI.query(f"band == '{band_sel}'"), x='oc_ratio', y='Pxx', hue='sujet', col='cond', col_order=['ctrl', 'MECA', 'CO2', 'BOTH'])
    plt.suptitle(band_sel)
    # plt.show()
    plt.savefig(os.path.join(outdir, f"{ROI_sel}_{band_sel}_reg_example.png"))

    df_median = df_plot_reg_one_ROI.drop(columns=['Unnamed: 0', 'chan_label', 'ROI']).groupby(['cond', 'band', 'phase_cycle', 'chan', 'sujet']).median().reset_index()
    sns.lmplot(df_median, x='oc_ratio', y='Pxx', hue='cond', col='band', hue_order=['ctrl', 'MECA', 'CO2', 'BOTH'])
    # plt.show()
    plt.savefig(os.path.join(outdir, f"{ROI_sel}_chanwise_reg_example.png"))

    plt.close('all')

    #### only one ROI LMM
    ROI_sel = 'Amygdala s(3)'
    ROI_title = 'AMYGDALA'
    band_list = list(freq_band_dict.keys())

    LMM_param_dict = {  'full' : ['oc_ratio', 'condMECA', 'condCO2', 'condBOTH', 
                                'condMECA:oc_ratio', 'condCO2:oc_ratio', 'condBOTH:oc_ratio'],
                        'COND' : ['condMECA', 'condCO2', 'condBOTH'],
                        'OC_RATIO' : ['oc_ratio', 'condMECA:oc_ratio', 'condCO2:oc_ratio', 'condBOTH:oc_ratio'],}

    for effect_type in LMM_param_dict:

        fig = make_subplots(
            rows=len(freq_band_dict),
            cols=1,
            shared_xaxes=False,
            subplot_titles=band_list
        )

        for r, band in enumerate(band_list, start=1):

            _df_plot = df_R_Pxx.query(f"band == '{band}' and ROI == '{ROI_sel}' and term in {LMM_param_dict[effect_type]}").copy()

            _df_plot["sig"] = _df_plot["pvalue"].apply(p_to_stars)

            for phase_cycle in _df_plot["phase_cycle"].unique():

                df_term = _df_plot[_df_plot["phase_cycle"] == phase_cycle]

                fig.add_bar(
                    x=df_term["term"],
                    y=df_term["estimate"],
                    name=phase_cycle,
                    text=df_term["sig"],
                    textposition="outside",
                    marker_color=color_map[phase_cycle],
                    row=r,
                    col=1
                )

                y_max = _df_plot["estimate"].abs().max()

                margin = y_max * 0.50   # 15% extra space

                fig.update_yaxes(
                    range=[-y_max - margin, y_max + margin],
                    row=r,
                    col=1
                )

                fig.update_xaxes(
                    categoryorder="array",
                    categoryarray=LMM_param_dict[effect_type],
                    row=r,
                    col=1
                )

            fig.update_layout(
            title=f"{ROI_title} all band",
            template="simple_white",
            barmode="group",
            yaxis_title="estimate",
            legend_title="phase_cycle",
            height=500 * len(LMM_param_dict[effect_type]),
            width=400,
            )

        fig.update_xaxes(tickangle=-45)

        # fig.show()

        outdir = os.path.join(path_results, "LMM", "MECACO2", "fig")
        fig.write_html(
            os.path.join(outdir, f"{ROI_title}_{effect_type}_REG_MECACO2_threshROI_LMM.html"),
            include_plotlyjs="cdn"
        )

        fig.write_image(
            os.path.join(path_paper_figure_export, f"{ROI_title}_{effect_type}_REG_MECACO2_threshROI_LMM.svg")
        )

    #### export df data
    df_export = df_R_Pxx.query(f"ROI == '{ROI_sel}' and term != '(Intercept)'")
    
    output_df = os.path.join(path_results, "LMM", "MECACO2", "fig")
    filename = os.path.join(output_df, f"df_Pxx_{ROI_title}_REG_MECACO2_LMM.xlsx")
    df_export.to_excel(filename)









################################
######## EXECUTE ########
################################


if __name__ == '__main__':

    export_res_Pxx_general()
    export_res_Pxx_MECACO2()
    export_res_Pxx_MECACO2_REG()
            
        


    























