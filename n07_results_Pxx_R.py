



from n00_config_params import *
from n00bis_config_analysis_functions import *
from n01_manip_data import *


import plotly.express as px










################################
######## GENERATE ALLDF ########
################################

def generate_all_df_Pxx():

    os.chdir(os.path.join(path_results, 'LMM', 'df'))

    df_loca_allsujet = get_df_loca_allsujet()
    ROI_list = [_loca for _loca in df_loca_allsujet['loca'].unique() if _loca.find('UNSORTED') == -1]

    df_R_Pxx = []

    for band in freq_band_dict:

        for ROI in ROI_list:

            _df = pd.read_excel(f"{band}_Pxx_lmm_{ROI}_res.xlsx")
            _df_add = pd.concat([pd.DataFrame({'band' : [band]*_df.shape[0], 'ROI' : [ROI]*_df.shape[0]}), _df[['term', 'estimate', 'p.value']]], axis=1)
            
            df_R_Pxx.append(_df_add)

    df_R_Pxx = pd.concat(df_R_Pxx)
    df_R_Pxx = df_R_Pxx.rename(columns={'p.value' : 'pvalue'})

    return df_R_Pxx




def generate_all_df_OLSa():

    os.chdir(os.path.join(path_results, 'LMM', 'df'))

    df_loca_allsujet = get_df_loca_allsujet()
    ROI_list = [_loca for _loca in df_loca_allsujet['loca'].unique() if _loca.find('UNSORTED') == -1]

    df_R_OLSa = []

    rf_metrics_tot_sel = ['cycle_duration', 'inspi_duration', 'expi_duration', 'cycle_freq', 'inspi_volume',
       'expi_volume', 'total_amplitude', 'inspi_amplitude', 'expi_amplitude',
       'total_volume']
    rf_metrics_whole_cycle = ['cycle_duration', 'cycle_freq', 'total_amplitude', 'total_volume']
    rf_metrics_phase_cycle = ['inspi_duration', 'expi_duration', 'inspi_volume', 'expi_volume', 'inspi_amplitude', 'expi_amplitude']

    for rf_metric in rf_metrics_whole_cycle:

        for band in freq_band_dict:

            for ROI in ROI_list:

                _df = pd.read_excel(f"{band}_{rf_metric}_OLS_a_lmm_{ROI}_res.xlsx")
                _df_add = pd.concat([pd.DataFrame({'band' : [band]*_df.shape[0], 'rf_metric' : [rf_metric]*_df.shape[0], 'ROI' : [ROI]*_df.shape[0]}), _df[['term', 'estimate', 'p.value']]], axis=1)
                
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


    #### load df
    df_R_Pxx = generate_all_df_Pxx()

    df_loca_allsujet = get_df_loca_allsujet()
    ROI_list = [_loca for _loca in df_loca_allsujet['loca'].unique() if _loca.find('UNSORTED') == -1]

    label_name_nsujet = {}
    
    for ROI in ROI_list:

        _nsujet = df_loca_allsujet.query(f"loca == '{ROI}'")['sujet'].unique().size
        label_name_nsujet[ROI] = f"{ROI} s({_nsujet})"

    for ROI in ROI_list:
        df_R_Pxx['ROI'] = df_R_Pxx['ROI'].replace({ROI: label_name_nsujet[ROI]})



    #### plot
        #### all ROI
    band = 'theta'

    for band in freq_band_dict:

        _df_plot = df_R_Pxx.query(f"band == '{band}' and term != '(Intercept)'")
        _df_plot['sig'] = _df_plot['pvalue'].apply(p_to_stars).copy()

        g = sns.catplot(kind='bar', data=_df_plot, x='ROI', y='estimate', hue='term')
        g.set_xticklabels(rotation=45, ha='right')

        for ax in g.axes.flat:
            for patch, (_, row) in zip(ax.patches, _df_plot.iterrows()):
                if row['sig'] != '':
                    height = patch.get_height()
                    ax.text(
                        patch.get_x() + patch.get_width() / 2,
                        height,
                        row['sig'],
                        ha='center',
                        va='bottom',
                        fontsize=12,
                        color='black'
                    )
        plt.title(band)

        # plt.show()

        os.chdir(os.path.join(path_results, 'LMM', 'fig', 'Pxx'))
        g.savefig(f"{band}_allROI_barplot_LMM.jpeg")
    

        #### signi ROI
    band = 'theta'

    for band in freq_band_dict:

        ROI_signi_sel = df_R_Pxx.query(f"band == '{band}' and term == 'respoc:statechl' and pvalue < 0.05")['ROI'].values
        _df_plot = df_R_Pxx.query(f"band == '{band}' and term != '(Intercept)' and ROI in {ROI_signi_sel.tolist()}")
        _df_plot['sig'] = _df_plot['pvalue'].apply(p_to_stars).copy()

        g = sns.catplot(kind='bar', data=_df_plot, x='ROI', y='estimate', hue='term')
        g.set_xticklabels(rotation=45, ha='right')

        for ax in g.axes.flat:
            for patch, (_, row) in zip(ax.patches, _df_plot.iterrows()):
                if row['sig'] != '':
                    height = patch.get_height()
                    ax.text(
                        patch.get_x() + patch.get_width() / 2,
                        height,
                        row['sig'],
                        ha='center',
                        va='bottom',
                        fontsize=12,
                        color='black'
                    )
        plt.title(band)

        # plt.show()

        os.chdir(os.path.join(path_results, 'LMM', 'fig', 'Pxx'))
        g.savefig(f"{band}_signiROI_barplot_LMM.jpeg")





    #### plotly

    for band in freq_band_dict:
        
        _df_plot = df_R_Pxx.query(f"band == '{band}' and term != '(Intercept)'").copy()
        _df_plot["sig"] = _df_plot["pvalue"].apply(p_to_stars)

        fig = px.bar(_df_plot, x="ROI", y="estimate", color="term", barmode="group", text="sig", title=band)

        fig.update_layout(xaxis_tickangle=-45, yaxis_title="estimate", xaxis_title="ROI", legend_title="term")
        fig.update_traces(textposition="outside", cliponaxis=False)

        # fig.show()

        outdir = os.path.join(path_results, "LMM", "fig", 'Pxx')
        fig.write_html(os.path.join(outdir, f"{band}_allROI_barplot_LMM.html"),
                    include_plotlyjs="cdn")
        
        
    for band in freq_band_dict:
        
        ROI_signi_sel = df_R_Pxx.query(f"band == '{band}' and term == 'respoc:statechl' and pvalue < 0.05")["ROI"].values
        _df_plot = df_R_Pxx.query(f"band == '{band}' and term != '(Intercept)' and ROI in {ROI_signi_sel.tolist()}").copy()

        _df_plot["sig"] = _df_plot["pvalue"].apply(p_to_stars)

        fig = px.bar(_df_plot, x="ROI", y="estimate", color="term", barmode="group", text="sig", title=band)

        fig.update_layout(xaxis_tickangle=-45, yaxis_title="estimate", xaxis_title="ROI", legend_title="term")
        fig.update_traces(textposition="outside", cliponaxis=False)

        # fig.show()

        outdir = os.path.join(path_results, "LMM", "fig", 'Pxx')
        fig.write_html(os.path.join(outdir, f"{band}_signiROI_barplot_LMM.html"),
                   include_plotlyjs="cdn")






########################
######## OLSa ########
########################




def export_res_OLSa():


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

    rf_metrics = df_R_OLSa['rf_metric'].unique()


    #### plot
        #### all ROI
    band = 'theta'

    for band in freq_band_dict:

        for rf_metric in rf_metrics:

            _df_plot = df_R_OLSa.query(f"band == '{band}' and term != '(Intercept)' and rf_metric == '{rf_metric}'")
            _df_plot['sig'] = _df_plot['pvalue'].apply(p_to_stars).copy()

            g = sns.catplot(kind='bar', data=_df_plot, x='ROI', y='estimate', hue='term')
            g.set_xticklabels(rotation=45, ha='right')

            for ax in g.axes.flat:
                for patch, (_, row) in zip(ax.patches, _df_plot.iterrows()):
                    if row['sig'] != '':
                        height = patch.get_height()
                        ax.text(
                            patch.get_x() + patch.get_width() / 2,
                            height,
                            row['sig'],
                            ha='center',
                            va='bottom',
                            fontsize=12,
                            color='black'
                        )
            plt.title(f"{band}, {rf_metric}")

            # plt.show()

            os.chdir(os.path.join(path_results, 'LMM', 'fig'))
            g.savefig(f"OLSa_{band}_allROI_barplot_LMM.jpeg")
    

        #### signi ROI
    band = 'theta'

    for band in freq_band_dict:

        for rf_metric in rf_metrics:

            ROI_signi_sel = df_R_OLSa.query(f"band == '{band}' and term == 'respoc:statechl' and pvalue < 0.05 and rf_metric == '{rf_metric}'")['ROI'].values
            _df_plot = df_R_OLSa.query(f"band == '{band}' and term != '(Intercept)' and ROI in {ROI_signi_sel.tolist()} and rf_metric == '{rf_metric}'")
            _df_plot['sig'] = _df_plot['pvalue'].apply(p_to_stars).copy()

            g = sns.catplot(kind='bar', data=_df_plot, x='ROI', y='estimate', hue='term')
            g.set_xticklabels(rotation=45, ha='right')

            for ax in g.axes.flat:
                for patch, (_, row) in zip(ax.patches, _df_plot.iterrows()):
                    if row['sig'] != '':
                        height = patch.get_height()
                        ax.text(
                            patch.get_x() + patch.get_width() / 2,
                            height,
                            row['sig'],
                            ha='center',
                            va='bottom',
                            fontsize=12,
                            color='black'
                        )
            plt.title(f"{band}, {rf_metric}")

            # plt.show()

            os.chdir(os.path.join(path_results, 'LMM', 'fig'))
            g.savefig(f"OLSa_{band}_signiROI_barplot_LMM.jpeg")





    #### plotly

    for band in freq_band_dict:

        for rf_metric in rf_metrics:
        
            _df_plot = df_R_OLSa.query(f"band == '{band}' and term != '(Intercept)' and rf_metric == '{rf_metric}'").copy()
            _df_plot["sig"] = _df_plot["pvalue"].apply(p_to_stars)

            fig = px.bar(_df_plot, x="ROI", y="estimate", color="term", barmode="group", text="sig", title=band)

            fig.update_layout(xaxis_tickangle=-45, yaxis_title="estimate", xaxis_title="ROI", legend_title="term")
            fig.update_traces(textposition="outside", cliponaxis=False)

            # fig.show()

            outdir = os.path.join(path_results, "LMM", "fig")
            fig.write_html(os.path.join(outdir, 'OLSa', f"OLSa_{band}_{rf_metric}_allROI_barplot_LMM.html"),
                        include_plotlyjs="cdn")
        
        
    for band in freq_band_dict:

        for rf_metric in rf_metrics:
        
            ROI_signi_sel = df_R_OLSa.query(f"band == '{band}' and term == 'respoc:statechl' and pvalue < 0.05 and rf_metric == '{rf_metric}'")["ROI"].values
            _df_plot = df_R_OLSa.query(f"band == '{band}' and term != '(Intercept)' and ROI in {ROI_signi_sel.tolist()}").copy()

            _df_plot["sig"] = _df_plot["pvalue"].apply(p_to_stars)

            fig = px.bar(_df_plot, x="ROI", y="estimate", color="term", barmode="group", text="sig", title=band)

            fig.update_layout(xaxis_tickangle=-45, yaxis_title="estimate", xaxis_title="ROI", legend_title="term")
            fig.update_traces(textposition="outside", cliponaxis=False)

            # fig.show()

            outdir = os.path.join(path_results, "LMM", "fig")
            fig.write_html(os.path.join(outdir, 'OLSa', f"OLSa_{band}_{rf_metric}_signiROI_barplot_LMM.html"),
                    include_plotlyjs="cdn")









################################
######## EXECUTE ########
################################


if __name__ == '__main__':



    export_res_Pxx()

    export_res_OLSa()
        


    























