



from A_config.n01_O_params import *
from A_config.n02_O_analysis_functions import *
from A_config.n03_X_manip_data import *







########################################
######## PLOT RESPI MEAN ########
########################################


def export_allcycles_fig_for_allcond():

    #### load
    df_resp_allsujet_list = []

    for sujet in sujet_list:

        if sujet == 'NS215':

            path_respfeatures_sujet = os.path.join(path_precompute, 'RESP', 'respfeatures', f"{sujet}_cycles_info_cleaned.xlsx")
            df_resp_cycle_cleaned_sujet = pd.read_excel(path_respfeatures_sujet)
            df_resp_cycle_cleaned_sujet['sujet'] = df_resp_cycle_cleaned_sujet.shape[0] * [sujet]
            max_val_load = df_resp_cycle_cleaned_sujet[['challengeLoadMagInh', 'challengeLoadMagExh']].values.max(axis=1)
            df_resp_cycle_cleaned_sujet['challengeLoadMag'] = max_val_load
            df_resp_cycle_cleaned_sujet = df_resp_cycle_cleaned_sujet.drop(columns=['challengeLoadMagInh', 'challengeLoadMagExh'])
            df_resp_allsujet_list.append(df_resp_cycle_cleaned_sujet)

        else:

            path_respfeatures_sujet = os.path.join(path_precompute, 'RESP', 'respfeatures', f"{sujet}_cycles_info_cleaned.xlsx")
            df_resp_cycle_cleaned_sujet = pd.read_excel(path_respfeatures_sujet)
            df_resp_cycle_cleaned_sujet['sujet'] = df_resp_cycle_cleaned_sujet.shape[0] * [sujet]
            df_resp_allsujet_list.append(df_resp_cycle_cleaned_sujet)

    df_resp_allsujet = pd.concat(df_resp_allsujet_list)
    df_resp_allsujet = df_resp_allsujet.fillna(21) # fill NaN challengeO2conc no values
    df_resp_allsujet = df_resp_allsujet.query(f"isGoodBreath == 1 and challengeO2conc == 21")

    #### sort 4 conds
    df_control_count = df_resp_allsujet.query(f"isControl == 1 and isChallenge == 0 ").groupby(['sujet', 'occlusionType']).size().reset_index(name='ctrl')
    df_chl_BOTH_count = df_resp_allsujet.query(f"isControl == 0 and isChallenge == 1 ").groupby(['sujet', 'occlusionType']).size().reset_index(name='chl_BOTH')
    df_chl_MECA_count = df_resp_allsujet.query(f"isControl == 0 and isChallenge == 1 and challengeLoadMag != 0 and challengeCO2conc == 0").groupby(['sujet', 'occlusionType']).size().reset_index(name='chl_MECA')
    df_chl_CO2_count = df_resp_allsujet.query(f"isControl == 0 and isChallenge == 1 and challengeLoadMag == 0 and challengeCO2conc != 0").groupby(['sujet', 'occlusionType']).size().reset_index(name='chl_CO2')
    
    df_plot = df_control_count
    keys = ['sujet', 'occlusionType']
    df_plot = df_plot.merge(df_chl_MECA_count, on=keys, how='outer')
    df_plot = df_plot.merge(df_chl_CO2_count, on=keys, how='outer')
    df_plot = df_plot.merge(df_chl_BOTH_count, on=keys, how='outer')

    df_plot = pd.melt(df_plot, id_vars=keys, value_vars=['ctrl', 'chl_BOTH', 'chl_MECA', 'chl_CO2'], value_name='count', var_name='cond')

    # occlusionType	:  0 no occlusion, 1 control occulusion, 2 occlusion
    df_plot = df_plot.query(f"occlusionType != 1")

    df_plot_allOCtogether = df_plot.groupby(['sujet', 'cond']).sum('count').reset_index()

    #### plot 
    path_export_plot = os.path.join(path_results, 'respi', 'count_cycles')

    g = sns.catplot(kind='bar', data=df_plot, x='sujet', y="count", hue='occlusionType', col='cond', sharey=False)
    for ax in g.axes.flat:
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
        ax.axhline(10, color='red', linestyle='-', linewidth=1)
    g.savefig(os.path.join(path_export_plot, f"count_cycle_allcond_allpatient_allOC.png"))
    # plt.show()

    plt.close('all')

    






#sujet = sujet_list[0]
def export_resp_mean(sujet):

    #### load
    os.chdir(os.path.join(path_precompute, 'RESP', 'session'))
    
    time_phase_list = ['pre', 'post']
    cycles_allcond = {}
    
    for cond in conditions:

        cycles_allcond[cond] = {}

        for time_phase in time_phase_list:

            cycles_allcond[cond][time_phase] = np.load(f"{sujet}_{cond}_stretch_resp_{time_phase}.npy")
    
    #### plot all cycles
    os.chdir(os.path.join(path_results, 'respi', 'session'))

    time_vec = np.arange(stretch_point_TF)

    for cond in conditions:

        for time_phase in time_phase_list: 
    
            data_plot = cycles_allcond[cond][time_phase]
            fig_resp_control, ax = plt.subplots()
            for cycle_i in range(data_plot.shape[0]):  
                plt.plot(time_vec, data_plot[cycle_i])
            plt.title(f"{sujet} {time_phase} / {cond} / n:{data_plot.shape[0]}")
            plt.tight_layout()
            # plt.show()
            fig_resp_control.savefig(f"{sujet}_{cond}_{time_phase}_allcycles.jpeg")

            plt.close('all')

    
    #### plot mean
    os.chdir(os.path.join(path_results, 'respi', 'session'))

    time_vec = np.arange(stretch_point_TF)

    label_line_cond_dict = {'rsp_ctrl' : {'linestyle' : '--', 'color' : 'r'}, 
                            'rsp_chl' : {'linestyle' : '-', 'color' : 'r'},
                            'oc_ctrl' : {'linestyle' : '--', 'color' : 'b'},
                            'oc_chl' : {'linestyle' : '-', 'color' : 'b'}
                            }

    for time_phase in time_phase_list: 

        fig_summary_mean, ax = plt.subplots()
        for cond in conditions:
            data_plot = cycles_allcond[cond][time_phase]
            plt.plot(time_vec, np.median(data_plot, axis=0), label=f"{cond} {data_plot.shape[0]}", linestyle=label_line_cond_dict[cond]['linestyle'], color=label_line_cond_dict[cond]['color'])
        plt.title(f"{sujet} {time_phase} median allcond")
        plt.legend()
        plt.tight_layout()
        # plt.show()
        fig_summary_mean.savefig(f"summary_{sujet}_{time_phase}.jpeg")

    plt.close('all')





    


def export_respi_allsujet():

    resp_data = np.zeros((len(sujet_list), len(conditions), stretch_point_TF))

    path_extract_data = os.path.join(path_precompute, 'RESP', 'session')

    for sujet_i, sujet in enumerate(sujet_list):

        _sujet_allcycles = []
                
        for cond_i, cond in enumerate(conditions):

            file_name = f"{sujet}_{cond}_stretch_resp_post.npy"
            _sujet_cond_allcycles = np.load(os.path.join(path_extract_data, file_name))
            _sujet_allcycles.append(_sujet_cond_allcycles)

        _sujet_allcycles = np.concat(_sujet_allcycles)
        _med = np.median(_sujet_allcycles)
        _mad = np.median(np.abs(_sujet_allcycles - _med))

        for cond_i, cond in enumerate(conditions):

            file_name = f"{sujet}_{cond}_stretch_resp_post.npy"
            _sujet_cond_allcycles = np.load(os.path.join(path_extract_data, file_name))
            _rzscore = (_sujet_cond_allcycles - _med) * 0.6745 / _mad
            resp_data[sujet_i,cond_i] = np.median(_rzscore, axis=0)

    #### plot
    time_vec = np.arange(stretch_point_TF)
    label_line_cond_dict = {'rsp_ctrl' : {'linestyle' : '--', 'color' : 'r'}, 
                            'rsp_chl' : {'linestyle' : '-', 'color' : 'r'},
                            'oc_ctrl' : {'linestyle' : '--', 'color' : 'b'},
                            'oc_chl' : {'linestyle' : '-', 'color' : 'b'}
                            }

    lmin_cond, lmax_cond = [], []

    for cond_i, cond in enumerate(conditions):

        lmin_cond.append(resp_data[:,cond_i].min())
        lmax_cond.append(resp_data[:,cond_i].max()) 

    lmin, lmax = np.array(lmin_cond).min(), np.array(lmax_cond).max()

    #### allsujet
    fig_cond_allsujet, axs = plt.subplots(ncols=len(conditions), figsize=(25,5))

    #cond = 'VS'
    for cond_i, cond in enumerate(conditions):

        ax = axs[cond_i]

        for sujet_i, sujet in enumerate(sujet_list):

            ax.plot(time_vec, resp_data[sujet_i,cond_i])

        ax.vlines(stretch_point_TF/2, ymin=lmin, ymax=lmax, color='r')
        ax.set_title(f'{cond}')

    plt.suptitle('allsujet')

    # fig_cond_allsujet.show()

    os.chdir(os.path.join(path_results, 'respi', 'plot_median'))
    fig_cond_allsujet.savefig(f"allsujet_respi_median.png")

    fig_cond_allsujet.savefig(os.path.join(path_paper_figure_export, f"fig03c_allsujet_respi_median.svg"))

    plt.close('all')

    #### median
    ymin, ymax = [], []

    for cond_i, cond in enumerate(conditions):

        ymin.append(np.median(resp_data[:,cond_i], axis=0).min())
        ymax.append(np.median(resp_data[:,cond_i], axis=0).max())

    ymin, ymax = np.array(ymin).min(), np.array(ymax).max()

    fig_median, ax = plt.subplots()

    #cond = 'VS'
    for cond_i, cond in enumerate(conditions):

        ax.plot(time_vec, np.median(resp_data[:,cond_i], axis=0), color=label_line_cond_dict[cond]['color'], linestyle=label_line_cond_dict[cond]['linestyle'], label=cond)

    ax.vlines(stretch_point_TF/2, ymin=ymin, ymax=ymax, color='r')
    plt.ylim(ymin, ymax)
    plt.title('median')
    plt.legend()

    # fig_median.show()

    os.chdir(os.path.join(path_results, 'respi', 'plot_median'))
    fig_median.savefig(f"respi_median.png")

    fig_median.savefig(os.path.join(path_paper_figure_export, f"fig03b_respi_median.svg"))

    
    






################################
######## EXECUTE ########
################################


if __name__ == '__main__':

    export_allcycles_fig_for_allcond()

    for sujet in sujet_list:

        export_resp_mean(sujet)

    export_respi_allsujet()
    



                        