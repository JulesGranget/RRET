



from n00_config_O_params import *
from n00bis_config_O_analysis_functions import *
from n00ter_X_manip_data import *







########################################
######## PLOT RESPI MEAN ########
########################################

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

    
    




    


def export_count_cycle_allsujet():

    #### load data
    input_dir = os.path.join(path_results, 'respi', 'count_cycles')
    filename = os.path.join(input_dir, "ALLSUJET_count_cycle.xlsx")
    df_count_allsujet = pd.read_excel(filename)

    df_count_allsujet = df_count_allsujet.drop(columns=['Unnamed: 0'])
    df_count_allsujet = df_count_allsujet.rename(columns={'resp_control' : 'rsp_ctrl', 'resp_chall' : 'rsp_chl', 'oc_control' : 'oc_ctrl', 'oc_chall' : 'oc_chl'})
    df_count_allsujet = df_count_allsujet.query(f"sujet in {sujet_list}")
    df_count_allsujet.melt()

    df_count_allsujet = df_count_allsujet.melt(id_vars="sujet", value_vars=["rsp_ctrl", "rsp_chl", "oc_ctrl", "oc_chl"],
                var_name="cond", value_name="count")
    
    #### plot
    fig, ax = plt.subplots()
    sns.barplot(df_count_allsujet, x='sujet', y='count', hue='cond', ax=ax)

    # plt.show()

    output_dir = os.path.join(path_results, 'respi', 'count_cycles')
    filename = os.path.join(output_dir, "ALLSUJET_count_cycle.png")
    fig.savefig(filename)



################################
######## EXECUTE ########
################################


if __name__ == '__main__':

    
    for sujet in sujet_list:

        export_resp_mean(sujet)

    export_respi_allsujet()
    export_count_cycle_allsujet()



                        