



from n00_config_params import *
from n00bis_config_analysis_functions import *
from n01_manip_data import *







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

    




################################
######## EXECUTE ########
################################


if __name__ == '__main__':

    
    for sujet in sujet_list:

        export_resp_mean(sujet)



                        