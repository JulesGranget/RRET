



from n00_config_O_params import *
from n00bis_config_O_analysis_functions import *
from n00ter_X_manip_data import *



####################################
######## EXTRACT OC SIZE ########
####################################

#sujet = sujet_list[0]
def extract_oc_size(sujet):

    print(sujet)

    #### load
    os.chdir(os.path.join(path_precompute, 'RESP', 'respfeatures')) 
    respfeatures = pd.read_excel(f'{sujet}_respfeatures_cleaned_label.xlsx')

    #### extract
    df_oc = []

    #cond = conditions[2]
    for cond in conditions:

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

                plt.plot(resp_stretch[cycle_i])
            
            plt.show()

            for cycle_i in range(resp_stretch.shape[0]):

                plt.plot(np.diff(resp_stretch[cycle_i]))
            
            plt.show()

            cycle_i = 0

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
                plt.show()

        _df_oc = pd.DataFrame({'sujet' : [sujet]*len(oc_ratio), 'cond' : [cond]*len(oc_ratio), 'cycle_i' : np.arange(len(oc_ratio)).tolist(), 'oc_ratio' : oc_ratio, 'oc_val' : oc_val})

        if _respfeature.shape[0] != _df_oc.shape[0]:

            raise ValueError('!!! NOT SAME CYCLE NUMBER !!!')

        df_oc.append(_df_oc)

    df_oc = pd.concat(df_oc)

    #### plot
    df_plot = pd.melt(df_oc, id_vars=[_col for _col in df_oc.columns if _col not in ['oc_ratio', 'oc_val']], value_vars=['oc_ratio', 'oc_val'], var_name="metric_type", value_name="val",)
    g = sns.catplot(df_plot, kind='strip', x='cond', y='val', col='metric_type', jitter=True, sharey=False)
    plt.suptitle(f"{sujet}")
    # plt.show()


    #### save
    os.chdir(os.path.join(path_results, 'respi', 'oc_ratio'))
    g.savefig(f"{sujet}_oc.png")
    plt.close('all')

    os.chdir(os.path.join(path_precompute, 'RESP', 'respfeatures')) 
    df_oc.to_excel(f"{sujet}_df_oc.xlsx")





################################
######## EXECUTE ########
################################


if __name__ == '__main__':

    
    for sujet in sujet_list:

        extract_oc_size(sujet)



                        