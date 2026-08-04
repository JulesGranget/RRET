



from A_config.n01_O_params import *
from A_config.n02_O_analysis_functions import *
from A_config.n03_X_manip_data import *

import joblib




################################
######## PXX RELATIVE ########
################################

def Pxx_relative_patientwise_allpatient():

    #### load Pxx
    path_load_data = os.path.join(path_precompute, 'TF', 'MECACO2_INTER')
    
    df_Pxx = []

    for sujet in sujet_list_interaction_analysis:
        path_openfile = os.path.join(path_load_data, f"{sujet}_xr_Pxx_MECACO2INTER.nc")
        da = xr.load_dataarray(path_openfile)
        da = da.expand_dims(sujet=[sujet])
        _roi_list = da['ROI'].values

        mask_ROI = np.array([True if _ROI.find('UNSORTED') == -1 and _ROI not in ['Inf-Lat-Vent', 'WM', 'unknown'] else False for _ROI in _roi_list ])

        da = da.isel(chan=mask_ROI)
        da = da.dropna(dim="cycle")

        df_Pxx.append(da.to_dataframe().reset_index())

    df_Pxx = pd.concat(df_Pxx)

    #### plot
    path_export_plot = os.path.join(path_results, 'Pxx', 'MECACO2_BOTH')
    df_diff = pd.pivot_table(df_Pxx, values='Pxx', columns='cond', index=[col for col in df_Pxx.columns if col not in ['Pxx', 'cond']]).reset_index()
    df_diff['MECA_diff'] = df_diff['MECA'] - df_diff['ctrl']
    df_diff['CO2_diff'] = df_diff['CO2'] - df_diff['ctrl']
    df_diff['BOTH_MC_diff'] = df_diff['BOTH_MC'] - df_diff['ctrl']
    df_diff = df_diff.drop(columns=['MECA', 'CO2', 'BOTH_MC', 'ctrl'])
    col_names_diff = ['MECA_diff', 'CO2_diff', 'BOTH_MC_diff']
    df_diff = pd.melt(df_diff, id_vars=[col for col in df_diff.columns if col not in col_names_diff], value_vars=col_names_diff, value_name='Pxx', col_level='cond')

    for band in df_diff['band'].unique():
    
        df_plot = df_diff.query(f"band == '{band}' and ROI in {ROI_short_list}")
        sns.catplot(df_plot, kind='strip', x='cond', y='Pxx', hue='sujet', row='phase', col='ROI', sharey=False, size=3)
        # plt.show()

        plt.savefig(os.path.join(path_export_plot, f"{band}_allpatient_allphase_allcond.jpg"))






def tf_absolute_allpatient():

    #### load
    os.chdir(os.path.join(path_results, 'anatomy'))
    df_loca_allpatient = pd.read_excel('df_loca_allsujet_removetag.xlsx').query(f"SOZ != 1 and Spike != 1 and Out != 1 and NotFound != 1 and BAD != 1 and sujet in {sujet_list_interaction_analysis}")
    localist_unique = ROI_short_list_MECACO2
    oc_cond_list = ['noc', 'oc']

    ticks_freq = [2, 8, 12, 30, 50, 60, 150]
    percentile_plot = [1,99]
    time_vec = np.arange(stretch_point_TF)

    #loca_sel_i, loca_sel = 0, ROI_short_list_MECACO2[0]
    # for loca_sel_i, loca_sel in enumerate(ROI_short_list_MECACO2):
    def plot_allsujet_loca(loca_sel_i, loca_sel):

        print_advancement(loca_sel_i, len(localist_unique), [25,50,75])

        #### load sujet
        os.chdir(os.path.join(path_precompute, 'TF', 'MECACO2_INTER'))
        loca_sel_vec = df_loca_allpatient['loca'].values == loca_sel
        _patient_list_unique = np.unique(df_loca_allpatient['sujet'].values[loca_sel_vec])
        _tf_allsujet = np.zeros((_patient_list_unique.size, len(oc_cond_list), len(cond_list_interaction), frex.size, stretch_point_TF))
        _ncycle_tot = np.zeros((len(oc_cond_list), len(cond_list_interaction)))

        for _sujet_i, _sujet in enumerate(_patient_list_unique):

            _sujet_chanlist, _sujet_localist = get_chanlist(_sujet)
            _sujet_chan_sel_vec = _sujet_localist == loca_sel
            
            for oc_cond_i, oc_cond in enumerate(oc_cond_list):
                
                for cond_i, cond in enumerate(cond_list_interaction):

                    _tf = np.load(f'{_sujet}_{oc_cond}_{cond}_tf_allchan_stretch_MECACO2INTER.npy')[_sujet_chan_sel_vec]
                    _tf_allsujet[_sujet_i, oc_cond_i, cond_i, :, :] = np.median(np.median(_tf, axis=1), axis=0)
                    _ncycle_tot[oc_cond_i, cond_i] += _tf.shape[1]

        _tf_allsujet_median = np.median(_tf_allsujet, axis=0)

        #### DIV
            #### plot TF
        _vlim = np.abs([np.percentile(np.stack(_tf_allsujet_median, axis=0), percentile_plot[0]), np.percentile(np.stack(_tf_allsujet_median, axis=0), percentile_plot[1])]).max()

        fig, axs = plt.subplots(ncols=4, nrows=2, figsize=(16,8))

        for oc_cond_sel_i, oc_cond_sel in enumerate(oc_cond_list):
        
            for cond_i, cond in enumerate(cond_list_interaction):
            
                ax = axs[oc_cond_sel_i, cond_i]
                pcm = ax.pcolormesh(time_vec, frex, _tf_allsujet_median[oc_cond_sel_i, cond_i], vmin=-_vlim, vmax=_vlim, cmap='seismic')
                ax.vlines(x=int(stretch_point_TF/2), ymin=frex[0], ymax=frex[-1], color="k", linewidth=2)
                ax.set_yscale('log')
                
                ax.set_yticks(ticks_freq)
                ax.set_yticklabels([str(t) for t in ticks_freq])

                if cond_i == 0:
                    ax.set_ylabel(oc_cond_sel)
                if oc_cond_sel_i == 1:
                    ax.set_xlabel('Phase')
                if oc_cond_sel_i == 0:
                    ax.set_title(f"{cond} c{int(_ncycle_tot[oc_cond_sel_i, cond_i])}")

        cbar = fig.colorbar(pcm, ax=axs, orientation="vertical", fraction=0.02, pad=0.04)
        plt.suptitle(f"{loca_sel} s({len(_patient_list_unique)})")

        # plt.show()

            #### save
        os.chdir(os.path.join(path_results, 'TF', 'MECACO2', 'absolute'))
        fig.savefig(f"DIV_allsujet_{loca_sel}_TF.jpg")
        plt.close('all')

            #### plot lateral median
        med_lines = np.stack([np.median(_tf_allsujet_median[:,:,:,:int(stretch_point_TF/2)], axis=-1), np.median(_tf_allsujet_median[:,:,:,int(stretch_point_TF/2):], axis=-1)])
        _vlim = np.abs([np.min(med_lines), np.max(med_lines)]).max()

        for phase_i, phase in enumerate(['inspi', 'expi']):

            for oc_cond_sel_i, oc_cond_sel in enumerate(oc_cond_list):
                    
                for cond_i, cond in enumerate(conditions):

                    med_lines[phase_i, oc_cond_sel_i, cond_i] = scipy.signal.savgol_filter(med_lines[phase_i, oc_cond_sel_i, cond_i], window_length=10, polyorder=4)

        fig, axs = plt.subplots(ncols=4, nrows=2, figsize=(20,8))

        for oc_cond_sel_i, oc_cond_sel in enumerate(oc_cond_list):
        
            for cond_i, cond in enumerate(conditions):
            
                ax = axs[oc_cond_sel_i, cond_i]
                ax.plot(frex, med_lines[0, oc_cond_sel_i, cond_i], color='b')
                ax.plot(frex, med_lines[1, oc_cond_sel_i, cond_i], color='r')
                ax.set_xscale('log')
                ax.set_ylim(-_vlim, _vlim)
                
                ax.set_xticks(ticks_freq)
                ax.set_xticklabels([str(t) for t in ticks_freq], rotation=45, ha="right")

                if cond_i == 0:
                    ax.set_xlabel(oc_cond_sel)
                if oc_cond_sel_i == 1:
                    ax.set_ylabel('Power Median')
                if oc_cond_sel_i == 0:
                    ax.set_title(f"{cond} c{int(_ncycle_tot[oc_cond_sel_i, cond_i])}")

        plt.suptitle(f"{loca_sel} s({len(_patient_list_unique)})")

        # plt.show()

            #### save
        os.chdir(os.path.join(path_results, 'TF', 'MECACO2', 'absolute'))
        fig.savefig(f"MEDLINE_allsujet_{loca_sel}.jpg")
        plt.close('all')

        #### JET
            #### plot TF
        _vlim = np.abs([np.percentile(np.stack(_tf_allsujet_median, axis=0), percentile_plot[0]), np.percentile(np.stack(_tf_allsujet_median, axis=0), percentile_plot[1])]).max()
        
        fig, axs = plt.subplots(ncols=4, nrows=2, figsize=(16,8))

        for oc_cond_sel_i, oc_cond_sel in enumerate(oc_cond_list):
        
            for cond_i, cond in enumerate(cond_list_interaction):
            
                ax = axs[oc_cond_sel_i, cond_i]
                pcm = ax.pcolormesh(time_vec, frex, _tf_allsujet_median[oc_cond_sel_i, cond_i], vmin=-_vlim, vmax=_vlim, cmap='jet')
                ax.vlines(x=int(stretch_point_TF/2), ymin=frex[0], ymax=frex[-1], color="k", linewidth=2)
                ax.set_yscale('log')
                
                ax.set_yticks(ticks_freq)
                ax.set_yticklabels([str(t) for t in ticks_freq])

                if cond_i == 0:
                    ax.set_ylabel(oc_cond_sel)
                if oc_cond_sel_i == 1:
                    ax.set_xlabel('Phase')
                if oc_cond_sel_i == 0:
                    ax.set_title(f"{cond} c{int(_ncycle_tot[oc_cond_sel_i, cond_i])}")

        cbar = fig.colorbar(pcm, ax=axs, orientation="vertical", fraction=0.02, pad=0.04)
        plt.suptitle(f"{loca_sel} s({len(_patient_list_unique)})")

        # plt.show()

            #### save
        os.chdir(os.path.join(path_results, 'TF', 'MECACO2', 'absolute'))
        fig.savefig(f"JET_allsujet_{loca_sel}_TF.jpg")
        plt.close('all')

    joblib.Parallel(n_jobs = n_core, prefer = 'processes')(joblib.delayed(plot_allsujet_loca)(loca_sel_i, loca_sel) for loca_sel_i, loca_sel in enumerate(localist_unique))
    
    plt.close('all')



    



################################
######## EXECUTE ########
################################


if __name__ == '__main__':

    Pxx_relative_patientwise_allpatient()

    tf_absolute_allpatient()

    
    






















    


    


















