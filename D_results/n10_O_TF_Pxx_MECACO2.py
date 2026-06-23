



from A_config.n01_O_params import *
from A_config.n02_O_analysis_functions import *
from A_config.n03_X_manip_data import *






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





    



################################
######## EXECUTE ########
################################


if __name__ == '__main__':

    Pxx_relative_patientwise_allpatient()

    
    






















    


    


















