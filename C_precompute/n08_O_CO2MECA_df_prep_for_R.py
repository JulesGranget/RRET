

import joblib


from A_config.n01_O_params import *
from A_config.n02_O_analysis_functions import *
from A_config.n03_X_manip_data import *

import joblib
import statsmodels.api as sm
from io import StringIO
from statsmodels.tools.sm_exceptions import ValueWarning
import warnings





################################
######## EXPORT DF ########
################################

def export_Pxx_df():

    df_allband = []

    #band = 'gamma'    
    for band in freq_band_dict:

        print(band)
        
        #### load
        os.chdir(os.path.join(path_precompute, 'TF', 'MECACO2_INTER'))

        df_band_Pxx_allsujet = []

        for sujet in sujet_list_interaction_analysis:

            _xr = xr.open_dataarray(f"{sujet}_xr_Pxx_MECACO2INTER.nc").sel(band=band)
            df_band_Pxx_allsujet.append(_xr.to_dataframe().reset_index().dropna())

        df_band_Pxx_allsujet = pd.concat(df_band_Pxx_allsujet)
        df_band_Pxx_allsujet = df_band_Pxx_allsujet.query(f"ROI in {ROI_short_list}")

        df_allband.append(df_band_Pxx_allsujet)

    df_allband = pd.concat(df_allband)

    #### export
    os.chdir(os.path.join(path_precompute, 'TF', 'session', 'df_R'))
    df_allband.to_excel(f"df_R_allband_MECACO2INTER.xlsx")

    #### plot
    band = 'gamma'

    #### plot general
    df_band_Pxx_allsujet_plot = df_allband.groupby(['ROI', 'oc_cond', 'cond', 'phase', 'band']).median('Pxx').reset_index(drop=False)
    df_band_Pxx_allsujet_plot = df_band_Pxx_allsujet_plot.drop(columns=['cycle'])

    df_plot = df_band_Pxx_allsujet_plot
    sns.catplot(data=df_plot, kind='bar', x='cond', y='Pxx', hue='oc_cond', row='phase', col='ROI', sharey=False, 
                order=['ctrl', 'MECA', 'CO2', 'BOTH'], row_order=['inspi', 'expi'])
    plt.show()

    #### plot diff
    df_plot_diff = df_band_Pxx_allsujet_plot.pivot(index=[col for col in df_band_Pxx_allsujet_plot.columns if col not in ['oc_cond', 'Pxx']], columns='oc_cond', values='Pxx').reset_index(drop=False)

    df_plot_diff[f"Pxx"] = df_plot_diff['oc'] - df_plot_diff['noc']

    df_plot_diff = df_plot_diff[['ROI', 'cond', 'phase', 'band', 'Pxx']]

    sns.catplot(data=df_plot_diff, kind='bar', x='cond', y='Pxx', hue='phase',
                order=['ctrl', 'MECA', 'CO2', 'BOTH'], hue_order=['inspi', 'expi'], col='ROI', sharey=False)
    plt.show()


    






################################
######## REG ########
################################



def get_reg_allsujet():

    if os.path.exists(os.path.join(path_precompute, 'TF', 'session', 'df_reg', f"df_reg_MECACO2.xlsx")):
        
        print('ALREADY COMPUTED')
    
    #### load Pxx
    das = []
    chan_list_shape = []

    for sujet in sujet_list_interaction_analysis:
        
        fp = os.path.join(path_precompute, "TF", "MECACO2_INTER", f"{sujet}_xr_Pxx_MECACO2INTER.nc")
        da = xr.load_dataarray(fp)

        da = da.expand_dims(sujet=[sujet])

        chan_list_shape.append(da.shape[3])
        das.append(da)

    max_nchan = np.max(np.array(chan_list_shape))

    das_padded = []
    for da in das:

        nchan = da.sizes["chan"]
        
        _da_padded = da.assign_coords(chan=np.arange(nchan)).reindex(chan=np.arange(max_nchan)).assign_coords(
            chan_label=("chan", list(da["chan"].values) + [np.nan] * (max_nchan - nchan)),
            ROI=("chan", list(da["ROI"].values) + [np.nan] * (max_nchan - nchan)))

        das_padded.append(_da_padded)

    da_Pxx_allsujet_phase_cycle = xr.concat(das_padded, dim="sujet")
    da_Pxx_allsujet_phase_cycle = da_Pxx_allsujet_phase_cycle.rename({'phase' : 'phase_cycle'})
    da_whole = da_Pxx_allsujet_phase_cycle.median('phase_cycle').expand_dims(phase_cycle=["whole"])
    da_Pxx_allphase_cycle = xr.concat([da_Pxx_allsujet_phase_cycle, da_whole],dim="phase_cycle")
    
    #### load resp features
    phase_cycle_list = ['inspi', 'expi']
    rf_metric_allphase_cycle = ['oc_ratio', 'oc_val']
    
    resp_vars = {}

    for sujet in sujet_list_interaction_analysis:

        # load respiration features
        os.chdir(os.path.join(path_precompute, 'RESP', 'respfeatures')) 
        _oc_respfeatures = pd.read_excel(f'{sujet}_MECACO2_df_oc.xlsx')

        _cond_da = {}
        
        for cond_i, cond in enumerate(cond_list_interaction):

            resp_export = _oc_respfeatures.query(f"cond == '{cond}'")[['oc_ratio', 'oc_val']]
            _vals = resp_export.values

            ds_resp = xr.Dataset(
                data_vars={
                    "respfeatures": (("cycle", "feature"),
                                    _vals)
                },
                coords={
                    "cycle": np.arange(_vals.shape[0]),
                    "feature": rf_metric_allphase_cycle,
                    "sujet": sujet,
                    "cond": cond,
                }
            ).expand_dims(("sujet", "cond"))

            # collect
            _cond_da[cond] = ds_resp

        resp_vars[sujet] = xr.concat(list(_cond_da.values()), dim="cond").sortby("cond")

    ds_respfeatures_all = xr.concat(list(resp_vars.values()), dim="sujet").sortby("sujet")

    #### merge
    ds_whole_oc = xr.merge([da_Pxx_allphase_cycle, ds_respfeatures_all]).sel(oc_cond='oc')

    #### df
    df_whole = ds_whole_oc.to_dataframe()
    df_whole = df_whole[~df_whole.isna()['Pxx'].values]
    df_whole = df_whole.query(f"ROI in {ROI_short_list_MECACO2} and phase_cycle != 'whole'")
    df_whole = df_whole.reset_index(drop=False).drop(columns='oc_cond')
    df_whole = df_whole.pivot(index=[col for col in df_whole.columns if col not in ['feature', 'respfeatures']], columns='feature', values='respfeatures').reset_index(drop=False)

    filename = os.path.join(path_precompute, 'TF', 'session', 'df_R', f"df_reg_MECACO2_R.xlsx")
    df_whole.to_excel(filename)









################################
######## EXECUTE ########
################################


if __name__ == '__main__':

    export_Pxx_df()
    get_reg_allsujet()




                        