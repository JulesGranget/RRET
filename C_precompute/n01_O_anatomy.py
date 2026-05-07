



from A_config.n01_O_params import *
from A_config.n02_O_analysis_functions import *
from A_config.n03_X_manip_data import *






def export_anat_allpatient_df():

    #### import alldata
    filename = os.path.join(path_results, 'anatomy', "anatomy_allpatient.xlsx")
    df_allpatient = pd.read_excel(filename)

    #### generate remove tag df
    df_loca_remove_allsujet = []

    #sujet = sujet_list[0]
    for sujet in sujet_list:

        df_loca_sujet = df_allpatient.query(f"sujet == '{sujet}'")

        filepath_remove_elec = os.path.join(path_general, 'Data', 'anatomy', f"{sujet}_elec_recon", 'elec_recon', f"{sujet}_Electrodes_Natus_TDT_correspondence.xlsx")
        df_remove_chan_sujet = pd.read_excel(filepath_remove_elec)

        df_loca_sel_sujet = []
        for _row_i, _row_val in df_loca_sujet.iterrows():
            _chan, _loca = _row_val['chan_name'], _row_val['loca']
            if _loca.find('UNSORTED') == -1 or _loca != 'WM' or _loca != 'unknown':
                _df_chan = df_remove_chan_sujet.query(f"Label == '{_chan}'")
                if _df_chan.shape[0] == 0:
                    _df = pd.DataFrame({'sujet' : [sujet], 'chan' : [_chan], 'loca' : [_loca], 'SOZ' : [0],
                    'Spike' : [0], 
                    'Out' : [0], 
                    'NotFound' : [1], 
                    'BAD' : [0]})
                else:
                    _df = pd.DataFrame({'sujet' : [sujet], 'chan' : [_chan], 'loca' : [_loca], 'SOZ' : [0 if np.isnan(_val) else _val for _val in _df_chan['SOZ (1, 0/empty)'].values],
                                        'Spike' : [0 if np.isnan(_val) else _val for _val in _df_chan['Spikey (1, 0/empty)'].values], 
                                        'Out' : [0 if np.isnan(_val) else _val for _val in _df_chan['Out (1, 0/empty)'].values], 
                                        'NotFound' : [0], 
                                        'BAD' : [0 if np.isnan(_val) else _val for _val in _df_chan['BAD (1, 0/empty)'].values]})
            df_loca_sel_sujet.append(_df)

        df_loca_sel_sujet = pd.concat(df_loca_sel_sujet)

        df_loca_remove_allsujet.append(df_loca_sel_sujet)     

    df_loca_remove_allsujet = pd.concat(df_loca_remove_allsujet)
    df_loca_remove_allsujet = df_loca_remove_allsujet.reset_index().drop(columns=['index'])

    #### export df
    filepath_export_df_remove = os.path.join(path_results, 'anatomy', 'df_loca_allsujet_removetag.xlsx')
    df_loca_remove_allsujet.to_excel(filepath_export_df_remove)

    



################################
######## EXECUTE ########
################################


if __name__ == '__main__':

    export_anat_allpatient_df()



                        