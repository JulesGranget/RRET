

from A_config.n01_O_params import *
from A_config.n02_O_analysis_functions import *
from A_config.n03_X_manip_data import *






################################
######## EXPORT DF ########
################################

def export_Pxx_df_psycho_order():

    path_psycho_file = os.path.join(path_precompute, 'PSYCHO', 'df_export')
    path_export_df_R = os.path.join(path_precompute, 'PSYCHO', 'df_R')

    df_allsujet_psycho_order = []
        
    for sujet in sujet_list:
            
        _df = pd.read_excel(os.path.join(path_psycho_file, f"{sujet}_df_psycho.xlsx"))
        
        _df = _df[['sujet', 'cond', 'trial', 'trialUnpleasantness', 'trialAnxiety']]

        _df_trial = _df.groupby(["sujet", "cond", "trial"], as_index=False).first().query(f"cond in ['rsp_ctrl', 'rsp_chl']")

        trial_chl_list_i = []
        trial_ctrl_list_i = []

        for trial_i in _df_trial['trial'].unique():

            _cond_extracted = _df_trial.query(f"trial == {trial_i}")['cond'].values

            if 'rsp_chl' not in _cond_extracted:
                trial_ctrl_list_i.append(int(trial_i))
            else:
                trial_chl_list_i.append(int(trial_i))

        _df_trial = pd.concat([_df_trial.query(f"cond == 'rsp_ctrl' and trial in {trial_ctrl_list_i}"), _df_trial.query(f"cond == 'rsp_chl' and trial in {trial_chl_list_i}")])

        df_allsujet_psycho_order.append(_df_trial)

    df_allsujet_psycho_order = pd.concat(df_allsujet_psycho_order)

    #### export
    df_allsujet_psycho_order.to_excel(os.path.join(path_export_df_R, f"df_PSYCHO_ORDER_R.xlsx"))

    






################################
######## EXECUTE ########
################################


if __name__ == '__main__':

    export_Pxx_df_psycho_order()



                        