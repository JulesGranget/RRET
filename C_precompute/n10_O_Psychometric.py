


from A_config.n01_O_params import *
from A_config.n02_O_analysis_functions import *
from A_config.n03_X_manip_data import *







################################
######## EXTRACT STRESS ########
################################

def get_df_unpl():

    #### config
    col_to_extract = ['sujet', 'cond', 'trialUnpleasantness', 'trialAnxiety', 'trialAnxietySource']
    path_import_df_psycho = os.path.join(path_precompute, 'PSYCHO', 'df_export')
    path_export_res = os.path.join(path_results, 'Psycho')

    #### load_df
    df_psycho = pd.concat([pd.read_excel(os.path.join(path_import_df_psycho, f'{sujet}_df_psycho.xlsx')) for sujet in sujet_list_interaction_analysis])

    #### trial num
    df_plot = df_psycho[["sujet", "cond", "trial"]].drop_duplicates().reset_index(drop=True).groupby(['sujet', 'cond']).count().reset_index()

    sns.catplot(data=df_plot, kind='bar', x='sujet', y='trial', hue='cond', hue_order=['ctrl', 'MECA', 'CO2', 'BOTH'])
    # plt.show()
    plt.savefig(os.path.join(path_export_res, f"trial_num.png"))

    #### plot
    df_plot = df_psycho.groupby(['sujet', 'cond', 'trial']).median().reset_index()[col_to_extract]
    df_plot = pd.melt(df_psycho, id_vars=['sujet', 'cond'], value_vars=['trialUnpleasantness', 'trialAnxiety', 'trialAnxietySource'], var_name='psychometric', value_name='val')
        
    sns.catplot(data=df_plot, kind='bar', x='sujet', y='val', hue='cond', col='psychometric', hue_order=['ctrl', 'MECA', 'CO2', 'BOTH'])
    # plt.show()
    plt.savefig(os.path.join(path_export_res, f"psychometric_values.png"))
        
    plt.close('all')        










################################
######## EXECUTE ########
################################


if __name__ == '__main__':
    
    get_df_unpl()

    


                        