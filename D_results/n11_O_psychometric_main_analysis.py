



from A_config.n01_O_params import *
from A_config.n02_O_analysis_functions import *
from A_config.n03_X_manip_data import *






################################
######## EXPORT FIG ########
################################

def export_fig_psychometric():

    #### load
    os.chdir(os.path.join(path_precompute, 'PSYCHO', 'df_export')) 
    df_psycho = pd.concat([pd.read_excel(f"{sujet}_df_psycho.xlsx") for sujet in sujet_list])[['sujet', 'cond', 'trial', 'trialUnpleasantness', 'trialAnxiety', 'trialAnxietySource']]
    df_psycho = df_psycho.groupby(['cond', 'sujet', 'trial']).median().reset_index()
    df_psycho = df_psycho.query(f"cond in ['rsp_ctrl', 'rsp_chl']")
    df_psycho_tidy = pd.melt(df_psycho, id_vars=['cond', 'sujet', 'trial'], value_vars=['trialUnpleasantness', 'trialAnxiety', 'trialAnxietySource'], var_name='psychometric', value_name='val')

    df_psycho_2states = df_psycho.copy()
    df_psycho_2states['resp'] = df_psycho_2states['cond'].str.split('_').str[1]
    df_psycho_2states = df_psycho_2states.drop(columns='cond')
    df_psycho_2states = pd.melt(df_psycho_2states, id_vars=['sujet', 'trial', 'resp'], value_vars=['trialUnpleasantness', 'trialAnxiety', 'trialAnxietySource'], var_name='psychometric', value_name='val')

    #### plot
    os.chdir(os.path.join(path_results, 'Psycho')) 

    df_plot = df_psycho_tidy
    g = sns.catplot(data=df_plot, kind='strip', x='sujet', y='val', hue='cond', col='psychometric', dodge=True)
    g.fig.suptitle("alltrial")
    g.fig.subplots_adjust(top=0.90)
    # plt.show()
    g.fig.savefig("psycho_alltrial.png")

    df_plot = df_psycho_tidy
    g = sns.catplot(data=df_plot, kind='bar', dodge=True, x='cond', y='val', col='psychometric')
    g.fig.suptitle("median_trial_allsujet")
    g.fig.subplots_adjust(top=0.90)
    # plt.show()
    g.fig.savefig("psycho_median_trial_allsujet.png")

    df_plot = df_psycho_2states
    g = sns.catplot(data=df_plot, kind='bar', dodge=True, x='sujet', y='val', hue='resp', col='psychometric', hue_order=['ctrl', 'chl'])
    g.fig.suptitle("median trial")
    g.fig.subplots_adjust(top=0.90)
    # plt.show()
    g.fig.savefig("psycho_median_trial.png")

    df_plot = df_psycho_2states.query(f"psychometric == 'trialUnpleasantness'").drop(columns=['psychometric', 'val']).groupby(['sujet', 'resp']).count()
    g = sns.catplot(data=df_plot, kind='bar', dodge=True, x='sujet', y='trial', hue='resp', hue_order=['ctrl', 'chl'])
    g.fig.suptitle("trial number")
    g.fig.subplots_adjust(top=0.90)
    # plt.show()
    g.fig.savefig("psycho_trial_number.png")


    



    



################################
######## EXECUTE ########
################################


if __name__ == '__main__':

    export_fig_psychometric()

    
    






















    


    


















