



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
    df_psycho_med_patientwise = df_psycho_2states.groupby(['sujet', 'resp', 'psychometric']).median('trial').drop(columns='trial').reset_index()

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


    #### stats
    psychometric_list = ['trialAnxiety', 'trialAnxietySource', 'trialUnpleasantness']

    stat_dict = {}

    for psychometric in psychometric_list:

        _df_stats = df_psycho_med_patientwise.query(f"psychometric == '{psychometric}'")
        stat_val, pval = scipy.stats.wilcoxon(x=_df_stats.query(f"resp == 'chl'")['val'].values, y=_df_stats.query(f"resp == 'ctrl'")['val'].values, alternative='two-sided')

        # _df_stats = df_psycho_tidy.query(f"psychometric == '{psychometric}'")
        # stat_val, pval = scipy.stats.mannwhitneyu(x=_df_stats.query(f"cond == 'rsp_chl'")['val'].values, y=_df_stats.query(f"cond == 'rsp_ctrl'")['val'].values, alternative='two-sided')

        stat_dict[psychometric] = {'stat_val' : stat_val, 'pval' : pval}

    #### patient linked
    df_plot = df_psycho_med_patientwise

    cond_order = ['ctrl', 'chl']
    psychometric_order = ['trialUnpleasantness', 'trialAnxiety', 'trialAnxietySource']

    df_plot['resp'] = pd.Categorical(df_plot['resp'], categories=cond_order, ordered=True)
    df_plot['psychometric'] = pd.Categorical(df_plot['psychometric'], categories=psychometric_order, ordered=True)

    fig, ax = plt.subplots(figsize=(7, 6))

    palette = {'ctrl': 'tab:blue', 'chl': 'tab:orange'}

    sns.barplot(data=df_plot, x='psychometric', y='val', hue='resp', order=psychometric_order, hue_order=cond_order,
        palette=palette, errorbar=None, alpha=0.4, dodge=True, ax=ax)

    psycho_x = {phase: i for i, phase in enumerate(psychometric_order)}

    cond_offset = {'ctrl': -0.2, 'chl': 0.2}

    for (sujet, phase), df_sub in df_plot.groupby(['sujet', 'psychometric'],observed=True):

        df_sub = df_sub.sort_values('resp')

        if df_sub['resp'].nunique() != len(cond_order):
            continue

        x = [psycho_x[phase] + cond_offset[cond] for cond in df_sub['resp']]

        y = df_sub['val'].to_numpy()

        ax.plot(x, y, color='gray', alpha=0.45, linewidth=1, zorder=2)

    for cond in cond_order:

        df_cond = df_plot[df_plot['resp'] == cond]

        x = [psycho_x[phase] + cond_offset[cond] for phase in df_cond['psychometric']]

        ax.scatter(x, df_cond['val'], color=palette[cond], edgecolor='black', linewidth=0.4, s=45, zorder=3)

    ax.set_xlabel('cond')
    ax.set_ylabel('val')
    ax.legend(title='Condition')

    plt.suptitle(f"psychometric")
    plt.tight_layout()

    # plt.show()

    os.chdir(os.path.join(path_results, 'Psycho')) 
    fig.savefig(f"psycho_patientlinked.png")

    ########

    cond_order = ["ctrl", "chl"]

    palette = {
        "ctrl": "tab:blue",
        "chl": "tab:orange"
    }

    # Create median bars, automatically faceted by psychometric
    g = sns.catplot(
        data=df_plot,
        kind="bar",
        x="resp",
        y="val",
        hue="resp",
        col="psychometric",
        order=cond_order,
        estimator=np.median,
        errorbar=None,
        palette=palette,
        alpha=0.4,
        legend=False,
        height=5,
        aspect=0.5
    )

    # Add one paired line per subject to every panel
    for psychometric, ax in zip(g.col_names, g.axes.flat):

        sns.pointplot(
            data=df_plot.query("psychometric == @psychometric"),
            x="resp",
            y="val",
            hue="sujet",
            order=cond_order,
            estimator=np.mean,
            errorbar=None,
            palette=["gray"] * df_plot["sujet"].nunique(),
            markers="o",
            linestyles="-",
            linewidth=1,
            alpha=0.7,
            legend=False,
            ax=ax
        )

    g.set_axis_labels("", "Value")
    g.set_titles("{col_name}")

    plt.show()



################################
######## EXECUTE ########
################################


if __name__ == '__main__':

    export_fig_psychometric()

    
    






















    


    


















