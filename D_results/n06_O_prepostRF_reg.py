



from A_config.n01_O_params import *
from A_config.n02_O_analysis_functions import *
from A_config.n03_X_manip_data import *


import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go





################################
######## FUNCTIONS ########
################################





def p_to_stars(p):
    if p < 0.001:
        return '***'
    elif p < 0.01:
        return '**'
    elif p < 0.05:
        return '*'
    else:
        return ''







########################################
######## PREPOST EXAMPLES ########
########################################





def export_reg_PREPOST_example():

    #### load
    filename = os.path.join(path_precompute, 'RESP', 'df_R', f"df_R_RFonly_selOC.xlsx")
    df_reg_prepost = pd.read_excel(filename)
    cond_list_sel = ['oc_ctrl', 'oc_chl']

    #### plot all patient
    for cond in cond_list_sel:
    
        df_plot = df_reg_prepost.query(f"cond == '{cond}'")[['sujet', 'pre_total_amplitude', 'post_oc_ratio']]

        sns.lmplot(data=df_plot, x="pre_total_amplitude", y="post_oc_ratio", hue="sujet", ci=None,
            scatter_kws={"alpha": 0.2},  # make points transparent
            line_kws={"linewidth": 3}  # optional: thicker regression line
        )

        plt.subplots_adjust(top=0.92)
        plt.suptitle(f"{cond}")

        # plt.show()

        filename = os.path.join(path_results, "respi", "pre_amp", f"lmplot_amp_OCratio_{cond}.png")
        plt.savefig(filename)

        plt.close('all')

    #### plot example pre post
    cond_prepost_list = ['oc_ctrl', 'oc_chl']

    i_to_plot_paper = {'oc_ctrl' : 23, 'oc_chl' : 27}
    index_search = {'oc_ctrl' : 20, 'oc_chl' : 20}

    sujet_sel_ex_reg = ['NS203', 'NS211']

    #cond = 'oc_chl'
    for cond in cond_prepost_list:

        df_prepost_example = df_reg_prepost.query(f"cond == '{cond}' and sujet in {sujet_sel_ex_reg}").reset_index(drop=True)

        cycle_i_list = []

        for row_i, row in df_prepost_example.iterrows():

            if row_i == 0:
                pre_sujet = row['sujet']
                add_cycle_i = 0
            else:
                if pre_sujet != row['sujet']:
                    add_cycle_i = 0
                    pre_sujet = row['sujet']
                else:
                    add_cycle_i += 1

            cycle_i_list.append(add_cycle_i)

        df_prepost_example['cycle_i'] = cycle_i_list

        prepost_norm = df_prepost_example['pre_total_amplitude'] / df_prepost_example['pre_total_amplitude'].max()
        ocratio_norm = 1 - df_prepost_example['post_oc_ratio'] / df_prepost_example['post_oc_ratio'].max()

        score_prepost = (prepost_norm + ocratio_norm).sort_values(ascending=False)
        score_prepost_i_list = np.concat([score_prepost.index[:index_search[cond]].values, score_prepost.index[-index_search[cond]:].values])

        load_df = df_prepost_example.loc[score_prepost_i_list][['sujet', 'cond']]
        score_prepost_i_list = load_df.sort_values('sujet', ascending=True).index.values

        load_dict = {}

        for row_i, row in load_df.iterrows():

            _sujet, _cond = row["sujet"], row["cond"]

            _ = load_dict.setdefault(_sujet, {})
            _ = load_dict[_sujet].setdefault(_cond, {})

        path_load_resp = os.path.join(path_precompute, 'RESP', 'session')

        vlim = []

        for sujet in load_dict:

            for cond in load_dict[sujet]:

                _pre_stretch = np.load(os.path.join(path_load_resp, f"{sujet}_{cond}_stretch_resp_pre.npy"))
                _post_stretch = np.load(os.path.join(path_load_resp, f"{sujet}_{cond}_stretch_resp_post.npy"))
                load_dict[sujet][cond]['pre'] = _pre_stretch
                load_dict[sujet][cond]['post'] = _post_stretch

                vlim.append([_pre_stretch.min(), _pre_stretch.max(), _post_stretch.min(), _post_stretch.max()])

        vlim = np.abs(np.array(vlim)).max()

        for cycle_i, cycle_i_val in enumerate(score_prepost_i_list):

            if cycle_i in [0,1]:
                continue

            _row_0 = df_prepost_example.loc[score_prepost_i_list[cycle_i-1]]
            _row_1 = df_prepost_example.loc[cycle_i_val]
            
            _sujet_0, _cond_0, cycle_i_mat0 = _row_0['sujet'], _row_0['cond'], _row_0['cycle_i']
            _sujet_1, _cond_1, cycle_i_mat1 = _row_1['sujet'], _row_1['cond'], _row_1['cycle_i']

            stretch_cycle_0 = np.concat([load_dict[_sujet_0][_cond_0]['pre'], load_dict[_sujet_0][_cond_0]['post']], axis=1)[cycle_i_mat0]
            stretch_cycle_1 = np.concat([load_dict[_sujet_1][_cond_1]['pre'], load_dict[_sujet_1][_cond_1]['post']], axis=1)[cycle_i_mat1]

            plt.plot(stretch_cycle_0)
            plt.plot(stretch_cycle_1)
            plt.title(f"{_sujet_0}, {_sujet_1}, {_cond}, iteration{cycle_i}")
            plt.ylim(-vlim, vlim)
            plt.show()

        #### save paper
        cycle_i, cycle_i_val = i_to_plot_paper[cond], score_prepost_i_list[i_to_plot_paper[cond]]

        _row_0 = df_prepost_example.loc[score_prepost_i_list[cycle_i-1]]
        _row_1 = df_prepost_example.loc[cycle_i_val]
        
        _sujet_0, _cond_0, cycle_i_mat0 = _row_0['sujet'], _row_0['cond'], _row_0['cycle_i']
        _sujet_1, _cond_1, cycle_i_mat1 = _row_1['sujet'], _row_1['cond'], _row_1['cycle_i']

        stretch_cycle_0 = np.concat([load_dict[_sujet_0][_cond_0]['pre'], load_dict[_sujet_0][_cond_0]['post']], axis=1)[cycle_i_mat0]
        stretch_cycle_1 = np.concat([load_dict[_sujet_1][_cond_1]['pre'], load_dict[_sujet_1][_cond_1]['post']], axis=1)[cycle_i_mat1]

        vlim = np.abs([stretch_cycle_0, stretch_cycle_1]).max()*1.2

        fig, ax = plt.subplots()

        ax.plot(stretch_cycle_0)
        ax.plot(stretch_cycle_1)

        ax.vlines([int(stretch_point_TF/2), stretch_point_TF+int(stretch_point_TF/2)], ymin=-vlim, ymax=vlim, color='k')

        ax.set_title(f"{_sujet_0}, {_sujet_1}, {_cond}, iteration{cycle_i}")
        ax.set_ylim(-vlim, vlim)

        ax.set_ylabel('r-zscore')
        ax.set_xlabel('Phase')

        ax.xaxis.label.set_size(15)
        ax.yaxis.label.set_size(15)
    
        ax.tick_params(axis="x", labelsize=15)
        ax.tick_params(axis="y", labelsize=15)

        ax.set_xticks(np.arange(0,stretch_point_TF*2,50), 
                      labels=np.concat([np.arange(0,stretch_point_TF,50), np.arange(0,stretch_point_TF,50)]))

        plt.tight_layout()

        # plt.show()

        fig.savefig(os.path.join(path_paper_figure_export, f"fig06a_example_{cond}.svg"))

    #### plot lm patient wise
    for sujet in sujet_list:

        df_plot = df_reg_prepost.query(f"sujet == '{sujet}' and cond in ['oc_ctrl', 'oc_chl']")
        df_plot['cond'] = df_plot['cond'].replace({'oc_ctrl' : 'O', 'oc_chl' : 'O+Ch'})
        g = sns.lmplot(df_plot, x='pre_total_amplitude', y='post_oc_ratio', hue='cond', ci=False)
        plt.title(f"{sujet}")
        plt.tight_layout()

        if sujet == 'NS211':

            ax = g.ax
            ax.set_xlabel("preA", fontsize=15)
            ax.set_ylabel("OR", fontsize=15)

            ax.tick_params(axis="x", labelsize=15)
            ax.tick_params(axis="y", labelsize=15)

            plt.tight_layout()

            plt.savefig(os.path.join(path_paper_figure_export, f"fig06b_patient_reg_cond.svg"))

        plt.show()

    #### plot allpatient
    df_plot = df_reg_prepost.query(f"cond in ['oc_ctrl', 'oc_chl']")
    df_plot['cond'] = df_plot['cond'].replace({'oc_ctrl' : 'O', 'oc_chl' : 'O+Ch'})
    g = sns.lmplot(df_plot, x='pre_total_amplitude', y='post_oc_ratio', hue='cond', ci=False, scatter_kws={"alpha": 0.3})
    ax = g.ax
    ax.set_xlabel("preA", fontsize=15)
    ax.set_ylabel("OR", fontsize=15)

    ax.tick_params(axis="x", labelsize=15)
    ax.tick_params(axis="y", labelsize=15)

    plt.tight_layout()
    plt.title(f"allpatient")

    plt.savefig(os.path.join(path_paper_figure_export, f"fig06c_allpatient_reg_cond.svg"))

    plt.show()





def export_res_LMM_prepost():

    #### generate df
    cond_list_sel = ['oc_ctrl', 'oc_chl']

    #### load df
    folder_import_dir = os.path.join(path_results, 'LMM', 'RESP', 'df')

    df_LMM_prepost = []

    for cond in cond_list_sel:
    
        _filename = os.path.join(folder_import_dir, f"OC_PRE_RES_{cond}_post_oc_ratio_LMM.xlsx")
        df_LMM_prepost.append(pd.read_excel(_filename))

    df_LMM_prepost = pd.concat(df_LMM_prepost)

    #### plot prepost
    val_signi = 0.05

    df_plot = df_LMM_prepost.query("term == 'pre_total_amplitude' and cond in @cond_list_sel").copy()

    df_plot["cond"] = df_plot["cond"].astype(str)
    df_plot = df_plot.set_index("cond").reindex(cond_list_sel).reset_index()


    x = df_plot["cond"].tolist()
    y = df_plot["estimate"].astype(float).to_numpy()
    p = df_plot["p.value"].astype(float).to_numpy()

    signi = (p < val_signi) & np.isfinite(y)

    # Star offset
    y_max = np.nanmax(np.abs(y)) if np.isfinite(y).any() else 1.0
    y_off = max(0.05 * y_max, 0.05)
    y_star = np.where(y >= 0, y + y_off, y - y_off)

    fig = make_subplots(rows=1, cols=1)

    # Bars
    fig.add_trace(
        go.Bar(
            x=x,
            y=y,
            name='pre_total_amplitude',
            marker=dict(line=dict(width=0)),
            text=[f"p={pv:.2g}" if np.isfinite(pv) else "" for pv in p],
            textposition="outside",
            cliponaxis=False,
            showlegend=False
        ),
        row=1, col=1
    )

    # Stars (only where significant)
    if signi.any():
        fig.add_trace(
            go.Scatter(
                x=[x[i] for i in np.where(signi)[0]],
                y=y_star[signi],
                mode="text",
                text=["★"] * int(signi.sum()),
                textfont=dict(size=18),
                showlegend=False
            ),
            row=1, col=1
        )

    # Zero line
    fig.add_hline(y=0, line_dash="dash", line_color="red", row=1, col=1)

    # A bit more y-range so stars aren't clipped
    y_min = np.nanmin(np.r_[y, y_star[signi]]) if np.isfinite(y).any() else -1
    y_max2 = np.nanmax(np.r_[y, y_star[signi]]) if np.isfinite(y).any() else 1
    margin = 0.15 * (y_max2 - y_min) if (y_max2 > y_min) else 0.2
    fig.update_yaxes(range=[y_min - margin, y_max2 + margin], row=1, col=1)

    fig.update_layout(
        template="simple_white",
        height=450,
        width=600,
        title=f"post_oc_ratio ~ pre_total_amplitude, by condition",
        bargap=0.3
    )

    # fig.show()

    filepath = os.path.join(path_results, 'Pxx', 'reg_with_RF', 'pre_amp', f"REG_RF_preamp.html")
    fig.write_html(filepath)

    #### export df values
    df_export = df_LMM_prepost.query(f"term != '(Intercept)'")
    
    filepath = os.path.join(path_results, 'respi', 'pre_amp')
    filename = os.path.join(filepath, "df_REG_RF_preamp.xlsx")
    df_export.to_excel(filename)







################################
######## EXECUTE ########
################################


if __name__ == '__main__':


    export_reg_PREPOST_example()
    export_res_LMM_prepost()






