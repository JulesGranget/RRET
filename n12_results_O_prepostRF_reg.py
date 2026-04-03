



from n00_config_O_params import *
from n00bis_config_O_analysis_functions import *
from n00ter_X_manip_data import *


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

    #### plot

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






