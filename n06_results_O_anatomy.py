



from n00_config_O_params import *
from n00bis_config_O_analysis_functions import *
from n00ter_X_manip_data import *

import plotly.express as px





########################################
######## EXPORT ANAT PLOT ########
########################################

def export_anat_plot():

    df_loca_allsujet = get_df_loca_allsujet()
    
    ROI_allsujet = [_loca for _loca in df_loca_allsujet['loca'].unique() if _loca.find('UNSORTED') == -1]
    df_sujet_count = df_loca_allsujet.query(f"loca in {ROI_allsujet}")
    df_sujet_count = df_sujet_count.groupby("loca")["sujet"].nunique().reset_index(name="count")

    for sujet_n_thresh in np.arange(1,5):
        
        df_plot = df_sujet_count.query(f"count >= {sujet_n_thresh}")
        fig = px.bar(df_plot, x="loca", y="count", title=f"Sujet count thresh:{sujet_n_thresh}")
        # fig.show()

        os.chdir(os.path.join(path_results, 'anatomy'))
        fig.write_html(f"sujet_count_{sujet_n_thresh}thresh.html")





################################
######## EXECUTE ########
################################


if __name__ == '__main__':

    
    export_anat_plot()



                        