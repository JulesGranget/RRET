

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
######## ENVELOPPE ########
################################




def export_enveloppe():

    #ROI_sel = ROI_short_list[0]
    for ROI_sel in ROI_short_list:

        print(ROI_sel)

        #### get chan to load
        localist_allsujet = get_df_loca_allsujet()
        localist_allsujet_ROI_sel = localist_allsujet.query(f"loca == '{ROI_sel}'")

        sujet_list_ROI_sel = []
        chan_sel_list_ROI_sel = {}

        for sujet_ROI_sel in localist_allsujet_ROI_sel['sujet'].unique():

            localist_sujet_ROI_sel = localist_allsujet_ROI_sel.query(f"sujet == '{sujet_ROI_sel}'").copy()
            sujet_list_ROI_sel.append(sujet_ROI_sel)
            chan_sel_list_ROI_sel[sujet_ROI_sel] = []

            for chan in localist_sujet_ROI_sel['chan'].values:

                chan_sel_list_ROI_sel[sujet_ROI_sel].append(chan)

        #### load Pxx
        os.chdir(os.path.join(path_precompute, 'TF', 'session'))

        ROI_sel_data_allcond = []

        #sujet = sujet_list_ROI_sel[0]
        for sujet_i, sujet in enumerate(sujet_list_ROI_sel):

            print_advancement(sujet_i, len(sujet_list_ROI_sel), [25,50,75])

            #### params
            sujet_chanlist, sujet_localist = get_chanlist(sujet)

            #### load
            tf_stretch_cond = []
            
            for cond in conditions:
                
                _tf = np.load(f'{sujet}_{cond}_tf_allchan_stretch_post.npy')  
                frex_sel = (frex > freq_band_dict['gamma'][0]) & (frex < freq_band_dict['gamma'][1])
                _tf_gamma = np.median(np.median(_tf[:,:,frex_sel], axis=2), axis=1)
                tf_stretch_cond.append(_tf_gamma)

            tf_stretch_cond = np.stack(tf_stretch_cond)

            #### diff
            tf_diff = []
            tf_baseline = tf_stretch_cond[0]

            for cond_i, cond in enumerate(conditions):

                if cond == 'rsp_ctrl':
                    continue

                tf_cond = tf_stretch_cond[cond_i]
                _tf_diff = tf_cond - tf_baseline
                tf_diff.append(_tf_diff)

            tf_diff = np.stack(tf_diff)

            #### inspcet
            if debug:

                chan_i = 0
                
                for cond_i, cond in enumerate(conditions[1:]):
                
                    plt.plot(tf_diff[cond_i][chan_i], label=cond)
                
                plt.legend()
                plt.show()

            #### extract
            sujet_chan_i_sel_list = np.array([np.where(sujet_chanlist == _chan)[0][0] for _chan in chan_sel_list_ROI_sel[sujet]])

            tf_diff_chan_sel = tf_diff[:, sujet_chan_i_sel_list]

            ROI_sel_data_allcond.append(tf_diff_chan_sel)

        ROI_sel_data_allcond = np.concat(ROI_sel_data_allcond, axis=1)

        fig, ax = plt.subplots()

        for cond_i, cond in enumerate(conditions[1:]):

            sig_med = np.median(ROI_sel_data_allcond[cond_i], axis=0) 
            ax.plot(sig_med, label=cond)

        plt.suptitle(f"{ROI_sel} s({len(sujet_list_ROI_sel)}), c({ROI_sel_data_allcond[cond_i].shape[0]})")
        plt.legend()

        # plt.show()

        os.chdir(os.path.join(path_results, 'Pxx', 'ERP_gamma_enveloppe'))
        fig.savefig(f"{ROI_sel}.png")

        plt.close('all')





################################
######## EXECUTE ########
################################


if __name__ == '__main__':

    export_enveloppe()


                        