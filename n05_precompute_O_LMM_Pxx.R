# Required libraries install.packages("lme4")

library(readxl)
library(Matrix)
library(lme4)
library(lmerTest)  # For p-values
library(emmeans)
library(ggplot2)
library(e1071)  # For skewness and kurtosis
library(pbkrtest)
library(sjPlot)
library(broom.mixed)
library(writexl)
library(dplyr)



########################
######## Pxx ########
########################

root = "/home/jules/Documents/RRET_JULES/Analyses/precompute/TF/session/df_R"
outputdir_fig = "/home/jules/Documents/RRET_JULES/Analyses/results/LMM/diagnosis"
outputdir_df_lmm = "/home/jules/Documents/RRET_JULES/Analyses/results/LMM/df"

band_list <- c("theta", "alpha", "beta", "gamma")
phase_list <- c("whole", "inspi", "expi")

phase_cycle_sel = phase_list[1]
band = band_list[1]

for (band in band_list) {
  
  # Load the Excel data
  filename = paste0("/df_R_", band, ".xlsx")
  df_raw <- read_excel(paste(root, filename, sep  = "/"))
  df_raw$chan <- paste0(df_raw$sujet, "_", df_raw$chan)
  
  ROI_list_raw <- unique(df_raw$ROI)
  ROI_list <- ROI_list_raw[!grepl("UNSORTED", ROI_list_raw)]
  
  for (phase_cycle_sel in phase_list) {
    
    df_phase <- subset(df_raw, phase_cycle == phase_cycle_sel)
    
    ROI_sel = ROI_list[1]
    
    for (ROI_sel in ROI_list) {
      
      withCallingHandlers({
        
        print(ROI_sel)
        
        filename_export_diagnostic = paste(band, phase_cycle_sel, ROI_sel, sep = "_")
        
        df_oneROI <- subset(df_phase, ROI == ROI_sel)
        
        df_count <- df_oneROI %>%
          group_by(sujet) %>%
          summarise(
            n_chan  = n_distinct(chan),
            n_cycle = n_distinct(cycle),
            .groups = "drop"
          )
        
        print(df_count)
        
        df <- df_oneROI[c("sujet", "chan", "cycle", "Pxx", "resp", "state")]
        
        # Convert categorical variables to factors
        df$sujet <- factor(df$sujet)
        df$chan <- factor(df$chan)
        df$resp  <- factor(df$resp)
        df$state <- factor(df$state)
        
        df$resp  <- relevel(df$resp,  ref = "rsp")   # baseline for resp
        df$state <- relevel(df$state, ref = "ctrl")  # baseline for state
        
        #### FIG 1
        p <- ggplot(df, aes(x = sujet, y = Pxx, color = sujet, fill = sujet)) +
          geom_boxplot(width = .2, alpha = .5, outlier.alpha = 0,
                       position = position_dodge(.9)) +
          stat_summary(fun = median, geom = "point", size = 2,
                       position = position_dodge(.9), color = "white")+
          labs(
            title    = paste(filename_export_diagnostic, "Pxx", sep = "_")
          ) +
          theme(
            plot.title    = element_text(hjust = 0.5),
          )
        
        p
        
        file_boxplot_subjectwise = paste("BOXPLOT", filename_export_diagnostic, "Pxx_subjectwise.png", sep = "_")
        # then explicitly:
        ggsave(paste(outputdir_fig, file_boxplot_subjectwise, sep = "/"), plot = p, width = 8, height = 5)
        
        #### MODEL
        #simple_form <- Pxx ~ resp * state + (resp | sujet/chan)
        #simple_form  <- Pxx ~ resp * state + (1 | sujet/chan)
        simple_form  <- Pxx ~ resp * state + (1 | sujet/chan)
        simple_form_refit  <- Pxx ~ resp * state
        
        model <- tryCatch({
          
          warn_triggered <- FALSE  # will catch if any warning is raised
          
          mod_attempt <- withCallingHandlers(
            expr = {
              lmer(
                simple_form,
                data = df,
                control = lmerControl(optCtrl = list(maxfun = 2e5))
              )
            },
            warning = function(w) {
              message("⚠️ Warning during lmer(): ", conditionMessage(w))
              warn_triggered <<- TRUE
              invokeRestart("muffleWarning")  # suppress so execution continues
            }
          )
          
          # Force fallback if warning was raised or model is singular
          if (warn_triggered || isSingular(mod_attempt, tol = 1e-4)) {
            message("⚠️️ Fallback to lm() due to warning or singular fit")
            stop("Trigger fallback to lm")
          }
          
          mod_attempt  # return valid model if all checks passed
          
        }, error = function(e) {
          lm(simple_form_refit, data = df)
        })
        
        summary(model)
        
        
        #### FIG 2
        filename_hist = paste(band, "histogram", ROI_sel, phase_cycle_sel, "Pxx.png", sep = "_")
        
        skew_chan = round(skewness(df$Pxx), 2)
        kurt_chan = round(kurtosis(df$Pxx), 2)
        
        png(
          filename = paste(outputdir_fig, filename_hist, sep = "/"),
          width    = 800,    # width in pixels
          height   = 600,    # height in pixels
          res      = 100     # resolution (pixels per inch)
        )
        
        hist(
          df$Pxx,
          breaks = 30,
          main   = "",          # leave main blank for now
          xlab   = "Pxx values",
          ylab   = "Resp",
          col    = "lightblue",
          border = "white"
        )
        
        title(
          main     = paste(band, ROI_sel, phase_cycle_sel, "Pxx", "kurtosis:", kurt_chan, "skewness", skew_chan),
          adj      = 0.5,       # 0.5 = center
          cex.main = 1.5,       # main title size
          font.main= 2,         # bold
          cex.sub  = 1.0        # subtitle size
        )
        
        dev.off()
        
        #### FIG 3
        filename_qqplot = paste(band, "qqplot", ROI_sel, phase_cycle_sel, "Pxx.png", sep = "_")
        
        png(
          filename = paste(outputdir_fig, filename_qqplot, sep = "/"),
          width    = 800,    # width in pixels
          height   = 600,    # height in pixels
          res      = 100     # resolution (pixels per inch)
        )
        
        qqnorm(resid(model))
        qqline(resid(model))  # points fall nicely onto the line - good!
        
        title(
          sub     = paste(band, ROI_sel, "qqplot"),
          adj      = 0.5,       # 0.5 = center
          cex.main = 1.5,       # main title size
          font.main= 2,         # bold
          cex.sub  = 1.0        # subtitle size
        )
        
        dev.off()
        
        #### EXPORT MODEL RES
        tab_model(model, show.re.var = TRUE, show.icc = TRUE, show.r2 = TRUE, show.se = TRUE)
        
        model_df <- broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE)
        
        filesxlsx_ROI = paste(band, "Pxx_lmm", ROI_sel, phase_cycle_sel, "res.xlsx", sep = "_")
        writexl::write_xlsx(model_df, paste(outputdir_df_lmm, filesxlsx_ROI, sep = "/"))
        
      }, warning = function(w) {
        message("Warning: ", conditionMessage(w))  # shows immediately
        invokeRestart("muffleWarning")
      })
      
      
    }
    
  }
  
}









################
#### OLS_a ####
################

root = "/home/jules/Documents/RRET_JULES/Analyses/precompute/TF/session/df_R"
outputdir_fig = "/home/jules/Documents/RRET_JULES/Analyses/results/LMM/diagnosis"
outputdir_df_lmm = "/home/jules/Documents/RRET_JULES/Analyses/results/LMM/df"

band_list <- c("theta", "alpha", "beta", "gamma")
phase_cycle_list <- c("whole", "inspi", "expi")
pre_post_list <- c("pre", "post") 

# Load the Excel data
filename = paste0("/df_reg_R.xlsx")
df_raw <- read_excel(paste(root, filename, sep  = "/"))
df_raw$chan <- paste0(df_raw$sujet, "_", df_raw$chan)

ROI_list_raw <- unique(df_raw$ROI)
ROI_list <- ROI_list_raw[!grepl("UNSORTED", ROI_list_raw)]

rf_metrics <- unique(df_raw$rf_metric)

phase_cycle_sel = phase_cycle_list[1]
rf_metric_sel = rf_metrics[1]
band_sel = band_list[1]
ROI_sel = ROI_list[1]
pre_post_sel = pre_post_list[1]

for (pre_post_sel in pre_post_list) {
  
  df_pre_post <- subset(df_raw, phase_protocol == pre_post_sel)

  for (phase_cycle_sel in phase_cycle_list) {
    
    df_phase_cycle <- subset(df_pre_post, phase_cycle == phase_cycle_sel)
  
    for (rf_metric_sel in rf_metrics) {
      
      df_rf_metric <- subset(df_phase_cycle, rf_metric == rf_metric_sel)
    
      for (band_sel in band_list) {
        
        df_band <- subset(df_rf_metric, band == band_sel)
        
        for (ROI_sel in ROI_list) {
          
          withCallingHandlers({
            
            print(ROI_sel)
            
            df_oneROI <- subset(df_band, ROI == ROI_sel)
      
            print(df_oneROI %>% group_by(sujet) %>% summarise(n_chan = n_distinct(chan)))
            
            df <- df_oneROI[c("sujet", "chan", "OLS_a", "resp", "state")]
            
            # Convert categorical variables to factors
            df$sujet <- as.factor(df$sujet)
            df$chan <- as.factor(df$chan)
            #df$cond <- as.factor(df$cond)
            
            # ---- Set baselines (reference levels) ----
            df$resp  <- factor(df$resp)   # ensure factor
            df$state <- factor(df$state)
            
            df$resp  <- relevel(df$resp,  ref = "rsp")   # baseline for resp
            df$state <- relevel(df$state, ref = "ctrl")  # baseline for state
            
            #### FIG 1
            p <- ggplot(df, aes(x = sujet, y = OLS_a, color = sujet, fill = sujet)) +
              geom_boxplot(width = .2, alpha = .5, outlier.alpha = 0,
                           position = position_dodge(.9)) +
              stat_summary(fun = median, geom = "point", size = 2,
                           position = position_dodge(.9), color = "white")+
              labs(
                title    = paste(ROI_sel, band_sel, phase_cycle_sel, rf_metric_sel, pre_post_sel, "OLS_a", sep = "_")
              ) +
              theme(
                plot.title    = element_text(hjust = 0.5),
              )
            
            p
            
            file_boxplot_subjectwise = paste(band_sel, rf_metric_sel, "boxplot", ROI_sel, phase_cycle_sel, pre_post_sel, "OLS_a_subjectwise.png", sep = "_")
            # then explicitly:
            ggsave(paste(outputdir_fig, file_boxplot_subjectwise, sep = "/"), plot = p, width = 8, height = 5)
            
            #### MODEL
            #complex_form <- OLS_a ~ resp * state + (resp | sujet/chan)
            #simple_form  <- OLS_a ~ resp * state + (1 | sujet/chan)
            simple_form  <- OLS_a ~ resp * state + (1 | sujet)
            simple_form_refit  <- OLS_a ~ resp * state
            
            model <- tryCatch({
              
              warn_triggered <- FALSE  # will catch if any warning is raised
              
              mod_attempt <- withCallingHandlers(
                expr = {
                  lmer(
                    simple_form,
                    data = df,
                    control = lmerControl(optCtrl = list(maxfun = 2e5))
                  )
                },
                warning = function(w) {
                  message("⚠️ Warning during lmer(): ", conditionMessage(w))
                  warn_triggered <<- TRUE
                  invokeRestart("muffleWarning")  # suppress so execution continues
                }
              )
              
              # Force fallback if warning was raised or model is singular
              if (warn_triggered || isSingular(mod_attempt, tol = 1e-4)) {
                message("⚠️️ Fallback to lm() due to warning or singular fit")
                stop("Trigger fallback to lm")
              }
              
              mod_attempt  # return valid model if all checks passed
              
            }, error = function(e) {
              lm(simple_form_refit, data = df)
            })
            
            summary(model)
            
            
            #### FIG 2
            filename_hist = paste(band_sel, rf_metric_sel, "histogram", ROI_sel, phase_cycle_sel, pre_post_sel, "OLS_a.png", sep = "_")
            
            skew_chan = round(skewness(df$OLS_a), 2)
            kurt_chan = round(kurtosis(df$OLS_a), 2)
            
            png(
              filename = paste(outputdir_fig, filename_hist, sep = "/"),
              width    = 800,    # width in pixels
              height   = 600,    # height in pixels
              res      = 100     # resolution (pixels per inch)
            )
            
            hist(
              df$OLS_a,
              breaks = 30,
              main   = "",          # leave main blank for now
              xlab   = "OLS_a values",
              ylab   = "Resp",
              col    = "lightblue",
              border = "white"
            )
            
            title(
              main     = paste(band_sel, rf_metric_sel, ROI_sel, phase_cycle_sel, pre_post_sel, "OLS_a", "kurtosis:", kurt_chan, "skewness", skew_chan),
              adj      = 0.5,       # 0.5 = center
              cex.main = 1.5,       # main title size
              font.main= 2,         # bold
              cex.sub  = 1.0        # subtitle size
            )
            
            dev.off()
            
            #### FIG 3
            filename_qqplot = paste(band_sel, rf_metric_sel, "qqplot", ROI_sel, phase_cycle_sel, pre_post_sel, "OLS_a.png", sep = "_")
            
            png(
              filename = paste(outputdir_fig, filename_qqplot, sep = "/"),
              width    = 800,    # width in pixels
              height   = 600,    # height in pixels
              res      = 100     # resolution (pixels per inch)
            )
            
            qqnorm(resid(model))
            qqline(resid(model))  # points fall nicely onto the line - good!
            
            title(
              sub     = paste(band, rf_metric, ROI_sel, "qqplot"),
              adj      = 0.5,       # 0.5 = center
              cex.main = 1.5,       # main title size
              font.main= 2,         # bold
              cex.sub  = 1.0        # subtitle size
            )
            
            dev.off()
            
            #### EXPORT MODEL RES
            tab_model(model, show.re.var = TRUE, show.icc = TRUE, show.r2 = TRUE, show.se = TRUE)
            
            model_df <- broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE)
            
            filesxlsx_ROI = paste(band_sel, rf_metric_sel, "OLS_a_lmm", ROI_sel, phase_cycle_sel, pre_post_sel, "res.xlsx", sep = "_")
            writexl::write_xlsx(model_df, paste(outputdir_df_lmm, filesxlsx_ROI, sep = "/"))
            
          }, warning = function(w) {
            message("Warning: ", conditionMessage(w))  # shows immediately
            invokeRestart("muffleWarning")
          })
          
          
        }
        
      }
    
    }
    
  }

}






