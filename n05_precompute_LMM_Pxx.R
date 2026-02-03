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



############################
#### Pxx WHOLE CYCLE ####
############################

root = "/home/jules/Documents/RRET_JULES/Analyses/precompute/TF/session/df_R"
outputdir_fig = "/home/jules/Documents/RRET_JULES/Analyses/results/LMM/diagnosis"
outputdir_df_lmm = "/home/jules/Documents/RRET_JULES/Analyses/results/LMM/df"

band_list <- c("theta", "alpha", "beta", "gamma")

for (band in band_list) {

  # Load the Excel data
  filename = paste0("/df_R_", band, ".xlsx")
  df_raw <- read_excel(paste(root, filename, sep  = "/"))
  df_raw$chan <- paste0(df_raw$sujet, "_", df_raw$chan)
  
  ROI_list_raw <- unique(df_raw$ROI)
  ROI_list <- ROI_list_raw[!grepl("UNSORTED", ROI_list_raw)]
  
  ROI_sel = ROI_list[1]
  
  
  for (ROI_sel in ROI_list) {
    
    withCallingHandlers({
      
      print(ROI_sel)
      
      df_oneROI <- subset(df_raw, ROI == ROI_sel)
      
      print(subset(df_oneROI) %>% count(sujet))
      
      df <- df_oneROI[c("sujet", "chan", "cycle", "Pxx", "resp", "state")]
      
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
      p <- ggplot(df, aes(x = sujet, y = Pxx, color = sujet, fill = sujet)) +
        geom_boxplot(width = .2, alpha = .5, outlier.alpha = 0,
                     position = position_dodge(.9)) +
        stat_summary(fun = median, geom = "point", size = 2,
                     position = position_dodge(.9), color = "white")+
        labs(
          title    = paste(ROI_sel, band, "Pxx", sep = "_")
        ) +
        theme(
          plot.title    = element_text(hjust = 0.5),
        )
      
      p
      
      file_boxplot_subjectwise = paste(band, "boxplot", ROI_sel, "Pxx_subjectwise.png", sep = "_")
      # then explicitly:
      ggsave(paste(outputdir_fig, file_boxplot_subjectwise, sep = "/"), plot = p, width = 8, height = 5)
      
      #### MODEL
      #complex_form <- Pxx ~ resp * state + (resp | sujet/chan)
      #simple_form  <- Pxx ~ resp * state + (1 | sujet/chan)
      simple_form  <- Pxx ~ resp * state + (1 | sujet)
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
      filename_hist = paste(band, "histogram", ROI_sel, "Pxx.png", sep = "_")
      
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
        main     = paste(band, ROI_sel, "Pxx", "kurtosis:", kurt_chan, "skewness", skew_chan),
        adj      = 0.5,       # 0.5 = center
        cex.main = 1.5,       # main title size
        font.main= 2,         # bold
        cex.sub  = 1.0        # subtitle size
      )
      
      dev.off()
      
      #### FIG 3
      filename_qqplot = paste(band, "qqplot", ROI_sel, "Pxx.png", sep = "_")
      
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
      
      filesxlsx_ROI = paste(band, "Pxx_lmm", ROI_sel, "res.xlsx", sep = "_")
      writexl::write_xlsx(model_df, paste(outputdir_df_lmm, filesxlsx_ROI, sep = "/"))
      
    }, warning = function(w) {
      message("Warning: ", conditionMessage(w))  # shows immediately
      invokeRestart("muffleWarning")
    })
    
    
  }
  
}





####################
#### Pxx PHASE ####
####################

root = "/home/jules/Documents/RRET_JULES/Analyses/precompute/TF/session/df_R"
outputdir_fig = "/home/jules/Documents/RRET_JULES/Analyses/results/LMM/diagnosis"
outputdir_df_lmm = "/home/jules/Documents/RRET_JULES/Analyses/results/LMM/df"

band_list <- c("theta", "alpha", "beta", "gamma")
phase_list <- c("inspi", "expi")



for (band in band_list) {
  
  # Load the Excel data
  filename = paste0("/df_R_", band, "_phase.xlsx")
  df_raw <- read_excel(paste(root, filename, sep  = "/"))
  df_raw$chan <- paste0(df_raw$sujet, "_", df_raw$chan)
  
  ROI_list_raw <- unique(df_raw$ROI)
  ROI_list <- ROI_list_raw[!grepl("UNSORTED", ROI_list_raw)]
  
  phase = phase_list[1]
  
  for (phase in phase_list) {
  
    ROI_sel = ROI_list[1]
    
    
    df_phase <- subset(df_raw, phase == phase)
    
    for (ROI_sel in ROI_list) {
      
      withCallingHandlers({
        
        print(ROI_sel)
        
        df_oneROI <- subset(df_phase, ROI == ROI_sel)
        
        print(subset(df_oneROI) %>% count(sujet))
        
        df <- df_oneROI[c("sujet", "chan", "cycle", "Pxx", "resp", "state")]
        
        df$resp  <- factor(df$resp)   # ensure factor
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
            title    = paste(ROI_sel, band, phase, "Pxx", sep = "_")
          ) +
          theme(
            plot.title    = element_text(hjust = 0.5),
          )
        
        p
        
        file_boxplot_subjectwise = paste(band, "boxplot", ROI_sel, phase, "Pxx_subjectwise.png", sep = "_")
        # then explicitly:
        ggsave(paste(outputdir_fig, file_boxplot_subjectwise, sep = "/"), plot = p, width = 8, height = 5)
        
        #### MODEL
        #complex_form <- Pxx ~ resp * state + (resp | sujet/chan)
        #simple_form  <- Pxx ~ resp * state + (1 | sujet/chan)
        simple_form  <- Pxx ~ resp * state + (1 | sujet)
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
        filename_hist = paste(band, "histogram", ROI_sel, phase, "Pxx.png", sep = "_")
        
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
          main     = paste(band, ROI_sel, phase, "Pxx", "kurtosis:", kurt_chan, "skewness", skew_chan),
          adj      = 0.5,       # 0.5 = center
          cex.main = 1.5,       # main title size
          font.main= 2,         # bold
          cex.sub  = 1.0        # subtitle size
        )
        
        dev.off()
        
        #### FIG 3
        filename_qqplot = paste(band, "qqplot", ROI_sel, phase, "Pxx.png", sep = "_")
        
        png(
          filename = paste(outputdir_fig, filename_qqplot, sep = "/"),
          width    = 800,    # width in pixels
          height   = 600,    # height in pixels
          res      = 100     # resolution (pixels per inch)
        )
        
        qqnorm(resid(model))
        qqline(resid(model))  # points fall nicely onto the line - good!
        
        title(
          sub     = paste(band, ROI_sel, phase, "qqplot"),
          adj      = 0.5,       # 0.5 = center
          cex.main = 1.5,       # main title size
          font.main= 2,         # bold
          cex.sub  = 1.0        # subtitle size
        )
        
        dev.off()
        
        #### EXPORT MODEL RES
        tab_model(model, show.re.var = TRUE, show.icc = TRUE, show.r2 = TRUE, show.se = TRUE)
        
        model_df <- broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE)
        
        filesxlsx_ROI = paste(band, "Pxx_lmm", ROI_sel, phase, "res.xlsx", sep = "_")
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

# Load the Excel data
filename = paste0("/df_reg_allsujet_whole_cycle.xlsx")
df_raw <- read_excel(paste(root, filename, sep  = "/"))
df_raw$chan <- paste0(df_raw$sujet, "_", df_raw$chan)

ROI_list_raw <- unique(df_raw$loca)
ROI_list <- ROI_list_raw[!grepl("UNSORTED", ROI_list_raw)]

rf_metrics <- unique(df_raw$rf_metric)

for (rf_metric in rf_metrics) {
  
  df_rf_metric <- subset(df_raw, rf_metric == rf_metric)

  for (band in band_list) {
    
    df_band <- subset(df_rf_metric, band == band)
    
    ROI_sel = ROI_list[1]
    
    for (ROI_sel in ROI_list) {
      
      withCallingHandlers({
        
        print(ROI_sel)
        
        df_oneROI <- subset(df_band, loca == ROI_sel)
  
        print(subset(df_oneROI) %>% count(sujet))
        
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
            title    = paste(ROI_sel, band, "OLS_a", sep = "_")
          ) +
          theme(
            plot.title    = element_text(hjust = 0.5),
          )
        
        p
        
        file_boxplot_subjectwise = paste(band, rf_metric, "boxplot", ROI_sel, "OLS_a_subjectwise.png", sep = "_")
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
        filename_hist = paste(band, rf_metric, "histogram", ROI_sel, "OLS_a.png", sep = "_")
        
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
          main     = paste(band, rf_metric, ROI_sel, "OLS_a", "kurtosis:", kurt_chan, "skewness", skew_chan),
          adj      = 0.5,       # 0.5 = center
          cex.main = 1.5,       # main title size
          font.main= 2,         # bold
          cex.sub  = 1.0        # subtitle size
        )
        
        dev.off()
        
        #### FIG 3
        filename_qqplot = paste(band, rf_metric, "qqplot", ROI_sel, "OLS_a.png", sep = "_")
        
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
        
        filesxlsx_ROI = paste(band, rf_metric, "OLS_a_lmm", ROI_sel, "res.xlsx", sep = "_")
        writexl::write_xlsx(model_df, paste(outputdir_df_lmm, filesxlsx_ROI, sep = "/"))
        
      }, warning = function(w) {
        message("Warning: ", conditionMessage(w))  # shows immediately
        invokeRestart("muffleWarning")
      })
      
      
    }
    
  }

}