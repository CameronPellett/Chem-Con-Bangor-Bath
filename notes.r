library(tidyverse)
library(lme4)
library(ggrepel)
library(ggbeeswarm)
source("C:/Jessica_data_analysis/R/Inventory/AccuracyAssessment.R")
source("C:/Jessica_data_analysis/R/Disertation/rescalechange.r")

#function
#####################
pvalue_fun <- function(p.value){
  if(is.na(p.value)){
    NA_character_
  }else if(p.value > 0.05){
    "> 0.05"
  }else if(p.value < 0.05 & p.value > 0.01){
    "< 0.05"
  }else if(p.value < 0.01 & p.value > 0.001){
    "< 0.01"
  }else if(p.value < 0.001){
    "< 0.001"
  }else {
    "something went wrong"
  }
}

pvalue_star <- function(p.value){
  if(is.na(p.value)){
    NA_character_
  }else if(p.value > 0.05){
    "."
  }else if(p.value < 0.05 & p.value > 0.01){
    "*"
  }else if(p.value < 0.01 & p.value > 0.001){
    "**"
  }else if(p.value < 0.001){
    "***"
  }else {
    "something went wrong"
  }
}

stat_paste_fun <- function(stat, log = FALSE){
  if(log != TRUE){
    paste0(str_remove_all(stat[["method"]][[1]], "\n\t"), ": ", 
           str_remove_all(names(stat[["statistic"]]), "Kruskal-Wallis "), " = ",
           round(stat[["statistic"]][[1]],1), ", ",
           if(length(stat[["parameter"]]) == 0)
           {""
           }else if(is.na(stat[["parameter"]])){
             ""
           }else if(length(stat[["parameter"]]) > 0){
             paste0("df = ", round(stat[["parameter"]], 0), ", ")
           } ,
           "p-value ",  pvalue_fun(stat[["p.value"]]))
  }else{
    paste0(str_remove_all(stat[["method"]][[1]], "\n\t"), " (log y): ", 
           str_remove_all(names(stat[["statistic"]]), "Kruskal-Wallis "), " = ",
           round(stat[["statistic"]][[1]],1), ", ",
           if(length(stat[["parameter"]]) == 0)
           {""
           }else if(is.na(stat[["parameter"]])){
             ""
           }else if(length(stat[["parameter"]]) > 0){
             paste0("df = ", round(stat[["parameter"]], 0), ", ")
           } ,
           "p-value ",  pvalue_fun(stat[["p.value"]]))
  }
}


#####################

#data import
###########################
data <- read_csv("Barbara_full_dataset_Bangor.csv")

flow <- read_csv("flow.csv") 
##########################


#clean
##########################
long <- data|>
  select(-contains("RSV"), -contains("phi6"))|>
  group_by(sample_site_name, timestamp_sample_collected, timestamp_sample_received)|>
  mutate(biorep = row_number())|>
  pivot_longer(cols = contains("_gc_l"), names_to = "con_vir", values_to = "gc_l")|>
  separate(con_vir, into = c("method", "virus"), sep = "_", extra = "drop")|>
  filter(method != "SN")#, !str_detect(virus, "Flu|MeV|RoV|EVD68"))

long_sc <- long|>
  ungroup()|>
  mutate(gc_l_log = log10(gc_l),
         sample_ph_pre_ansis = scale(sample_ph_pre_ansis)[,1], 
         conductivity_ms_cm = scale(conductivity_ms_cm)[,1],
         ammonia_mg_l = scale(ammonia_mg_l)[,1],
         ophosph_mg_l = scale(ophosph_mg_l)[,1],
         Turbidity = scale(Turbidity)[,1])
  
diff <- long_sc|>
  select(-gc_l_log)|>
  pivot_wider(names_from = "method", values_from = "gc_l")|>
  group_by(virus)|>
  mutate(
      #PEG = rescale01(PEG),
      #AM = rescale01(AM),
      PEG_divby_AM = PEG / AM,
        AM_divby_PEG = AM / PEG)|>
  filter(PEG != 0, AM != 0,
         PEG_divby_AM < 10 & AM_divby_PEG < 10)

outliers <- diff|>
  filter(PEG_divby_AM > 10 | AM_divby_PEG > 10)|>
  select(-PEG_divby_AM)|>
  pivot_wider(names_from = "virus", values_from = c("PEG", "AM"))|>
  select("sample_id", "sample_site_name", "sample_site_code", "grab_compo_ind" ,
         "timestamp_sample_collected", "timestamp_sample_received" ,
         contains("CrAss"), contains("N1"), contains("Flu"),
         contains("EVD68"), contains("EV"), contains("NoVGII"),  contains("NoV"))


flow_clean <- flow|> #not done....
  mutate(site = str_replace_all(site, c("Barnoldswick STW" = "BARNOLDSWICK STW",
                                        "Cambridge STW" = "CAMBRIDGE STW"
                                        )))
#################################

diff|>
  ggplot(aes(sample = log(PEG_divby_AM)))+
  geom_qq()+
  geom_qq_line()

diff|>
  ggplot(aes(log(PEG_divby_AM)))+
  geom_histogram(bins = 70)
  
diff|>
  ggplot(aes(sample = log(AM_divby_PEG)))+
  geom_qq()+
  geom_qq_line()

diff|>
  ggplot(aes(log(AM_divby_PEG)))+
  geom_histogram(bins = 70)


#mean comp
###########################
meth_comp_data <- semi_join(long_sc, diff,
                            by = c("sample_site_name",
                                   "timestamp_sample_collected",
                                   "timestamp_sample_received", "biorep", "virus"))

#mean_adj_gc_l <- meth_comp_data|>
#  group_by(virus)|>
#  mutate(mean = mean(gc_l, na.rm = TRUE),
#         gc_l = gc_l - mean,
#         virus = "Combined mean adj.")|>
#  select(-mean)

#meth_comp_data_wadj <- meth_comp_data|>
#  bind_rows(mean_adj_gc_l)

qq_meth_comp <- meth_comp_data|>
  group_by(virus)|>
  mutate(row = row_number())|>
  ggplot(aes(sample = log10(gc_l)))+
  geom_qq()+
  geom_qq_line()+
  facet_grid(method ~ virus, scales = "free")


viruses_t.test <- meth_comp_data|>
  group_by(virus)|>
  summarise(max = max(gc_l, na.rm = T),
            min = min(gc_l, na.rm = T),
            mean = mean(gc_l, na.rm = T))|>
  bind_cols(stat = NA_character_,
            row = 1)

for (i in 1:nrow(viruses_t.test )) {
  filt_virus <- viruses_t.test [i,1][[1]]
  
  viruses_t.test [i,5] <- stat_paste_fun(t.test(gc_l_log ~ method, data = filter(meth_comp_data, virus == filt_virus)), log = TRUE)
  
  
}



meth_comp_data|>
  group_by(virus)|>
  mutate(row = row_number())|>
  left_join(viruses_t.test, by = c("virus", "row"))|>
  ggplot(aes(method, gc_l))+
  geom_boxplot(outlier.shape = NA)+
  geom_beeswarm(cex = 1.5, size = 1.5)+
  geom_beeswarm(cex = 1.5, colour = "white", shape = ".")+
  geom_text(aes(label = stat, x = 1.5, y = max + (0.07 * (max - min))), size = 3)+
  facet_wrap(~virus, scale = "free")

############################

#peg_div_am
######################################
boxcox <- MASS::boxcox(PEG_divby_AM ~ virus + sample_ph_pre_ansis + ammonia_mg_l + ophosph_mg_l + conductivity_ms_cm + Turbidity, 
                       data = diff)

boxcoxlambda <- boxcox$x[which.max(boxcox$y)]

boxcox_trans <- function(x, lambda){
  (x ^ lambda - 1) / lambda
}

lmm_bc <- lmer(PEG_divby_AM ~ (1|virus) + 
                sample_ph_pre_ansis + ammonia_mg_l + ophosph_mg_l + conductivity_ms_cm + Turbidity,
              data = mutate(diff, PEG_divby_AM = boxcox_trans(PEG_divby_AM, boxcoxlambda)))

plot_assumptions(lmm_bc,
                 mutate(diff, PEG_divby_AM = boxcox_trans(PEG_divby_AM, boxcoxlambda))|>
                   filter(!is.na(PEG_divby_AM))|>
                   pull(PEG_divby_AM))

as_tibble(summary(lmm_bc)[["coefficients"]], rownames = "variables")|>
  mutate(variables = str_replace_all(variables, c("methodPEG" = "PEG",
                                                  "methodAM" = "AM",
                                                  "sample_ph_pre_ansis" = "pH",
                                                  "ammonia_mg_l" = "Ammonia",
                                                  "ophosph_mg_l" = "Orthophosphate",
                                                  "conductivity_ms_cm" = "Conductivity")),
         CI90 = (`Std. Error`*1.65),
         CI95 = (`Std. Error`*1.96),
         CI99 = (`Std. Error`*2.59),
         label99 = c("99% CI", rep(NA_character_, 5)),
         label95 = c( "95% CI", rep(NA_character_, 5)),
         label90 = c( "90% CI", rep(NA_character_, 5)))|>
  #filter(!str_detect(variables, "ntercept"))|>
  ggplot(aes( Estimate, variables))+
  geom_point(shape = 16, fill = "red")+
  geom_errorbar(aes(xmin = Estimate - CI99,
                    xmax = Estimate + CI99), width = 0.2)+
  geom_errorbar(aes(xmin = Estimate - CI95,
                    xmax = Estimate + CI95), colour = "grey40", width = 0.2)+
  geom_errorbar(aes(xmin = Estimate - CI90,
                    xmax = Estimate + CI90), colour = "grey50", width = 0.2)+
  geom_text_repel(aes(x = (Estimate + CI99), label = label99), nudge_y = 0.5)+
  geom_text_repel(aes(x = (Estimate + CI95), label = label95), nudge_y = 0.5, colour = "grey40")+
  geom_text_repel(aes(x = (Estimate + CI90), label = label90), nudge_y = 0.5, colour = "grey50")+
  geom_vline(xintercept = 0)
#########################


#am div peg
##########################
boxcox2 <- MASS::boxcox(AM_divby_PEG ~ virus + sample_ph_pre_ansis + ammonia_mg_l + ophosph_mg_l + conductivity_ms_cm + Turbidity, 
                       data = diff)

boxcoxlambda2 <- boxcox2$x[which.max(boxcox2$y)]

boxcox_trans <- function(x, lambda){
  (x ^ lambda - 1) / lambda
}


lmm_bc2 <- lmer(AM_divby_PEG ~ (1|virus) + 
                 sample_ph_pre_ansis + ammonia_mg_l + ophosph_mg_l + conductivity_ms_cm + Turbidity,
               data = mutate(diff, AM_divby_PEG = boxcox_trans(AM_divby_PEG, boxcoxlambda2)))

plot_assumptions(lmm_bc2,
                 mutate(diff, AM_divby_PEG = boxcox_trans(AM_divby_PEG, boxcoxlambda2))|>
                   filter(!is.na(AM_divby_PEG))|>
                   pull(AM_divby_PEG))

#plot(lm(AM_divby_PEG ~ virus + sample_ph_pre_ansis + ammonia_mg_l + ophosph_mg_l + conductivity_ms_cm + Turbidity, 
#        data = mutate(diff, AM_divby_PEG = boxcox_trans(AM_divby_PEG, boxcoxlambda2))))

as_tibble(summary(lmm_bc2)[["coefficients"]], rownames = "variables")|>
  mutate(variables = str_replace_all(variables, c("methodPEG" = "PEG",
                                                  "methodAM" = "AM",
                                                  "sample_ph_pre_ansis" = "pH",
                                                  "ammonia_mg_l" = "Ammonia",
                                                  "ophosph_mg_l" = "Orthophosphate",
                                                  "conductivity_ms_cm" = "Conductivity")),
         CI90 = (`Std. Error`*1.65),
         CI95 = (`Std. Error`*1.96),
         CI99 = (`Std. Error`*2.59),
         label99 = c("99% CI", rep(NA_character_, 5)),
         label95 = c( "95% CI", rep(NA_character_, 5)),
         label90 = c( "90% CI", rep(NA_character_, 5)))|>
  #filter(!str_detect(variables, "ntercept"))|>
  ggplot(aes( Estimate, variables))+
  geom_point(shape = 16, fill = "red")+
  geom_errorbar(aes(xmin = Estimate - CI99,
                    xmax = Estimate + CI99), width = 0.2)+
  geom_errorbar(aes(xmin = Estimate - CI95,
                    xmax = Estimate + CI95), colour = "grey40", width = 0.2)+
  geom_errorbar(aes(xmin = Estimate - CI90,
                    xmax = Estimate + CI90), colour = "grey50", width = 0.2)+
  geom_text_repel(aes(x = (Estimate + CI99), label = label99), nudge_y = 0.5)+
  geom_text_repel(aes(x = (Estimate + CI95), label = label95), nudge_y = 0.5, colour = "grey40")+
  geom_text_repel(aes(x = (Estimate + CI90), label = label90), nudge_y = 0.5, colour = "grey50")+
  geom_vline(xintercept = 0)
#######################

qt(1 - 0.001 / 2,  (round(348/2)) - 1)



glmm <- glmer(PEG_divby_AM ~ (1|virus) + 
              sample_ph_pre_ansis + ammonia_mg_l + ophosph_mg_l + conductivity_ms_cm + Turbidity, family = Gamma(link = vlog(boxcoxlambda)),
            data = diff)

lm <- lm(PEG_divby_AM ~ 1, data = diff)



boxcox$x[which.max(boxcox$y)]

plot_assumptions(lmm, filter(diff, !is.na(PEG_divby_AM))$PEG_divby_AM)

as_tibble(summary(lmm)[["coefficients"]], rownames = "variables")|>
  mutate(variables = str_replace_all(variables, c("methodPEG" = "PEG",
                                                  "methodAM" = "AM",
                                                  "sample_ph_pre_ansis" = "pH",
                                                  "ammonia_mg_l" = "Ammonia",
                                                  "ophosph_mg_l" = "Orthophosphate",
                                                  "conductivity_ms_cm" = "Conductivity")),
         CI90 = (`Std. Error`*1.65),
         CI95 = (`Std. Error`*1.96),
         CI99 = (`Std. Error`*2.59),
         label99 = c("99% CI", rep(NA_character_, 5)),
         label95 = c( "95% CI", rep(NA_character_, 5)),
         label90 = c( "90% CI", rep(NA_character_, 5)))|>
  #filter(!str_detect(variables, "ntercept"))|>
  ggplot(aes( Estimate, variables))+
  geom_point(shape = 16, fill = "red")+
  geom_errorbar(aes(xmin = Estimate - CI99,
                     xmax = Estimate + CI99), width = 0.2)+
  geom_errorbar(aes(xmin = Estimate - CI95,
                     xmax = Estimate + CI95), colour = "grey40", width = 0.2)+
  geom_errorbar(aes(xmin = Estimate - CI90,
                     xmax = Estimate + CI90), colour = "grey50", width = 0.2)+
  geom_text_repel(aes(x = (Estimate + CI99), label = label99), nudge_y = 0.5)+
  geom_text_repel(aes(x = (Estimate + CI95), label = label95), nudge_y = 0.5, colour = "grey40")+
  geom_text_repel(aes(x = (Estimate + CI90), label = label90), nudge_y = 0.5, colour = "grey50")+
  geom_vline(xintercept = 0)



lmm2 <- glmer(AM_divby_PEG ~ (1|virus) + 
               sample_ph_pre_ansis + ammonia_mg_l + ophosph_mg_l + conductivity_ms_cm + Turbidity, family = Gamma(link = "log"),
             data = diff)

plot_assumptions(lmm2, filter(diff, !is.na(PEG_divby_AM))$PEG_divby_AM)

as_tibble(summary(lmm2)[["coefficients"]], rownames = "variables")|>
  mutate(variables = str_replace_all(variables, c("methodPEG" = "PEG",
                                                  "methodAM" = "AM",
                                                  "sample_ph_pre_ansis" = "pH",
                                                  "ammonia_mg_l" = "Ammonia",
                                                  "ophosph_mg_l" = "Orthophosphate",
                                                  "conductivity_ms_cm" = "Conductivity")),
         CI90 = (`Std. Error`*1.65),
         CI95 = (`Std. Error`*1.96),
         CI99 = (`Std. Error`*2.59),
         label99 = c("99% CI", rep(NA_character_, 5)),
         label95 = c( "95% CI", rep(NA_character_, 5)),
         label90 = c( "90% CI", rep(NA_character_, 5)))|>
  #filter(!str_detect(variables, "ntercept"))|>
  ggplot(aes( Estimate, variables))+
  geom_point(shape = 16, fill = "red")+
  geom_errorbar(aes(xmin = Estimate - CI99,
                    xmax = Estimate + CI99), width = 0.2)+
  geom_errorbar(aes(xmin = Estimate - CI95,
                    xmax = Estimate + CI95), colour = "grey40", width = 0.2)+
  geom_errorbar(aes(xmin = Estimate - CI90,
                    xmax = Estimate + CI90), colour = "grey50", width = 0.2)+
  geom_text_repel(aes(x = (Estimate + CI99), label = label99), nudge_y = 0.5)+
  geom_text_repel(aes(x = (Estimate + CI95), label = label95), nudge_y = 0.5, colour = "grey40")+
  geom_text_repel(aes(x = (Estimate + CI90), label = label90), nudge_y = 0.5, colour = "grey50")+
  geom_vline(xintercept = 0)


library(ggeffects)
ggpredict(lmm, terms = c("ammonia_mg_l"))|>
  ggplot() + 
  geom_line(aes(x = x, y = predicted)) +          # slope
  geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), 
              fill = "lightgrey", alpha = 0.5) +  # error band
  geom_point(data = diff,                      # adding the raw data (scaled values)
             aes(x = ammonia_mg_l, y = PEG_divby_AM)) + #, colour = virus
  labs(x = "Ammonia (mg/l)", y = "log10 (PEG / AM)")+
  scale_y_log10()






lmm <- lmer(gc_l_log ~ 0 + (1|virus) + (1|method) +
              method:sample_ph_pre_ansis + method:ammonia_mg_l + method:ophosph_mg_l + method:conductivity_ms_cm + method:Turbidity,
                   data = filter(long, !is.infinite(gc_l_log)), REML = FALSE)
summary(lmm)
plot_assumptions(lmm, filter(long, !is.infinite(gc_l_log), !is.na(gc_l_log))$gc_l_log)

as_tibble(summary(lmm)[["coefficients"]], rownames = "variables")|>
  mutate(variables = str_replace_all(variables, c("methodPEG" = "PEG",
                                                  "methodAM" = "AM",
                                                  "sample_ph_pre_ansis" = "pH",
                                                  "ammonia_mg_l" = "Ammonia",
                                                  "ophosph_mg_l" = "Orthophosphate",
                                                  "conductivity_ms_cm" = "Conductivity")))|>
    #filter(!str_detect(variables, "ntercept"))|>
    ggplot(aes( Estimate, variables))+
    geom_point(shape = 16, fill = "red")+
    geom_linerange(aes(xmin = Estimate - (`Std. Error`*1.96),
                       xmax = Estimate + (`Std. Error`*1.96)))+
    geom_vline(xintercept = 0)
  


qt(1 - 0.01 / 2,  (round(973/2)) - 1)



    geom_text(aes(Variables, Estimate, label = paste(round(Estimate,2), p.value)), nudge_x = 0.2)+
    ylab("Coefficients")+
    coord_flip()+
    theme_bw())

(random_effect_varience_lmm_sur <- summary(lmm_sur_sc)[["varcor"]][["method"]][[1]])
