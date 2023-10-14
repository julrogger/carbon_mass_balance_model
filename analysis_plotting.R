## Analysis and plotting of carbon mass balance model output

# Set working directory 
setwd()

# libraries
library(ggplot2)
library(deeptime)
library(RColorBrewer)
library(viridis)
library(gridExtra)
library(ggpubr)
library(cowplot)
library(dplyr)


# Import results -> Change directories to output collection
out_A <- readRDS("../Version_14_Aug23/OUTPUT_collection/model_V14E/compiled_output_A.rds")
out_B <- readRDS("../Version_14_Aug23/OUTPUT_collection/model_V14E/compiled_output_B.rds")
out_C <- readRDS("../Version_14_Aug23/OUTPUT_collection/model_V14E/compiled_output_C.rds")
out_D <- readRDS("../Version_14_Aug23/OUTPUT_collection/model_V14E/compiled_output_D.rds")
out_E <- readRDS("../Version_14_Aug23/OUTPUT_collection/model_V14E/compiled_output_E.rds")
out_F <- readRDS("../Version_14_Aug23/OUTPUT_collection/model_V14E/compiled_output_F.rds")
out_G <- readRDS("../Version_14_Aug23/OUTPUT_collection/model_V14E/compiled_output_G.rds")
out_H <- readRDS("../Version_14_Aug23/OUTPUT_collection/model_V14E/compiled_output_H.rds")
out_I <- readRDS("../Version_14_Aug23/OUTPUT_collection/model_V14E/compiled_output_I.rds")
out_J <- readRDS("../Version_14_Aug23/OUTPUT_collection/model_V14E/compiled_output_J.rds")
out_K <- readRDS("../Version_14_Aug23/OUTPUT_collection/model_V14E/compiled_output_K.rds")
out_L <- readRDS("../Version_14_Aug23/OUTPUT_collection/model_V14E/compiled_output_L.rds")
out_M <- readRDS("../Version_14_Aug23/OUTPUT_collection/model_V14E/compiled_output_M.rds")
out_N <- readRDS("../Version_14_Aug23/OUTPUT_collection/model_V14E/compiled_output_N.rds")
out_O <- readRDS("../Version_14_Aug23/OUTPUT_collection/model_V14E/compiled_output_O.rds")
out_P <- readRDS("../Version_14_Aug23/OUTPUT_collection/model_V14E/compiled_output_P.rds")
out_Q <- readRDS("../Version_14_Aug23/OUTPUT_collection/model_V14E/compiled_output_Q.rds")
out_R <- readRDS("../Version_14_Aug23/OUTPUT_collection/model_V14E/compiled_output_R.rds")
out_S <- readRDS("../Version_14_Aug23/OUTPUT_collection/model_V14E/compiled_output_S.rds")
out_T <- readRDS("../Version_14_Aug23/OUTPUT_collection/model_V14E/compiled_output_T.rds")
out_U <- readRDS("../Version_14_Aug23/OUTPUT_collection/model_V14E/compiled_output_U.rds")
out_V <- readRDS("../Version_14_Aug23/OUTPUT_collection/model_V14E/compiled_output_V.rds")
out_W <- readRDS("../Version_14_Aug23/OUTPUT_collection/model_V14E/compiled_output_W.rds")
out_X <- readRDS("../Version_14_Aug23/OUTPUT_collection/model_V14E/compiled_output_X.rds")

out <- c(out_A, out_B, out_C, out_D, out_E, out_F, out_G, out_H, out_I, out_J, out_K, out_L, out_M, out_N, out_O, out_P, out_Q, out_R, out_S, out_T, out_U, out_V, out_W, out_X) 
out_df <- do.call(rbind.data.frame, out)



# 1) Plot imbalance and effect on CO2 for non-adaptation scenario 
# Calculate hypothetical effect on atmospheric CO2 if imbalance is  maintained over duration of timestep (5e+5 years)
A_0 <- 3e+18
CO2_0 <- 280
DT <- 5e+5
# RCO2 <- (A/A0)^2
# A <- sqrt(RCO2) * A0
out_df$A_reference <- sqrt(out_df$CO2/CO2_0) * 3e+18
out_df$deltaCO2 <- (((out_df$A_reference + (out_df$dA * DT))/out_df$A_reference)^2) * out_df$CO2
out_df$deltaCO2[(out_df$A_reference + (out_df$dA * DT)) <= 0] <- 0
out_df$imbalance_per <- 1 - ((out_df$silw + out_df$sfw + out_df$locb + out_df$mocb)/(out_df$degass + out_df$oxidw))

out_df_immediate <- out_df[out_df$disp_par == 4000 & out_df$therm_adapt == 60/1e+6, ]
mean(out_df_immediate$dA)
mean(out_df_immediate$imbalance_per)
sd(out_df_immediate$imbalance_per)

out_df_slow <- out_df[out_df$disp_par == 650 & out_df$therm_adapt == 0.025/1e+6, ]
mean(out_df_slow$dA)
mean(out_df_slow$imbalance_per)
sd(out_df_slow$imbalance_per)

out_df_intermediate <- out_df[out_df$disp_par == 650 & out_df$therm_adapt == 0.05/1e+6, ]
mean(out_df_intermediate$dA)
mean(out_df_intermediate$imbalance_per)
sd(out_df_intermediate$imbalance_per)

split_list <- split(out_df, f=list(out_df$therm_adapt, out_df$disp_par))
re_out <- lapply(split_list, FUN=function(input){
  min <- aggregate(input$dA, by = list(input$time), FUN = min)$x
  max <- aggregate(input$dA, by = list(input$time), FUN = max)$x
  median <- aggregate(input$dA, by = list(input$time), FUN = median)$x
  mean <- aggregate(input$dA, by = list(input$time), FUN = mean)$x
  min_CO2rec <- aggregate(input$CO2, by = list(input$time), FUN = min)$x
  max_CO2rec <- aggregate(input$CO2, by = list(input$time), FUN = max)$x
  min_CO2hyp <- aggregate(input$deltaCO2, by = list(input$time), FUN = min)$x
  max_CO2hyp <- aggregate(input$deltaCO2, by = list(input$time), FUN = max)$x
  median_CO2hyp <- aggregate(input$deltaCO2, by = list(input$time), FUN = median)$x
  mean_CO2hyp <- aggregate(input$deltaCO2, by = list(input$time), FUN = mean)$x
  min_imbalance_per <- aggregate(input$imbalance_per, by = list(input$time), FUN = min)$x
  max_imbalance_per <- aggregate(input$imbalance_per, by = list(input$time), FUN = max)$x
  median_imbalance_per <- aggregate(input$imbalance_per, by = list(input$time), FUN = median)$x
  time <- unique(input$time)
  therm_adapt <- input$therm_adapt[1]
  disp_par <- input$disp_par[1]
  outframe <- data.frame("time" = time, "min_dA" = min, "max_dA" = max, "median_dA" = median, "mean_dA" = mean,   
                         "therm_adapt" = therm_adapt, "disp_par" = disp_par, 
                         "min_CO2hyp" = min_CO2hyp, "max_CO2hyp" = max_CO2hyp, "median_CO2hyp"= median_CO2hyp, "mean_CO2hyp" = mean_CO2hyp, 
                         "min_CO2rec" = min_CO2rec, "max_CO2rec" = max_CO2rec, 
                         "min_imbalance_per" = min_imbalance_per, "max_imbalance_per" = max_imbalance_per, "median_imbalance_per" = median_imbalance_per)
  return(outframe)
})
re_out_df <- do.call(rbind.data.frame, re_out)
re_out_df$time_geol <- re_out_df$time/(-1e+6)
re_out_df$adaptation_capacity <- interaction(re_out_df$therm_adapt, re_out_df$disp_par)


# Imalance plot percentage
temp_imbalance_immediate_per <- ggplot() +
  geom_hline(yintercept = 0, linetype = 2, col = "black") +
  geom_ribbon(data = subset(re_out_df, disp_par == 4000 & therm_adapt == 60/1e+6),
              aes(x = time_geol, ymin = min_imbalance_per*100, ymax=max_imbalance_per*100, fill = adaptation_capacity), alpha = 0.8, fill = "darkgrey") +
  theme_bw() +
  ylab(expression(Carbon ~ Imabalance ~ "(%)"))  +
  coord_geo(xlim=c(390, 0), ylim=c(-40, 60), dat = list("periods", "eras"), pos = list("b", "b"), abbrv = list(TRUE, FALSE), alpha = 0.8, size = 1.5, lwd = 0.01, height = unit(0.4, "line")) +
  border() + theme(legend.position="left") +
  theme(text = element_text(size = 6))  +
  theme(legend.key.size = unit(0.35, "cm")) +
  scale_x_reverse() +
  theme(plot.margin = unit(c(0.25, 0.25, 0, 0), "cm")) +
  theme(legend.position = "left", legend.box.spacing = unit(0, "pt")) +
  xlim(c(390, 0)) + 
  xlab("Age (Ma)") + 
  theme(axis.title.x = element_text(vjust= 1))
temp_imbalance_immediate_per

# supplementary imalance plot percentage for best balanced scenarios
lab2 <- expression(0.025 ~ "°"*C ~ Myr^-1)
lab3 <- expression(0.05 ~ "°"*C ~ Myr^-1)
temp_imbalance_bestbalanced_per <- ggplot() +
  geom_hline(yintercept = 0, linetype = 2, col = "black") +
  geom_ribbon(data = subset(re_out_df, disp_par == 650 & therm_adapt == 0.025/1e+6),
              aes(x = time_geol, ymin = min_imbalance_per*100, ymax=max_imbalance_per*100, fill = adaptation_capacity), alpha = 0.5) +
  geom_ribbon(data = subset(re_out_df, disp_par == 650 & therm_adapt == 0.05/1e+6),
              aes(x = time_geol, ymin = min_imbalance_per*100, ymax=max_imbalance_per*100, fill = adaptation_capacity), alpha = 0.5) +
  scale_fill_manual("Thermal adaptation\ncapacity", values = c("darkgrey", "#FDE725FF"), labels = c(lab2, lab3)) +
  theme_bw() +
  ylab(expression(Carbon ~ Imabalance ~ "(%)"))  +
  coord_geo(xlim=c(390, 0), ylim=c(-65, 65), dat = list("periods", "eras"), pos = list("b", "b"), abbrv = list(TRUE, FALSE), alpha = 0.8, size = 1.5, lwd = 0.01, height = unit(0.4, "line")) +
  border() + theme(legend.position="left") +
  theme(text = element_text(size = 6))  +
  theme(legend.key.size = unit(0.35, "cm")) +
  scale_x_reverse("Age (Ma)",expand=c(0.2, 0)) +
  theme(plot.margin = unit(c(0.25, 0.25, 0, 0), "cm")) +
  theme(legend.position = "left", legend.box.spacing = unit(0, "pt")) +
  xlim(c(390, 0)) +
  xlab("Age (Ma)") +
  theme(axis.title.x = element_text(vjust= 1))
temp_imbalance_bestbalanced_per


# 2) Plot the temporal imbalance for different adaptation capacities --------------------------------------------------------------------------------------------------------------------------------------
# Sort by adaptation capacity (therm + disp) --> stats over all climate and degassing trajectories
split_list <- split(out_df, f=list(out_df$therm_adapt, out_df$disp_par))
re_out <- lapply(split_list, FUN=function(input){
  min <- aggregate(input$dA, by = list(input$time), FUN = min)$x
  max <- aggregate(input$dA, by = list(input$time), FUN = max)$x
  median <- aggregate(input$dA, by = list(input$time), FUN = median)$x
  mean <- aggregate(input$dA, by = list(input$time), FUN = mean)$x
  time <- unique(input$time)
  therm_adapt <- input$therm_adapt[1]
  disp_par <- input$disp_par[1]
  outframe <- data.frame("time" = time, "min_dA" = min, "max_dA" = max, "median_dA" = median, "mean_dA" = mean,   "therm_adapt" = therm_adapt, "disp_par" = disp_par)
  return(outframe)
})

ylab <- expression(Carbon ~ Imabalance ~ "("*mol ~ C ~ yr^-1*")")
lab1 <- expression(0.01 ~ "°"*C ~ Myr^-1)
lab2 <- expression(0.025 ~ "°"*C ~ Myr^-1)
lab3 <- expression(0.05 ~ "°"*C ~ Myr^-1)
lab4 <- expression(0.1 ~ "°"*C ~ Myr^-1)
lab5 <- expression(0.25 ~ "°"*C ~ Myr^-1)
lab6 <- expression(0.5 ~ "°"*C ~ Myr^-1)
lab7 <- expression(1 ~ "°"*C ~ Myr^-1)
lab8 <- "Immediate"
re_out_df <- do.call(rbind.data.frame, re_out)
re_out_df$time_geol <- re_out_df$time/(-1e+6)
re_out_df$adaptation_capacity <- interaction(re_out_df$therm_adapt, re_out_df$disp_par)
temp_imbalance <- ggplot(subset(re_out_df, disp_par == 650 & therm_adapt < 6.0e-05 | (disp_par == 4000 & therm_adapt == 6.0e-05))) + 
  geom_ribbon(aes(x = time_geol, ymin = min_dA, ymax=max_dA, fill = as.factor(adaptation_capacity)), alpha = 1) +
  theme_bw() + 
  geom_hline(yintercept = 0, linetype = 2, col = "darkgrey") + 
  ylab(ylab)  +
  coord_geo(xlim=c(390, 0), ylim=c(-2.75e+13, 2.5e+13), dat = list("periods", "eras"), pos = list("b", "b"), abbrv = list(TRUE, FALSE), alpha = 0.8, size = 1.5, lwd = 0.01, height = unit(0.4, "line")) +  #  ylim=c(-2.5e+13, 2.5e+13)
  scale_fill_viridis("Thermal adaptation\ncapacity", option = "A", discrete = TRUE, labels = c(lab1, lab2, lab3, lab4, lab5, lab6, lab7, lab8)) + 
  scale_color_viridis("Thermal adaptation\ncapacity", option = "A", discrete = TRUE) + 
  border() + theme(legend.position="left") +
  theme(text = element_text(size = 6.25))  + 
  theme(legend.key.size = unit(0.35, "cm")) +
  scale_x_reverse("Age (Ma)",expand=c(0.2, 0)) + 
  theme(axis.title.x = element_text(vjust= 9)) + 
  theme(plot.margin = unit(c(0.25, 0.25, 0, 0), "cm")) + 
  theme(legend.position = "left", legend.box.spacing = unit(0, "pt"))
temp_imbalance


out_df$adaptation_capacity <- interaction(out_df$therm_adapt, out_df$disp_par)
out_boxplot <- subset(re_out_df, disp_par == 650 & therm_adapt < 6.0e-05 | (disp_par == 4000 & therm_adapt == 6.0e-05))
boxplot <- ggplot(out_boxplot) + 
  geom_boxplot(aes(x = adaptation_capacity, y = median_dA, fill = adaptation_capacity), outlier.shape = NA,  lwd = 0.25) + 
  theme_bw() + 
  theme(legend.position = "none") + 
  geom_hline(yintercept = 0, linetype = 2) + 
  coord_cartesian(ylim =  c(-2.05e+13, 1.5e+13)) +
  xlab("") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  scale_x_discrete(labels = c(lab1, lab2, lab3, lab4, lab5, lab6, lab7, lab8)) +
  scale_fill_viridis("Thermal adaptation\ncapacity", option = "A", discrete = TRUE) + 
  ylab("") + 
  theme(text = element_text(size = 6.25)) + 
  theme(plot.margin = unit(c(0.25, 0.25, 0, 0), "cm"))
boxplot


# 2) For model with best fit - plot single fluxes and show mass balance  --------------------------------------------------------------------------------------------------------------------------------------
median_of_medianperTS <- aggregate(re_out_df$median_dA, by=list(re_out_df$adaptation_capacity), FUN = median)
median_of_medianperTS

re_out_median <- lapply(split_list, FUN=function(input){
  med_dA <- aggregate(input$dA, by = list(input$time), FUN = median)$x
  med_degass <- aggregate(input$degass, by = list(input$time), FUN = median)$x
  med_oxidw <- aggregate(input$oxidw, by = list(input$time), FUN = median)$x
  med_silw <- aggregate(input$silw, by = list(input$time), FUN = median)$x
  med_sfw <- aggregate(input$sfw, by = list(input$time), FUN = median)$x
  med_locb <- aggregate(input$locb, by = list(input$time), FUN = median)$x
  med_mocb <- aggregate(input$mocb, by = list(input$time), FUN = median)$x
  therm_adapt <- input$therm_adapt[1]
  disp_par <- input$disp_par[1]
  time <- unique(input$time)
  outframe <- data.frame("time" = time, "med_dA" = med_dA, "med_degass" = med_degass, "med_oxidw" = med_oxidw, "med_silw" = med_silw, "med_sfw" = med_sfw, "med_locb" = med_locb, "med_mocb" = med_mocb,  "therm_adapt" = therm_adapt, "disp_par" = disp_par)
  return(outframe)
})
medians_df <- do.call(rbind.data.frame, re_out_median)
medians_df$time_geol <- medians_df$time/(-1e+6)
medians_df$adaptation_capacity <- interaction(medians_df$therm_adapt, medians_df$disp_par)

unique(medians_df$adaptation_capacity)

median_pick <- subset(medians_df, adaptation_capacity == "6e-05.4000")
median_pick$imbalance <- 1 - ((median_pick$med_silw + median_pick$med_sfw + median_pick$med_locb + median_pick$med_mocb)/(median_pick$med_degass + median_pick$med_oxidw))
mean(median_pick$imbalance)* 100
sd(median_pick$imbalance) * 100

median_pick <- subset(medians_df, adaptation_capacity == "2.5e-08.650")
median_pick$imbalance <- 1 - ((median_pick$med_silw + median_pick$med_sfw + median_pick$med_locb + median_pick$med_mocb)/(median_pick$med_degass + median_pick$med_oxidw))
mean(median_pick$imbalance) * 100
sd(median_pick$imbalance) * 100

median_pick <- subset(medians_df, adaptation_capacity == "5e-08.650")
median_pick$imbalance <- 1 - ((median_pick$med_silw + median_pick$med_sfw + median_pick$med_locb + median_pick$med_mocb)/(median_pick$med_degass + median_pick$med_oxidw))
mean(median_pick$imbalance) * 100
sd(median_pick$imbalance) * 100

ylab <- expression(Carbon ~ Fluxes ~ "("*mol ~ C ~ yr^-1*")")
oxidw <- expression(F[oxidw])
degass <- expression(F[degass])
silw <- expression(F[silw])
sfw <- expression(F[sfw])
locb <- expression(F[locb])
mocb <- expression(F[mocb])
flux_plot_slow <- ggplot(subset(medians_df,adaptation_capacity == "2.5e-08.650")) + 
  geom_ribbon(aes(x = time_geol, ymin = 0, max = med_degass, fill = "degass")) + 
  geom_ribbon(aes(x = time_geol, ymin = med_degass, ymax = med_degass + med_oxidw, fill = "oxidw")) + 
  geom_ribbon(aes(x = time_geol, ymin = -med_silw, ymax = 0, fill ="silw")) +
  geom_ribbon(aes(x = time_geol, ymin = -med_silw-med_sfw, ymax = -med_silw, fill = "sfw")) +
  geom_ribbon(aes(x = time_geol, ymin = -med_silw-med_sfw-med_locb, ymax=-med_silw-med_sfw, fill="locb")) + 
  geom_ribbon(aes(x = time_geol, ymin=-med_silw-med_sfw-med_locb-med_mocb, ymax=-med_silw-med_sfw-med_locb, fill="mocb")) + 
  geom_line(data = subset(medians_df, adaptation_capacity %in% c("2.5e-08.650")), aes(x = time_geol, y = med_dA, linetype = adaptation_capacity)) +
  scale_linetype_manual("Median C\nImbalance", values = c(1, 6, 3), 
                        labels = c("")) +
  theme_bw() + 
  scale_fill_manual("Fluxes", values = c("#B2182B" , "#FDDBC7",   "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC"), labels = c(oxidw, degass, silw, sfw, locb, mocb), breaks = c("oxidw", "degass", "silw", "sfw", "locb", "mocb")) +
  scale_color_manual("", values = c("#646F73")) +
  coord_geo(xlim=c(390, 0), ylim=c(-4.5e+13, 4.5e+13), dat = list("periods", "eras"), pos = list("b", "b"), abbrv = list(TRUE, FALSE), alpha = 0.8, size = 1.5, lwd = 0.01, height = unit(0.4, "line")) + scale_x_reverse("Age (Ma)") + 
  ylab(ylab) + 
  border() + theme(legend.position="right") + 
  theme(text = element_text(size = 6))  + 
  theme(legend.key.size = unit(0.35, "cm")) + 
  theme(legend.text.align = 0)
flux_plot_slow

flux_plot_immediate <- ggplot(subset(medians_df,adaptation_capacity == "6e-05.4000")) + 
  geom_ribbon(aes(x = time_geol, ymin = 0, max = med_degass, fill = "degass")) + 
  geom_ribbon(aes(x = time_geol, ymin = med_degass, ymax = med_degass + med_oxidw, fill = "oxidw")) + 
  geom_ribbon(aes(x = time_geol, ymin = -med_silw, ymax = 0, fill ="silw")) +
  geom_ribbon(aes(x = time_geol, ymin = -med_silw-med_sfw, ymax = -med_silw, fill = "sfw")) +
  geom_ribbon(aes(x = time_geol, ymin = -med_silw-med_sfw-med_locb, ymax=-med_silw-med_sfw, fill="locb")) + 
  geom_ribbon(aes(x = time_geol, ymin=-med_silw-med_sfw-med_locb-med_mocb, ymax=-med_silw-med_sfw-med_locb, fill="mocb")) + 
  geom_line(data = subset(medians_df, adaptation_capacity %in% c("6e-05.4000")), aes(x = time_geol, y = med_dA, linetype = adaptation_capacity)) +
  scale_linetype_manual("Median C\nImbalance", values = c(1, 6, 3), 
                        labels = c("")) +
  theme_bw() + 
  scale_fill_manual("Fluxes", values = c("#B2182B" , "#FDDBC7",   "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC"), labels = c(oxidw, degass, silw, sfw, locb, mocb), breaks = c("oxidw", "degass", "silw", "sfw", "locb", "mocb")) +
  scale_color_manual("", values = c("#646F73")) +
  coord_geo(xlim=c(390, 0), ylim=c(-4.5e+13, 4.5e+13), dat = list("periods", "eras"), pos = list("b", "b"), abbrv = list(TRUE, FALSE), alpha = 0.8, size = 1.5, lwd = 0.01, height = unit(0.4, "line")) + scale_x_reverse("Age (Ma)") + 
  ylab("") + 
  border() + theme(legend.position="right") + 
  theme(text = element_text(size = 6))  + 
  theme(legend.key.size = unit(0.35, "cm")) + 
  theme(legend.text.align = 0)
flux_plot_immediate

# Plot single fluxes - redo aggregation
re_out_2 <- lapply(split_list, FUN=function(input){
  min_dA <- aggregate(input$dA, by = list(input$time), FUN = min)$x
  max_dA <- aggregate(input$dA, by = list(input$time), FUN = max)$x
  median_dA <- aggregate(input$dA, by = list(input$time), FUN = median)$x
  min_locb <- aggregate(input$locb, by = list(input$time), FUN = min)$x
  max_locb <- aggregate(input$locb, by = list(input$time), FUN = max)$x
  median_locb <- aggregate(input$locb, by = list(input$time), FUN = median)$x
  min_silw <- aggregate(input$silw, by = list(input$time), FUN = min)$x
  max_silw <- aggregate(input$silw, by = list(input$time), FUN = max)$x
  median_silw <- aggregate(input$silw, by = list(input$time), FUN = median)$x
  min_oxidw <- aggregate(input$oxidw, by = list(input$time), FUN = min)$x
  max_oxidw <- aggregate(input$oxidw, by = list(input$time), FUN = max)$x
  median_oxidw <- aggregate(input$oxidw, by = list(input$time), FUN = median)$x
  min_mocb <- aggregate(input$mocb, by = list(input$time), FUN = min)$x
  max_mocb <- aggregate(input$mocb, by = list(input$time), FUN = max)$x
  median_mocb <- aggregate(input$mocb, by = list(input$time), FUN = median)$x
  min_sfw <- aggregate(input$sfw, by = list(input$time), FUN = min)$x
  max_sfw <- aggregate(input$sfw, by = list(input$time), FUN = max)$x
  median_sfw <- aggregate(input$sfw, by = list(input$time), FUN = median)$x
  min_degass <- aggregate(input$degass, by = list(input$time), FUN = min)$x
  max_degass <- aggregate(input$degass, by = list(input$time), FUN = max)$x
  median_degass <- aggregate(input$degass, by = list(input$time), FUN = median)$x
  min_DoA <- aggregate(input$DoA, by = list(input$time), FUN = min)$x
  max_DoA <- aggregate(input$DoA, by = list(input$time), FUN = max)$x
  median_DoA <- aggregate(input$DoA, by = list(input$time), FUN = median)$x
  min_GAST <- aggregate(input$GAST, by = list(input$time), FUN = min)$x
  max_GAST <- aggregate(input$GAST, by = list(input$time), FUN = max)$x
  median_GAST <- aggregate(input$GAST, by = list(input$time), FUN = median)$x
  min_CO2 <- aggregate(input$CO2, by = list(input$time), FUN = min)$x
  max_CO2 <- aggregate(input$CO2, by = list(input$time), FUN = max)$x
  median_CO2 <- aggregate(input$CO2, by = list(input$time), FUN = median)$x
  time <- unique(input$time)
  therm_adapt <- input$therm_adapt[1]
  disp_par <- input$disp_par[1]
  outframe <- data.frame("time" = time, "min_dA" = min_dA, "max_dA" = max_dA, "median_dA" = median_dA,
                         "min_locb" = min_locb, "max_locb" = max_locb, "median_locb"= median_locb,
                         "min_silw" = min_silw, "max_silw" = max_silw, "median_silw" = median_silw,
                         "min_oxidw" = min_oxidw, "max_oxidw" = max_oxidw, "median_oxidw" = median_oxidw,
                         "min_mocb" = min_mocb, "max_mocb" = max_mocb, "median_mocb" = median_mocb,
                         "min_sfw" = min_sfw, "max_sfw" = max_sfw, "median_sfw" = median_sfw,
                         "min_degass" = min_degass, "max_degass" = max_degass, "median_degass" = median_degass,
                         "min_DoA" = min_DoA, "max_DoA" = max_DoA, "median_DoA" = median_DoA,
                         "min_GAST" = min_GAST, "max_GAST" = max_GAST, "median_GAST" = median_GAST,
                         "min_CO2" = min_CO2, "max_CO2" = max_CO2, "median_CO2" = median_CO2,
                         "therm_adapt" = therm_adapt, "disp_par" = disp_par)
  return(outframe)
})
re_out_df <- do.call(rbind.data.frame, re_out_2)
re_out_df$time_geol <- re_out_df$time/(-1e+6)
re_out_df$adaptation_capacity <- interaction(re_out_df$therm_adapt, re_out_df$disp_par)


# Silw and locb flux plot for three scenarios
ylab <- expression(F[locb] ~ "("*mol ~ C ~ yr^-1*")")
locb_flux <- ggplot(subset(re_out_df, adaptation_capacity %in% c("2.5e-08.650", "5e-08.650", "6e-05.4000"))) + 
  geom_ribbon(aes(x = time_geol, ymin = min_locb, ymax=max_locb, fill = as.factor(adaptation_capacity)), alpha = 0.5) +
  geom_line(aes(x = time_geol, y = median_locb, col = as.factor(adaptation_capacity)), size = 0.85, alpha = 1) +
  theme_bw() + 
  ylab(ylab)  +
  coord_geo(xlim=c(390, 0), dat = list("periods", "eras"), pos = list("b", "b"), abbrv = list(TRUE, TRUE), alpha = 0.8, size = 0, lwd = 0.01, height = unit(0, "line")) + scale_x_reverse("") + 
  scale_fill_manual("Adaptation\ncapacity",
                    values = c("#231151FF", "#5F187FFF", "#FCFDBFFF"),
                    labels = c(expression(0.025 ~ "°"*C ~ Myr^-1), expression(0.05 ~ "°"*C ~ Myr^-1), "Immediate")) +
  scale_color_manual("Adaptation\ncapacity",
                     values = c("#231151FF", "#5F187FFF", "#FCFDBFFF"),
                     labels = c(expression(0.025 ~ "°"*C ~ Myr^-1), expression(0.05 ~ "°"*C ~ Myr^-1), "Immediate")) +
  border() + theme(legend.position="left") + 
  theme(text = element_text(size = 7))  + 
  theme(legend.key.size = unit(0.35, "cm"))
locb_flux

ylab <- expression(F[silw] ~ "("*mol ~ C ~ yr^-1*")")
silw_flux <- ggplot(subset(re_out_df, adaptation_capacity %in% c("2.5e-08.650", "5e-08.650", "6e-05.4000"))) + 
  geom_ribbon(aes(x = time_geol, ymin = min_silw, ymax=max_silw, fill = as.factor(adaptation_capacity)), alpha = 0.5) +
  geom_line(aes(x = time_geol, y = median_silw, col = as.factor(adaptation_capacity)), size = 0.85, alpha = 1) +
  theme_bw() + 
  ylab(ylab)  +
  coord_geo(xlim=c(390, 0), dat = list("periods", "eras"), pos = list("b", "b"), abbrv = list(TRUE, TRUE), alpha = 0.8, size = 1.5, lwd = 0.01, height = unit(0, "line")) + scale_x_reverse("Age (Ma)") + 
  scale_fill_manual("Adaptation\ncapacity",
                    values = c("#231151FF", "#5F187FFF", "#FCFDBFFF"),
                    labels = c(expression(0.025 ~ "°"*C ~ Myr^-1), expression(0.05 ~ "°"*C ~ Myr^-1), "Immediate")) +
  scale_color_manual("Adaptation\ncapacity",
                     values = c("#231151FF", "#5F187FFF", "#FCFDBFFF"),
                     labels = c(expression(0.025 ~ "°"*C ~ Myr^-1), expression(0.05 ~ "°"*C ~ Myr^-1), "Immediate")) +
  coord_geo(xlim=c(390, 0), dat = list("periods", "eras"), pos = list("b", "b"), abbrv = list(TRUE, FALSE), alpha = 0.8, size = 1.5, lwd = 0.01, height = unit(0.4, "line")) +
  border() + theme(legend.position="left") + 
  theme(text = element_text(size = 7))  + 
  theme(legend.key.size = unit(0.35, "cm"))
silw_flux


############################ Supplementary fluxes

lab1 <- expression(0.01 ~ "°"*C ~ Myr^-1)
lab2 <- expression(0.025 ~ "°"*C ~ Myr^-1)
lab3 <- expression(0.05 ~ "°"*C ~ Myr^-1)
lab4 <- expression(0.1 ~ "°"*C ~ Myr^-1)
lab5 <- expression(0.25 ~ "°"*C ~ Myr^-1)
lab6 <- expression(0.5 ~ "°"*C ~ Myr^-1)
lab7 <- expression(1 ~ "°"*C ~ Myr^-1)
lab8 <- "Immediate"
ylab <- expression(F[locb] ~ "("*mol ~ C ~ yr^-1*")")
locb_flux <- ggplot(subset(re_out_df, disp_par == 650 & therm_adapt < 6.0e-05 | (disp_par == 4000 & therm_adapt == 6.0e-05)))+ 
  geom_ribbon(aes(x = time_geol, ymin = min_locb, ymax=max_locb, fill = as.factor(adaptation_capacity)), alpha = 0.5) +
  geom_line(aes(x = time_geol, y = median_locb, col = as.factor(adaptation_capacity)), size = 2, alpha = 1) +
  theme_bw() + 
  ylab(ylab)  +
  coord_geo(xlim=c(390, 0), dat = list("periods", "eras"), pos = list("b", "b"), abbrv = list(TRUE, TRUE), alpha = 0.8, size = 1.5, lwd = 0.01, height = unit(0, "line")) + scale_x_reverse("Age (Ma)") + 
  scale_fill_viridis("Adaptation\ncapacity", option = "A", discrete = T,
                     labels = c(lab1, lab2, lab3, lab4, lab5, lab6, lab7, lab8)) +
  scale_color_viridis("Adaptation\ncapacity", option = "A", discrete = T,
                      labels = c(lab1, lab2, lab3, lab4, lab5, lab6, lab7, lab8)) +
  border() + theme(legend.position="left") + 
  theme(text = element_text(size = 7))  + 
  theme(legend.key.size = unit(0.35, "cm"))
locb_flux

ylab <- expression(F[silw] ~ "("*mol ~ C ~ yr^-1*")")
silw_flux <- ggplot(subset(re_out_df, disp_par == 650 & therm_adapt < 6.0e-05 | (disp_par == 4000 & therm_adapt == 6.0e-05)))+ 
  geom_ribbon(aes(x = time_geol, ymin = min_silw, ymax=max_silw, fill = as.factor(adaptation_capacity)), alpha = 0.5) +
  geom_line(aes(x = time_geol, y = median_silw, col = as.factor(adaptation_capacity)), size = 2, alpha = 1) +
  theme_bw() + 
  ylab(ylab)  +
  coord_geo(xlim=c(390, 0), dat = list("periods", "eras"), pos = list("b", "b"), abbrv = list(TRUE, TRUE), alpha = 0.8, size = 1.5, lwd = 0.01, height = unit(0, "line")) + scale_x_reverse("Age (Ma)") + 
  scale_fill_viridis("Adaptation\ncapacity", option = "A", discrete = T,
                    labels = c(lab1, lab2, lab3, lab4, lab5, lab6, lab7, lab8)) +
  scale_color_viridis("Adaptation\ncapacity", option = "A", discrete = T,
                     labels = c(lab1, lab2, lab3, lab4, lab5, lab6, lab7, lab8)) +
  border() + theme(legend.position="left") + 
  theme(text = element_text(size = 7))  + 
  theme(legend.key.size = unit(0.35, "cm"))
silw_flux

# The following are the same for all, not affected by eco-evolutionary dynamics
ylab <- expression(F[degass] ~ "("*mol ~ C ~ yr^-1*")")
degass_flux <- ggplot(subset(re_out_df, (disp_par == 4000 & therm_adapt == 60e-06))) + 
  geom_ribbon(aes(x = time_geol, ymin = min_degass, ymax=max_degass, fill = as.factor(adaptation_capacity)), alpha = 0.5, fill = "darkgrey") +
  geom_line(aes(x = time_geol, y = median_degass, col = as.factor(adaptation_capacity)), size = 2, alpha = 1, col = "darkgrey") +
  theme_bw() + 
  ylab(ylab)  +
  coord_geo(xlim=c(390, 0), dat = list("periods", "eras"), pos = list("b", "b"), abbrv = list(TRUE, TRUE), alpha = 0.8, size = 1.5, lwd = 0.01, height = unit(0.8, "line")) + scale_x_reverse("Age (Ma)") + 
  border() + theme(legend.position="left") + 
  theme(text = element_text(size = 7))  + 
  theme(legend.key.size = unit(0.35, "cm"))
degass_flux

ylab <- expression(F[mocb] ~ "("*mol ~ C ~ yr^-1*")")
mocb_flux <- ggplot(subset(re_out_df, (disp_par == 4000 & therm_adapt == 60e-06))) + 
  geom_ribbon(aes(x = time_geol, ymin = min_mocb, ymax=max_mocb, fill = as.factor(adaptation_capacity)), alpha = 0.5, fill = "darkgrey") +
  geom_line(aes(x = time_geol, y = median_mocb, col = as.factor(adaptation_capacity)), size = 2, alpha = 1, col = "darkgrey") +
  theme_bw() + 
  ylab(ylab)  +
  coord_geo(xlim=c(390, 0), dat = list("periods", "eras"), pos = list("b", "b"), abbrv = list(TRUE, TRUE), alpha = 0.8, size = 1.5, lwd = 0.01, height = unit(0, "line")) + scale_x_reverse("Age (Ma)") + 
  border() + theme(legend.position="left") + 
  theme(text = element_text(size = 7))  + 
  theme(legend.key.size = unit(0.35, "cm"))
mocb_flux


ylab <- expression(F[sfw] ~ "("*mol ~ C ~ yr^-1*")")
sfw_flux <- ggplot(subset(re_out_df, (disp_par == 4000 & therm_adapt == 60e-06))) + 
  geom_ribbon(aes(x = time_geol, ymin = min_sfw, ymax=max_sfw, fill = as.factor(adaptation_capacity)), alpha = 0.5, fill = "darkgrey") +
  geom_line(aes(x = time_geol, y = median_sfw, col = as.factor(adaptation_capacity)), size = 2, alpha = 1, col = "darkgrey") +
  theme_bw() + 
  ylab(ylab)  +
  coord_geo(xlim=c(390, 0), dat = list("periods", "eras"), pos = list("b", "b"), abbrv = list(TRUE, TRUE), alpha = 0.8, size = 1.5, lwd = 0.01, height = unit(0, "line")) + scale_x_reverse("Age (Ma)") + 
  border() + theme(legend.position="left") + 
  theme(text = element_text(size = 7))  + 
  theme(legend.key.size = unit(0.35, "cm"))
sfw_flux

ylab <- expression(F[oxidw] ~ "("*mol ~ C ~ yr^-1*")")
oxidw_flux <- ggplot(subset(re_out_df, (disp_par == 4000 & therm_adapt == 60e-06))) + 
  geom_ribbon(aes(x = time_geol, ymin = min_oxidw, ymax=max_oxidw, fill = as.factor(adaptation_capacity)), alpha = 0.5, fill = "darkgrey") +
  geom_line(aes(x = time_geol, y = median_oxidw, col = as.factor(adaptation_capacity)), size = 2, alpha = 1, col = "darkgrey") +
  theme_bw() + 
  ylab(ylab)  +
  coord_geo(xlim=c(390, 0), dat = list("periods", "eras"), pos = list("b", "b"), abbrv = list(TRUE, TRUE), alpha = 0.8, size = 1.5, lwd = 0.01, height = unit(0.8, "line")) + scale_x_reverse("Age (Ma)") + 
  border() + theme(legend.position="left") + 
  theme(text = element_text(size = 7))  + 
  theme(legend.key.size = unit(0.35, "cm"))
oxidw_flux


# GAST and CO2
ylab <- expression(GAST ~ "(°"*C*")")
GAST_curve <- ggplot(subset(re_out_df, (disp_par == 4000 & therm_adapt == 60e-06))) + 
  geom_ribbon(aes(x = time_geol, ymin = min_GAST, ymax=max_GAST, fill = as.factor(adaptation_capacity)), alpha = 0.5, fill = "darkgrey") +
  geom_line(aes(x = time_geol, y = median_GAST, col = as.factor(adaptation_capacity)), size = 0.75, alpha = 1, col = "darkgrey") +
  theme_bw() + 
  ylab(ylab)  +
  coord_geo(xlim=c(390, 0), dat = list("periods", "eras"), pos = list("b", "b"), abbrv = list(TRUE, FALSE), alpha = 0.8, size = 1.5, lwd = 0.01, height = unit(0.0, "line")) + scale_x_reverse("Age (Ma)") + 
  border() + theme(legend.position="left") + 
  theme(text = element_text(size = 7))  + 
  theme(legend.key.size = unit(0.35, "cm"))
GAST_curve

# DeltaT calcs
delta_GAST_low <- subset(out_df,(disp_par == 4000 & therm_adapt == 60e-06 & CO2_range == "low" & degass_range == "mid"))
delta_GAST_mid <- subset(out_df,(disp_par == 4000 & therm_adapt == 60e-06 & CO2_range == "mid" & degass_range == "mid"))
delta_GAST_high <- subset(out_df,(disp_par == 4000 & therm_adapt == 60e-06 & CO2_range == "high" & degass_range == "mid"))
delta_GAST_low$deltaT <- NA
for(i in c(2:length(delta_GAST_low$time))){
  delta_GAST_low$deltaT[i] <- abs( delta_GAST_low$GAST[i] - delta_GAST_low$GAST[i - 1] )
}
delta_GAST_mid$deltaT <- NA
for(i in c(2:length(delta_GAST_mid$time))){
  delta_GAST_mid$deltaT[i] <- abs( delta_GAST_mid$GAST[i] - delta_GAST_mid$GAST[i - 1] )
}
delta_GAST_high$deltaT <- NA
for(i in c(2:length(delta_GAST_high$time))){
  delta_GAST_high$deltaT[i] <- abs( delta_GAST_high$GAST[i] - delta_GAST_high$GAST[i - 1] )
}
delta_GAST_data <- rbind(delta_GAST_low, delta_GAST_mid, delta_GAST_high)
mean(delta_GAST_data$deltaT, na.rm=T)*2
delta_GAST_min <- aggregate(delta_GAST_data, deltaT ~ cbind(delta_GAST_data$time), FUN = min, na.rm=T)
colnames(delta_GAST_min) <- c("time", "min")
delta_GAST_max <- aggregate(delta_GAST_data, deltaT ~ cbind(delta_GAST_data$time), FUN = max, na.rm=T)
colnames(delta_GAST_max) <- c("time", "max")
delta_GAST_median <- aggregate(delta_GAST_data, deltaT ~ cbind(delta_GAST_data$time), FUN = median, na.rm=T)
colnames(delta_GAST_median) <- c("time", "median")
delta_GAST_data <- left_join(delta_GAST_min, delta_GAST_max, by = c("time"))
delta_GAST_data <- left_join(delta_GAST_data, delta_GAST_median, by = c("time"))
delta_GAST_data <- as.data.frame(delta_GAST_data)
delta_GAST_data$t_geol <- delta_GAST_data$time/(-1e+6)

deltaT_curve <- ggplot(delta_GAST_data) + 
  geom_ribbon(aes(x = t_geol, ymin = min*2, ymax=max*2), alpha = 0.5, fill = "darkgrey") +
  geom_line(aes(x = t_geol, y = median*2), size = 0.75, alpha = 1, col = "darkgrey") +
  theme_bw() + 
  ylab(expression(Delta*"GAST"~"(°C Myr"^-1*")"))  +
  coord_geo(xlim=c(390, 0), dat = list("periods", "eras"), pos = list("b", "b"), abbrv = list(TRUE, FALSE), alpha = 0.8, size = 1.5, lwd = 0.01, height = unit(0.4, "line")) + 
  scale_x_reverse("Age (Ma)") + 
  border() + theme(legend.position="left") + 
  theme(text = element_text(size = 7))  + 
  theme(legend.key.size = unit(0.35, "cm"))
deltaT_curve


ylab <- expression(CO[2] ~ "("*ppm*")")
CO2_curve <- ggplot(subset(re_out_df,  (disp_par == 4000 & therm_adapt == 60e-06))) + 
  geom_ribbon(aes(x = time_geol, ymin = min_CO2, ymax=max_CO2, fill = as.factor(adaptation_capacity)), alpha = 0.5, fill = "darkgrey") +
  geom_line(aes(x = time_geol, y = median_CO2, col = as.factor(adaptation_capacity)), size = 0.75, alpha = 1, col = "darkgrey") +
  theme_bw() + 
  ylab(ylab)  +
  coord_geo(xlim=c(390, 0), dat = list("periods", "eras"), pos = list("b", "b"), abbrv = list(TRUE, FALSE), alpha = 0.8, size = 1.5, lwd = 0.01, height = unit(0.0, "line")) + 
  scale_x_reverse("Age (Ma)") + 
  border() + theme(legend.position="left") + 
  theme(text = element_text(size = 7))  + 
  theme(legend.key.size = unit(0.35, "cm"))
CO2_curve


# Imbalance boxplots for different dispersal capacities 
re_out <- lapply(split_list, FUN=function(input){
  min <- aggregate(input$dA, by = list(input$time), FUN = min)$x
  max <- aggregate(input$dA, by = list(input$time), FUN = max)$x
  median <- aggregate(input$dA, by = list(input$time), FUN = median)$x
  mean <- aggregate(input$dA, by = list(input$time), FUN = mean)$x
  time <- unique(input$time)
  therm_adapt <- input$therm_adapt[1]
  disp_par <- input$disp_par[1]
  outframe <- data.frame("time" = time, "min_dA" = min, "max_dA" = max, "median_dA" = median, "mean_dA" = mean,   "therm_adapt" = therm_adapt, "disp_par" = disp_par)
  return(outframe)
})

ylab <- expression(Imabalance ~ "("*mol ~ C ~ yr^-1*")")
lab1 <- expression(0.01 ~ "°"*C ~ Myr^-1)
lab2 <- expression(0.025 ~ "°"*C ~ Myr^-1)
lab3 <- expression(0.05 ~ "°"*C ~ Myr^-1)
lab4 <- expression(0.1 ~ "°"*C ~ Myr^-1)
lab5 <- expression(0.25 ~ "°"*C ~ Myr^-1)
lab6 <- expression(0.5 ~ "°"*C ~ Myr^-1)
lab7 <- expression(1 ~ "°"*C ~ Myr^-1)
lab8 <- "Immediate"
re_out_df <- do.call(rbind.data.frame, re_out)
re_out_df$time_geol <- re_out_df$time/(-1e+6)
re_out_df$adaptation_capacity <- interaction(re_out_df$therm_adapt, re_out_df$disp_par)

out_df$adaptation_capacity <- interaction(out_df$therm_adapt, out_df$disp_par)
out_boxplot <- subset(re_out_df, disp_par == 550 & therm_adapt < 60e-6 | (disp_par == 4000 & therm_adapt == 60e-06))
boxplot550 <- ggplot(out_boxplot) + 
  geom_boxplot(aes(x = adaptation_capacity, y = median_dA, fill = adaptation_capacity), lwd = 0.25, outlier.size = 0.2) + 
  theme_classic() + 
  theme(legend.position = "none") + 
  geom_hline(yintercept = 0, linetype = 2) + 
  coord_cartesian(ylim=c(-2.6e+13, 1.5e+13)) +
  xlab("") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  scale_x_discrete(labels = c(lab1, lab2, lab3, lab4, lab5, lab6, lab7, lab8)) +
  scale_fill_viridis("Thermal adaptation\ncapacity", option = "A", discrete = TRUE) + 
  ylab(ylab) + 
  theme(text = element_text(size = 6)) + 
  theme(plot.margin = unit(c(0.25, 0.25, 0, 0), "cm"))
boxplot550

out_boxplot <- subset(re_out_df, disp_par == 650 & therm_adapt < 60e-6 | (disp_par == 4000 & therm_adapt == 60e-06))
boxplot650 <- ggplot(out_boxplot) + 
  geom_boxplot(aes(x = adaptation_capacity, y = median_dA, fill = adaptation_capacity), lwd = 0.25, outlier.size = 0.2) + 
  theme_classic() + 
  theme(legend.position = "none") + 
  geom_hline(yintercept = 0, linetype = 2) + 
  coord_cartesian(ylim=c(-2.6e+13, 1.5e+13)) +
  xlab("") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  scale_x_discrete(labels = c(lab1, lab2, lab3, lab4, lab5, lab6, lab7, lab8)) +
  scale_fill_viridis("Thermal adaptation\ncapacity", option = "A", discrete = TRUE) + 
  ylab("") + 
  theme(text = element_text(size = 6)) + 
  theme(plot.margin = unit(c(0.25, 0.25, 0, 0), "cm"))
boxplot650

out_boxplot <- subset(re_out_df, disp_par == 750 & therm_adapt < 60e-6 | (disp_par == 4000 & therm_adapt == 60e-06))
boxplot750 <- ggplot(out_boxplot) + 
  geom_boxplot(aes(x = adaptation_capacity, y = median_dA, fill = adaptation_capacity), lwd = 0.25, outlier.size = 0.2) + 
  theme_classic() + 
  theme(legend.position = "none") + 
  geom_hline(yintercept = 0, linetype = 2) + 
  coord_cartesian(ylim=c(-2.6e+13, 1.5e+13)) +
  xlab("") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  scale_x_discrete(labels = c(lab1, lab2, lab3, lab4, lab5, lab6, lab7, lab8)) +
  scale_fill_viridis("Thermal adaptation\ncapacity", option = "A", discrete = TRUE) + 
  ylab("") + 
  theme(text = element_text(size = 6)) + 
  theme(plot.margin = unit(c(0.25, 0.25, 0, 0), "cm"))
boxplot750

out_boxplot <- subset(re_out_df, disp_par == 850 & therm_adapt < 60e-6 | (disp_par == 4000 & therm_adapt == 60e-06))
boxplot850 <- ggplot(out_boxplot) + 
  geom_boxplot(aes(x = adaptation_capacity, y = median_dA, fill = adaptation_capacity), lwd = 0.25, outlier.size = 0.2) + 
  theme_classic() + 
  theme(legend.position = "none") + 
  geom_hline(yintercept = 0, linetype = 2) + 
  coord_cartesian(ylim=c(-2.6e+13, 1.5e+13)) + 
  xlab("") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  scale_x_discrete(labels = c(lab1, lab2, lab3, lab4, lab5, lab6, lab7, lab8)) +
  scale_fill_viridis("Thermal adaptation\ncapacity", option = "A", discrete = TRUE) + 
  ylab("") + 
  theme(text = element_text(size = 6)) + 
  theme(plot.margin = unit(c(0.25, 0.25, 0, 0), "cm"))
boxplot850


# Spatial plots  - comparison of productivity potentials modelled
##### Necessary libraries for spatial plots
library(tidyverse)
library(raster)
library(rnaturalearth)
library(sf)
library(rnaturalearthdata)
library(viridis)
library(ggpubr)
library(chronosphere)
library(scales)


robinson <- CRS("+proj=robin +over")
countries <- ne_countries(scale = 50, returnclass = c("sf"))
bb <- sf::st_union(sf::st_make_grid(
  st_bbox(c(xmin = -180,
            xmax = 180,
            ymax = 90,
            ymin = -90), crs = st_crs(4326)),
  n = 100))
bb_robinson <- st_transform(bb, as.character(robinson))
countries_robinson <- st_transform(countries, robinson)

# Choice of experiments 
unique(out_df[out_df$CO2_range == "mid" & out_df$degass_range == "mid", c("therm_adapt", "disp_par", "ID")])
# 210.1     2.5e-08      650 210
# 240.1     6.0e-05     4000 240
# 211.1     5.0e-08      650 211

det_out_immediate <- readRDS("../Version_14_Aug23/OUTPUT_collection/model_V14E/detailed_out_immediate_ID240.rds")
det_out_slow <- readRDS("../Version_14_Aug23/OUTPUT_collection/model_V14E/detailed_out_bestbalance_ID210.rds")
det_out_intermediate <- readRDS("../Version_14_Aug23/OUTPUT_collection/model_V14E/detailed_out_intermediate_ID211.rds")

det_out_slow_raster <- rasterFromXYZ(det_out_slow[[780]])
det_out_immediate_raster <- rasterFromXYZ(det_out_immediate[[780]])
det_out_intermediate_raster <- rasterFromXYZ(det_out_intermediate[[780]])
plot(det_out_slow_raster$productivity_adapt * (1/det_out_slow_raster$area) * (1/0.0007) * (12))
plot(det_out_immediate_raster$productivity_adapt * (1/det_out_immediate_raster$area) * (1/0.0007) * (12))
plot(det_out_intermediate_raster$productivity_adapt * (1/det_out_intermediate_raster$area) * (1/0.0007) * (12))

summary_out_slow <- readRDS("../Version_14_Aug23/OUTPUT_collection/model_V14E/summary_out_bestbalance_ID210.rds")
summary_out_immediate <- readRDS("../Version_14_Aug23/OUTPUT_collection/model_V14E/summary_out_immediate_ID240.rds")
summary_out_intermediate <- readRDS("../Version_14_Aug23/OUTPUT_collection/model_V14E/summary_out_intermediate_ID211.rds")

# Factor by which burial rate / (or energy conversion factor) is multiplied to obtian 3.5e+12 mol C locb per year (calibration is different for every run)
calibration_conv_factor_slow <- tail(summary_out_slow$locb, n = 1)/det_out_slow_raster$productivity_adapt_sum_global[1]
calibration_conv_factor_immediate <- tail(summary_out_immediate$locb, n = 1)/det_out_immediate_raster$productivity_adapt_sum_global[1]
calibration_conv_factor_intermediate <- tail(summary_out_intermediate$locb, n = 1)/det_out_intermediate_raster$productivity_adapt_sum_global[1]

plot(det_out_slow_raster$productivity_adapt * (1/det_out_slow_raster$area) * (1/0.0007) * (12) * calibration_conv_factor_slow)
plot(det_out_immediate_raster$productivity_adapt * (1/det_out_immediate_raster$area) * (1/0.0007) * (12) * calibration_conv_factor_immediate)
plot(det_out_intermediate_raster$productivity_adapt * (1/det_out_intermediate_raster$area) * (1/0.0007) * (12) * calibration_conv_factor_intermediate)


NPP_pot_slow <- det_out_slow_raster$weathering_limit_adapt
crs(NPP_pot_slow) <- "+proj=longlat +datum=WGS84 +no_defs"
NPP_pot_slow <- raster::projectRaster(NPP_pot_slow, crs = robinson)
NPP_pot_slow_df <- as.data.frame(NPP_pot_slow, xy = T)

# Productivity modelled - slow adaptation
title <- expression("Productivity\n& weathering potential")
p_slow <- ggplot() + 
  geom_raster(data = NPP_pot_slow_df, aes(x = x, y = y, fill = weathering_limit_adapt)) + 
  scale_fill_gradient2(title,  low = "#E5D0B1", mid = "#729E39", high = "#084113", na.value = "white", limits=c(0, 1), midpoint = 0.5, oob=squish) +
  geom_sf(data = bb_robinson, colour = 'black', linetype = 'solid', fill = NA, size = 0.2) + 
  theme_void() + 
  theme(legend.position = "bottom") + 
  geom_sf(data=countries_robinson,
          colour='#76767A',
          linetype='solid',
          fill= NA,
          size=0.1) + 
  theme(text = element_text(size = 6))
p_slow


NPP_pot_intermediate <- det_out_intermediate_raster$weathering_limit_adapt
crs(NPP_pot_intermediate) <- "+proj=longlat +datum=WGS84 +no_defs"
NPP_pot_intermediate <- raster::projectRaster(NPP_pot_intermediate, crs = robinson)
NPP_pot_intermediate_df <- as.data.frame(NPP_pot_intermediate, xy = T)

# Productivity modelled - intermediate adaptation
title <- expression("Productivity\n& weathering potential")
p_intermediate <- ggplot() + 
  geom_raster(data = NPP_pot_intermediate_df, aes(x = x, y = y, fill = weathering_limit_adapt)) + 
  scale_fill_gradient2(title,  low = "#E5D0B1", mid = "#729E39", high = "#084113", na.value = "white", limits=c(0, 1), midpoint = 0.5, oob=squish) +
  geom_sf(data = bb_robinson, colour = 'black', linetype = 'solid', fill = NA, size = 0.2) + 
  theme_void() + 
  theme(legend.position = "bottom") + 
  geom_sf(data=countries_robinson,
          colour='#76767A',
          linetype='solid',
          fill= NA,
          size=0.1) + 
  theme(text = element_text(size = 6))
p_intermediate


NPP_pot_immediate <- det_out_immediate_raster$weathering_limit_adapt
crs(NPP_pot_immediate) <- "+proj=longlat +datum=WGS84 +no_defs"
NPP_pot_immediate <- raster::projectRaster(NPP_pot_immediate, crs = robinson)
NPP_pot_immediate_df <- as.data.frame(NPP_pot_immediate, xy = T)

# Productivity modelled - immediate adaptation
title <- expression("Productivity\n& weathering potential")
p_immediate <- ggplot() + 
  geom_raster(data = NPP_pot_immediate_df, aes(x = x, y = y, fill = weathering_limit_adapt)) + 
  scale_fill_gradient2(title,  low = "#E5D0B1", mid = "#729E39", high = "#084113", na.value = "white", limits=c(0, 1), midpoint = 0.5, oob=squish) +
  geom_sf(data = bb_robinson, colour = 'black', linetype = 'solid', fill = NA, size = 0.2) + 
  theme_void() + 
  theme(legend.position = "bottom") + 
  geom_sf(data=countries_robinson,
          colour='#76767A',
          linetype='solid',
          fill= NA,
          size=0.1) + 
  theme(text = element_text(size = 6))
p_immediate



# Illustration of productivity and weathering potential for three timesteps of intermediate model 
det_out_intermediate <- readRDS("../Version_14_Aug23/OUTPUT_collection/model_V14E/detailed_out_intermediate_ID211.rds")
summary_out_intermediate <- readRDS("../Version_14_Aug23/OUTPUT_collection/model_V14E/summary_out_intermediate_ID211.rds")

# Target times: 250 Ma, 125 Ma, 50 Ma
summary_out_intermediate[summary_out_intermediate$time == -350e+6, ] # 81 timestep
summary_out_intermediate[summary_out_intermediate$time == -250e+6, ] # 281 timestep
summary_out_intermediate[summary_out_intermediate$time == -180e+6, ] # 421 timestep
summary_out_intermediate[summary_out_intermediate$time == -110e+6, ] # 561 timestep
summary_out_intermediate[summary_out_intermediate$time == -25e+6, ] # 731 timestep

first_raster <- rasterFromXYZ(det_out_intermediate[[81]])
second_raster <- rasterFromXYZ(det_out_intermediate[[281]])
third_raster <- rasterFromXYZ(det_out_intermediate[[421]])
fourth_raster <- rasterFromXYZ(det_out_intermediate[[561]])
fifth_raster <- rasterFromXYZ(det_out_intermediate[[731]])


NPP_pot_first <- first_raster$weathering_limit_adapt
crs(NPP_pot_first) <- "+proj=longlat +datum=WGS84 +no_defs"
NPP_pot_first <- raster::projectRaster(NPP_pot_first, crs = robinson)
NPP_pot_first_df <- as.data.frame(NPP_pot_first, xy = T)

title <- expression("Productivity\n& weathering potential")
p_first <- ggplot() + 
  theme(panel.background = element_rect(fill = "white")) +
  geom_sf(data = bb_robinson, colour = 'black', linetype = 'solid', fill = "lightblue", size = 0.2) +
  geom_raster(data = NPP_pot_first_df, aes(x = x, y = y, fill = weathering_limit_adapt)) + 
  scale_fill_gradient2(title,  low = "#E5D0B1", mid = "#729E39", high = "#084113", na.value = "transparent", limits=c(0, 1), midpoint = 0.5) +
  geom_sf(data = bb_robinson, colour = 'black', linetype = 'solid', fill = NA, size = 0.2) +
  theme_void() + 
  theme(legend.position = "none") + 
  theme(text = element_text(size = 6)) + 
  labs(subtitle = "350 Ma") + 
  theme(plot.subtitle = element_text(hjust = 0.25))
p_first

NPP_pot_second <- second_raster$weathering_limit_adapt
crs(NPP_pot_second) <- "+proj=longlat +datum=WGS84 +no_defs"
NPP_pot_second <- raster::projectRaster(NPP_pot_second, crs = robinson)
NPP_pot_second_df <- as.data.frame(NPP_pot_second, xy = T)

title <- expression("Productivity\n& weathering potential")
p_second <- ggplot() + 
  theme(panel.background = element_rect(fill = "white")) +
  geom_sf(data = bb_robinson, colour = 'black', linetype = 'solid', fill = "lightblue", size = 0.2) +
  geom_raster(data = NPP_pot_second_df, aes(x = x, y = y, fill = weathering_limit_adapt)) + 
  scale_fill_gradient2(title,  low = "#E5D0B1", mid = "#729E39", high = "#084113", na.value = "transparent", limits=c(0, 1), midpoint = 0.5) +
  geom_sf(data = bb_robinson, colour = 'black', linetype = 'solid', fill = NA, size = 0.2) +
  theme_void() + 
  theme(legend.position = "none") + 
  theme(text = element_text(size = 6)) + 
  labs(subtitle = "250 Ma") + 
  theme(plot.subtitle = element_text(hjust = 0.25))
p_second

NPP_pot_third <- third_raster$weathering_limit_adapt
crs(NPP_pot_third) <- "+proj=longlat +datum=WGS84 +no_defs"
NPP_pot_third <- raster::projectRaster(NPP_pot_third, crs = robinson)
NPP_pot_third_df <- as.data.frame(NPP_pot_third, xy = T)

title <- expression("Productivity\n& weathering potential")
p_third <- ggplot() + 
  theme(panel.background = element_rect(fill = "white")) +
  geom_sf(data = bb_robinson, colour = 'black', linetype = 'solid', fill = "lightblue", size = 0.2) +
  geom_raster(data = NPP_pot_third_df, aes(x = x, y = y, fill = weathering_limit_adapt)) + 
  scale_fill_gradient2(title,  low = "#E5D0B1", mid = "#729E39", high = "#084113", na.value = "transparent", limits=c(0, 1), midpoint = 0.5, oob=squish) +
  geom_sf(data = bb_robinson, colour = 'black', linetype = 'solid', fill = NA, size = 0.2) +
  theme_void() + 
  theme(legend.position = "none") + 
  theme(text = element_text(size = 6))  + 
  labs(subtitle = "180 Ma") + 
  theme(plot.subtitle = element_text(hjust = 0.25))
p_third

NPP_pot_fourth <- fourth_raster$weathering_limit_adapt
crs(NPP_pot_fourth) <- "+proj=longlat +datum=WGS84 +no_defs"
NPP_pot_fourth <- raster::projectRaster(NPP_pot_fourth, crs = robinson)
NPP_pot_fourth_df <- as.data.frame(NPP_pot_fourth, xy = T)

title <- expression("Productivity\n& weathering potential")
p_fourth <- ggplot() + 
  theme(panel.background = element_rect(fill = "white")) +
  geom_sf(data = bb_robinson, colour = 'black', linetype = 'solid', fill = "lightblue", size = 0.2) +
  geom_raster(data = NPP_pot_fourth_df, aes(x = x, y = y, fill = weathering_limit_adapt)) + 
  scale_fill_gradient2(title,  low = "#E5D0B1", mid = "#729E39", high = "#084113", na.value = "transparent", limits=c(0, 1), midpoint = 0.5) +
  geom_sf(data = bb_robinson, colour = 'black', linetype = 'solid', fill = NA, size = 0.2) +
  theme_void() + 
  theme(legend.position = "none") + 
  # geom_sf(data=countries_robinson,
  # colour='#76767A',
  # linetype='solid',
  # fill= NA,
  # size=0.1) + 
  theme(text = element_text(size = 6))  + 
  labs(subtitle = "110 Ma") + 
  theme(plot.subtitle = element_text(hjust = 0.25))
p_fourth

NPP_pot_fifth <- fifth_raster$weathering_limit_adapt
crs(NPP_pot_fifth) <- "+proj=longlat +datum=WGS84 +no_defs"
NPP_pot_fifth <- raster::projectRaster(NPP_pot_fifth, crs = robinson)
NPP_pot_fifth_df <- as.data.frame(NPP_pot_fifth, xy = T)

title <- expression("Productivity\n& weathering potential")
p_fifth <- ggplot() + 
  theme(panel.background = element_rect(fill = "white")) +
  geom_sf(data = bb_robinson, colour = 'black', linetype = 'solid', fill = "lightblue", size = 0.2) +
  geom_raster(data = NPP_pot_fifth_df, aes(x = x, y = y, fill = weathering_limit_adapt)) + 
  scale_fill_gradient2(title,  low = "#E5D0B1", mid = "#729E39", high = "#084113", na.value = "transparent", limits=c(0, 1), midpoint = 0.5) +
  geom_sf(data = bb_robinson, colour = 'black', linetype = 'solid', fill = NA, size = 0.2) +
  theme_void() + 
  theme(legend.position = "none") + 
  theme(text = element_text(size = 6)) + 
  labs(subtitle = "25 Ma") + 
  theme(plot.subtitle = element_text(hjust = 0.25))
p_fifth


# Climate plots
temp  <- (det_out_immediate_raster$temp)
crs(temp) <- "+proj=longlat +datum=WGS84 +no_defs"
temp <- raster::projectRaster(temp, crs = robinson)
temp_df <- as.data.frame(temp, xy = T)
title <- expression(Temperature ~ "(°"*C*")")
temp <- ggplot() + 
  geom_raster(data = temp_df, aes(x = x, y = y, fill = temp)) + 
  scale_fill_viridis(title, option = "H", na.value = "white", limits = c(-40, 40))+
  geom_sf(data = bb_robinson, colour = 'black', linetype = 'solid', fill = NA, size = 0.2) + 
  theme_void() + 
  theme(legend.position = "bottom") + 
  geom_sf(data=countries_robinson,
          colour='#76767A',
          linetype='solid',
          fill= NA,
          size=0.1) + 
  theme(text = element_text(size = 6))
temp

# Productivity modelled - immediate adaptation
aridity  <- (det_out_immediate_raster$aridity)
crs(aridity) <- "+proj=longlat +datum=WGS84 +no_defs"
aridity <- raster::projectRaster(aridity, crs = robinson)
aridity_df <- as.data.frame(aridity, xy = T)
title <- expression(Aridity ~ "("*"-"*")")
aridity <- ggplot() + 
  geom_raster(data = aridity_df, aes(x = x, y = y, fill = aridity)) + 
  scale_fill_viridis(title, option = "G", na.value = "white")+
  geom_sf(data = bb_robinson, colour = 'black', linetype = 'solid', fill = NA, size = 0.2) + 
  theme_void() + 
  theme(legend.position = "bottom") + 
  geom_sf(data=countries_robinson,
          colour='#76767A',
          linetype='solid',
          fill= NA,
          size=0.1) + 
  theme(text = element_text(size = 6))
aridity


# Productivity modelled - immediate adaptation
radiation  <- (det_out_immediate_raster$radiation)*0.0316887646154128 # Conversion W m-2
crs(radiation) <- "+proj=longlat +datum=WGS84 +no_defs"
radiation <- raster::projectRaster(radiation, crs = robinson)
radiation_df <- as.data.frame(radiation, xy = T)
title <- expression(Radiation ~ "("*W~m^-2*")")
radiation <- ggplot() + 
  geom_raster(data = radiation_df, aes(x = x, y = y, fill = radiation)) + 
  scale_fill_viridis(title, option = "F", na.value = "white")+
  geom_sf(data = bb_robinson, colour = 'black', linetype = 'solid', fill = NA, size = 0.2) + 
  theme_void() + 
  theme(legend.position = "bottom") + 
  geom_sf(data=countries_robinson,
          colour='#76767A',
          linetype='solid',
          fill= NA,
          size=0.1) + 
  theme(text = element_text(size = 6))
radiation

crs(det_out_immediate_raster) <- "+proj=longlat +datum=WGS84 +no_defs"
silw  <- (det_out_immediate_raster$silw)
silw_area_scaled <- silw * area(silw)*1e+6
conv_factor <- 1.325e+13/sum(values(silw_area_scaled), na.rm=T)
silw_mol_C <- silw * conv_factor #[mol C m-2 yr-1]
silw_mol_C <- raster::projectRaster(silw_mol_C, crs = robinson)
silw_mol_C_df <- as.data.frame(silw_mol_C, xy = T)
title <- expression(Silicate ~ Weathering ~ "("*mol~C~m^-2~yr^-1*")")
silw <- ggplot() + 
  geom_raster(data = silw_mol_C_df, aes(x = x, y = y, fill = silw)) + 
  scale_fill_viridis(title, option = "H", na.value = "white")+
  geom_sf(data = bb_robinson, colour = 'black', linetype = 'solid', fill = NA, size = 0.2) + 
  theme_void() + 
  theme(legend.position = "bottom") + 
  geom_sf(data=countries_robinson,
          colour='#76767A',
          linetype='solid',
          fill= NA,
          size=0.1) + 
  theme(text = element_text(size = 6))
silw



# Alternative Hypothesis Models
# A) Default - reference scenario 
# Import results 
out_A <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/default/compiled_output_default_A.rds")
out_B <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/default/compiled_output_default_B.rds")
out_C <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/default/compiled_output_default_C.rds")
out_D <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/default/compiled_output_default_D.rds")
out_E <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/default/compiled_output_default_E.rds")
out_F <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/default/compiled_output_default_F.rds")
out_G <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/default/compiled_output_default_G.rds")
out_H <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/default/compiled_output_default_H.rds")
out_C <- lapply(out_C, function(input){
  output <- input[, -c(14:25)] # some cleaning for script C
  return(output)
  })
out_default <- c(out_A, out_B, out_C, out_D, out_E, out_F, out_G, out_H)
out_default_df <- do.call(rbind.data.frame, out_default)


split_list <- split(out_default_df, f=list(out_default_df$therm_adapt, out_default_df$disp_par))
re_out_default <- lapply(split_list, FUN=function(input){
  min <- aggregate(input$dA, by = list(input$time), FUN = min)$x
  max <- aggregate(input$dA, by = list(input$time), FUN = max)$x
  median <- aggregate(input$dA, by = list(input$time), FUN = median)$x
  mean <- aggregate(input$dA, by = list(input$time), FUN = mean)$x
  time <- unique(input$time)
  therm_adapt <- input$therm_adapt[1]
  disp_par <- input$disp_par[1]
  outframe <- data.frame("time" = time, "min_dA" = min, "max_dA" = max, "median_dA" = median, "mean_dA" = mean,   "therm_adapt" = therm_adapt, "disp_par" = disp_par)
  return(outframe)
})

ylab <- expression(C ~ Cycle ~ Imabalance ~ "("*mol ~ C ~ yr^-1*")")
lab1 <- expression(0.01 ~ "°"*C ~ Myr^-1)
lab2 <- expression(0.025 ~ "°"*C ~ Myr^-1)
lab3 <- expression(0.05 ~ "°"*C ~ Myr^-1)
lab4 <- expression(0.1 ~ "°"*C ~ Myr^-1)
lab5 <- expression(0.25 ~ "°"*C ~ Myr^-1)
lab6 <- expression(0.5 ~ "°"*C ~ Myr^-1)
lab7 <- expression(1 ~ "°"*C ~ Myr^-1)
lab8 <- "Immediate"
re_out_default_df <- do.call(rbind.data.frame, re_out_default)
re_out_default_df$time_geol <- re_out_default_df$time/(-1e+6)
re_out_default_df$adaptation_capacity <- interaction(re_out_default_df$therm_adapt, re_out_default_df$disp_par)
temp_imbalance <- ggplot(subset(re_out_default_df, disp_par == 650 & therm_adapt < 60e-6 | (disp_par == 4000 & therm_adapt == 60e-06))) + 
  geom_ribbon(aes(x = time_geol, ymin = min_dA, ymax=max_dA, fill = as.factor(adaptation_capacity)), alpha = 1) +
  # geom_line(aes(x = time_geol, y = median_dA, col = as.factor(adaptation_capacity)), alpha = 1) +
  theme_bw() + 
  geom_hline(yintercept = 0, linetype = 2, col = "darkgrey") + 
  ylab(ylab)  +
  coord_geo(xlim=c(390, 0), ylim=c(-2.5e+13, 2.5e+13), dat = list("periods", "eras"), pos = list("b", "b"), abbrv = list(TRUE, FALSE), alpha = 0.8, size = 1.5, lwd = 0.01, height = unit(0.4, "line")) + 
  scale_fill_viridis("Thermal adaptation\ncapacity", option = "A", discrete = TRUE, labels = c(lab1, lab2, lab3, lab4, lab5, lab6, lab7, lab8)) + 
  scale_color_viridis("Thermal adaptation\ncapacity", option = "A", discrete = TRUE) + 
  border() + theme(legend.position="left") +
  theme(text = element_text(size = 6))  + 
  theme(legend.key.size = unit(0.35, "cm")) +
  scale_x_reverse("Age [Ma]",expand=c(0.2, 0)) + 
  theme(axis.title.x = element_text(vjust= 8)) + 
  theme(plot.margin = unit(c(0.25, 0.25, 0, 0), "cm")) + 
  theme(legend.position = "left", legend.box.spacing = unit(0, "pt"))
temp_imbalance


out_default_df$adaptation_capacity <- interaction(out_default_df$therm_adapt, out_default_df$disp_par)
out_default_boxplot <- subset(re_out_default_df, disp_par == 650 & therm_adapt < 60e-6 | (disp_par == 4000 & therm_adapt == 60e-06))
boxplot_default <- ggplot(out_default_boxplot) + 
  geom_boxplot(aes(x = adaptation_capacity, y = median_dA, fill = adaptation_capacity), lwd = 0.25, outlier.size = 0.2) + 
  theme_bw() + 
  theme(legend.position = "none") + 
  geom_hline(yintercept = 0, linetype = 2) + 
  # ylim(c(-1.75e+13, 1.75e+13)) + 
  xlab("") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  scale_x_discrete(labels = c(lab1, lab2, lab3, lab4, lab5, lab6, lab7, lab8)) +
  scale_fill_viridis("Thermal adaptation\ncapacity", option = "A", discrete = TRUE) + 
  ylab("") + 
  theme(text = element_text(size = 6)) + 
  theme(plot.margin = unit(c(0.25, 0.25, 0, 0), "cm"))
boxplot_default



# B) Strong erosion locb feedback
out_A <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/high_locb_erosion_feedback/compiled_output_high_erosion_A.rds")
out_B <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/high_locb_erosion_feedback/compiled_output_high_erosion_B.rds")
out_C <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/high_locb_erosion_feedback/compiled_output_high_erosion_C.rds")
out_D <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/high_locb_erosion_feedback/compiled_output_high_erosion_D.rds")
out_E <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/high_locb_erosion_feedback/compiled_output_high_erosion_E.rds")
out_F <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/high_locb_erosion_feedback/compiled_output_high_erosion_F.rds")
out_G <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/high_locb_erosion_feedback/compiled_output_high_erosion_G.rds")
out_H <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/high_locb_erosion_feedback/compiled_output_high_erosion_H.rds")
out_highlocberosion <- c(out_A, out_B, out_C, out_D, out_E, out_F, out_G, out_H)

out_highlocberosion_df <- do.call(rbind.data.frame, out_highlocberosion)

split_list <- split(out_highlocberosion_df, f=list(out_highlocberosion_df$therm_adapt, out_highlocberosion_df$disp_par))
re_out_highlocberosion <- lapply(split_list, FUN=function(input){
  min <- aggregate(input$dA, by = list(input$time), FUN = min)$x
  max <- aggregate(input$dA, by = list(input$time), FUN = max)$x
  median <- aggregate(input$dA, by = list(input$time), FUN = median)$x
  mean <- aggregate(input$dA, by = list(input$time), FUN = mean)$x
  time <- unique(input$time)
  therm_adapt <- input$therm_adapt[1]
  disp_par <- input$disp_par[1]
  outframe <- data.frame("time" = time, "min_dA" = min, "max_dA" = max, "median_dA" = median, "mean_dA" = mean,   "therm_adapt" = therm_adapt, "disp_par" = disp_par)
  return(outframe)
})

ylab <- expression(C ~ Cycle ~ Imabalance ~ "("*mol ~ C ~ yr^-1*")")
lab1 <- expression(0.01 ~ "°"*C ~ Myr^-1)
lab2 <- expression(0.025 ~ "°"*C ~ Myr^-1)
lab3 <- expression(0.05 ~ "°"*C ~ Myr^-1)
lab4 <- expression(0.1 ~ "°"*C ~ Myr^-1)
lab5 <- expression(0.25 ~ "°"*C ~ Myr^-1)
lab6 <- expression(0.5 ~ "°"*C ~ Myr^-1)
lab7 <- expression(1 ~ "°"*C ~ Myr^-1)
lab8 <- "Immediate"
re_out_highlocberosion_df <- do.call(rbind.data.frame, re_out_highlocberosion)
re_out_highlocberosion_df$time_geol <- re_out_highlocberosion_df$time/(-1e+6)
re_out_highlocberosion_df$adaptation_capacity <- interaction(re_out_highlocberosion_df$therm_adapt, re_out_highlocberosion_df$disp_par)
temp_imbalance <- ggplot(subset(re_out_highlocberosion_df, disp_par == 650 & therm_adapt < 60e-6 | (disp_par == 4000 & therm_adapt == 60e-06))) + 
  geom_ribbon(aes(x = time_geol, ymin = min_dA, ymax=max_dA, fill = as.factor(adaptation_capacity)), alpha = 1) +
  # geom_line(aes(x = time_geol, y = median_dA, col = as.factor(adaptation_capacity)), alpha = 1) +
  theme_bw() + 
  geom_hline(yintercept = 0, linetype = 2, col = "darkgrey") + 
  ylab(ylab)  +
  coord_geo(xlim=c(390, 0), ylim=c(-2.5e+13, 2.5e+13), dat = list("periods", "eras"), pos = list("b", "b"), abbrv = list(TRUE, FALSE), alpha = 0.8, size = 1.5, lwd = 0.01, height = unit(0.4, "line")) + 
  scale_fill_viridis("Thermal adaptation\ncapacity", option = "A", discrete = TRUE, labels = c(lab1, lab2, lab3, lab4, lab5, lab6, lab7, lab8)) + 
  scale_color_viridis("Thermal adaptation\ncapacity", option = "A", discrete = TRUE) + 
  border() + theme(legend.position="left") +
  theme(text = element_text(size = 6))  + 
  theme(legend.key.size = unit(0.35, "cm")) +
  scale_x_reverse("Age [Ma]",expand=c(0.2, 0)) + 
  theme(axis.title.x = element_text(vjust= 8)) + 
  theme(plot.margin = unit(c(0.25, 0.25, 0, 0), "cm")) + 
  theme(legend.position = "left", legend.box.spacing = unit(0, "pt"))

temp_imbalance


out_highlocberosion_df$adaptation_capacity <- interaction(out_highlocberosion_df$therm_adapt, out_highlocberosion_df$disp_par)
out_highlocberosion_boxplot <- subset(re_out_highlocberosion_df, disp_par == 650 & therm_adapt < 60e-6 | (disp_par == 4000 & therm_adapt == 60e-06))
boxplot_highlocberosion <- ggplot(out_highlocberosion_boxplot) + 
  geom_boxplot(aes(x = adaptation_capacity, y = median_dA, fill = adaptation_capacity), lwd = 0.25, outlier.size = 0.2) + 
  theme_bw() + 
  theme(legend.position = "none") + 
  geom_hline(yintercept = 0, linetype = 2) + 
  # ylim(c(-1.75e+13, 1.75e+13)) + 
  xlab("") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  scale_x_discrete(labels = c(lab1, lab2, lab3, lab4, lab5, lab6, lab7, lab8)) +
  scale_fill_viridis("Thermal adaptation\ncapacity", option = "A", discrete = TRUE) + 
  ylab("") + 
  theme(text = element_text(size = 6)) + 
  theme(plot.margin = unit(c(0.25, 0.25, 0, 0), "cm"))
boxplot_highlocberosion



# C) weak erosion locb feedback
out_A <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/low_locb_erosion_feedback/compiled_output_low_erosion_A.rds")
out_B <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/low_locb_erosion_feedback/compiled_output_low_erosion_B.rds")
out_C <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/low_locb_erosion_feedback/compiled_output_low_erosion_C.rds")
out_D <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/low_locb_erosion_feedback/compiled_output_low_erosion_D.rds")
out_E <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/low_locb_erosion_feedback/compiled_output_low_erosion_E.rds")
out_F <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/low_locb_erosion_feedback/compiled_output_low_erosion_F.rds")
out_G <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/low_locb_erosion_feedback/compiled_output_low_erosion_G.rds")
out_H <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/low_locb_erosion_feedback/compiled_output_low_erosion_H.rds")
out_lowlocberosion <- c(out_A, out_B, out_C, out_D, out_E, out_F, out_G, out_H)

out_lowlocberosion_df <- do.call(rbind.data.frame, out_lowlocberosion)

split_list <- split(out_lowlocberosion_df, f=list(out_lowlocberosion_df$therm_adapt, out_lowlocberosion_df$disp_par))
re_out_lowlocberosion <- lapply(split_list, FUN=function(input){
  min <- aggregate(input$dA, by = list(input$time), FUN = min)$x
  max <- aggregate(input$dA, by = list(input$time), FUN = max)$x
  median <- aggregate(input$dA, by = list(input$time), FUN = median)$x
  mean <- aggregate(input$dA, by = list(input$time), FUN = mean)$x
  time <- unique(input$time)
  therm_adapt <- input$therm_adapt[1]
  disp_par <- input$disp_par[1]
  outframe <- data.frame("time" = time, "min_dA" = min, "max_dA" = max, "median_dA" = median, "mean_dA" = mean,   "therm_adapt" = therm_adapt, "disp_par" = disp_par)
  return(outframe)
})

ylab <- expression(C ~ Cycle ~ Imabalance ~ "("*mol ~ C ~ yr^-1*")")
lab1 <- expression(0.01 ~ "°"*C ~ Myr^-1)
lab2 <- expression(0.025 ~ "°"*C ~ Myr^-1)
lab3 <- expression(0.05 ~ "°"*C ~ Myr^-1)
lab4 <- expression(0.1 ~ "°"*C ~ Myr^-1)
lab5 <- expression(0.25 ~ "°"*C ~ Myr^-1)
lab6 <- expression(0.5 ~ "°"*C ~ Myr^-1)
lab7 <- expression(1 ~ "°"*C ~ Myr^-1)
lab8 <- "Immediate"
re_out_lowlocberosion_df <- do.call(rbind.data.frame, re_out_lowlocberosion)
re_out_lowlocberosion_df$time_geol <- re_out_lowlocberosion_df$time/(-1e+6)
re_out_lowlocberosion_df$adaptation_capacity <- interaction(re_out_lowlocberosion_df$therm_adapt, re_out_lowlocberosion_df$disp_par)
temp_imbalance <- ggplot(subset(re_out_lowlocberosion_df, disp_par == 650 & therm_adapt < 60e-6 | (disp_par == 4000 & therm_adapt == 60e-06))) + 
  geom_ribbon(aes(x = time_geol, ymin = min_dA, ymax=max_dA, fill = as.factor(adaptation_capacity)), alpha = 1) +
  # geom_line(aes(x = time_geol, y = median_dA, col = as.factor(adaptation_capacity)), alpha = 1) +
  theme_bw() + 
  geom_hline(yintercept = 0, linetype = 2, col = "darkgrey") + 
  ylab(ylab)  +
  coord_geo(xlim=c(390, 0), ylim=c(-2.5e+13, 2.5e+13), dat = list("periods", "eras"), pos = list("b", "b"), abbrv = list(TRUE, FALSE), alpha = 0.8, size = 1.5, lwd = 0.01, height = unit(0.4, "line")) + 
  scale_fill_viridis("Thermal adaptation\ncapacity", option = "A", discrete = TRUE, labels = c(lab1, lab2, lab3, lab4, lab5, lab6, lab7, lab8)) + 
  scale_color_viridis("Thermal adaptation\ncapacity", option = "A", discrete = TRUE) + 
  border() + theme(legend.position="left") +
  theme(text = element_text(size = 6))  + 
  theme(legend.key.size = unit(0.35, "cm")) +
  scale_x_reverse("Age [Ma]",expand=c(0.2, 0)) + 
  theme(axis.title.x = element_text(vjust= 8)) + 
  theme(plot.margin = unit(c(0.25, 0.25, 0, 0), "cm")) + 
  theme(legend.position = "left", legend.box.spacing = unit(0, "pt"))
temp_imbalance


out_lowlocberosion_df$adaptation_capacity <- interaction(out_lowlocberosion_df$therm_adapt, out_lowlocberosion_df$disp_par)
out_lowlocberosion_boxplot <- subset(re_out_lowlocberosion_df, disp_par == 650 & therm_adapt < 60e-6 | (disp_par == 4000 & therm_adapt == 60e-06))
boxplot_lowlocberosion <- ggplot(out_lowlocberosion_boxplot) + 
  geom_boxplot(aes(x = adaptation_capacity, y = median_dA, fill = adaptation_capacity), lwd = 0.25, outlier.size = 0.2) + 
  theme_bw() + 
  theme(legend.position = "none") + 
  geom_hline(yintercept = 0, linetype = 2) + 
  # ylim(c(-1.75e+13, 1.75e+13)) + 
  xlab("") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  scale_x_discrete(labels = c(lab1, lab2, lab3, lab4, lab5, lab6, lab7, lab8)) +
  scale_fill_viridis("Thermal adaptation\ncapacity", option = "A", discrete = TRUE) + 
  ylab("") + 
  theme(text = element_text(size = 6)) + 
  theme(plot.margin = unit(c(0.25, 0.25, 0, 0), "cm"))
boxplot_lowlocberosion


# D) low silicate weathering activation energy 
out_A <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/low_silw_clim_respons/compiled_output_lowresponse_A.rds")
out_B <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/low_silw_clim_respons/compiled_output_lowresponse_B.rds")
out_C <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/low_silw_clim_respons/compiled_output_lowresponse_C.rds")
out_D <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/low_silw_clim_respons/compiled_output_lowresponse_D.rds")
out_E <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/low_silw_clim_respons/compiled_output_lowresponse_E.rds")
out_F <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/low_silw_clim_respons/compiled_output_lowresponse_F.rds")
out_G <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/low_silw_clim_respons/compiled_output_lowresponse_G.rds")
out_H <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/low_silw_clim_respons/compiled_output_lowresponse_H.rds")
out_lowsilwclimresponse <- c(out_A, out_B, out_C, out_D, out_E, out_F, out_G, out_H)
out_lowsilwclimresponse_df <- do.call(rbind.data.frame, out_lowsilwclimresponse)



split_list <- split(out_lowsilwclimresponse_df, f=list(out_lowsilwclimresponse_df$therm_adapt, out_lowsilwclimresponse_df$disp_par))
re_out_lowsilwclimresponse <- lapply(split_list, FUN=function(input){
  min <- aggregate(input$dA, by = list(input$time), FUN = min)$x
  max <- aggregate(input$dA, by = list(input$time), FUN = max)$x
  median <- aggregate(input$dA, by = list(input$time), FUN = median)$x
  mean <- aggregate(input$dA, by = list(input$time), FUN = mean)$x
  time <- unique(input$time)
  therm_adapt <- input$therm_adapt[1]
  disp_par <- input$disp_par[1]
  outframe <- data.frame("time" = time, "min_dA" = min, "max_dA" = max, "median_dA" = median, "mean_dA" = mean,   "therm_adapt" = therm_adapt, "disp_par" = disp_par)
  return(outframe)
})

ylab <- expression(C ~ Cycle ~ Imabalance ~ "("*mol ~ C ~ yr^-1*")")
lab1 <- expression(0.01 ~ "°"*C ~ Myr^-1)
lab2 <- expression(0.025 ~ "°"*C ~ Myr^-1)
lab3 <- expression(0.05 ~ "°"*C ~ Myr^-1)
lab4 <- expression(0.1 ~ "°"*C ~ Myr^-1)
lab5 <- expression(0.25 ~ "°"*C ~ Myr^-1)
lab6 <- expression(0.5 ~ "°"*C ~ Myr^-1)
lab7 <- expression(1 ~ "°"*C ~ Myr^-1)
lab8 <- "Immediate"
re_out_lowsilwclimresponse_df <- do.call(rbind.data.frame, re_out_lowsilwclimresponse)
re_out_lowsilwclimresponse_df$time_geol <- re_out_lowsilwclimresponse_df$time/(-1e+6)
re_out_lowsilwclimresponse_df$adaptation_capacity <- interaction(re_out_lowsilwclimresponse_df$therm_adapt, re_out_lowsilwclimresponse_df$disp_par)
temp_imbalance <- ggplot(subset(re_out_lowsilwclimresponse_df, disp_par == 650 & therm_adapt < 60e-6 | (disp_par == 4000 & therm_adapt == 60e-06))) + 
  geom_ribbon(aes(x = time_geol, ymin = min_dA, ymax=max_dA, fill = as.factor(adaptation_capacity)), alpha = 1) +
  # geom_line(aes(x = time_geol, y = median_dA, col = as.factor(adaptation_capacity)), alpha = 1) +
  theme_bw() + 
  geom_hline(yintercept = 0, linetype = 2, col = "darkgrey") + 
  ylab(ylab)  +
  coord_geo(xlim=c(390, 0), ylim=c(-2.5e+13, 2.5e+13), dat = list("periods", "eras"), pos = list("b", "b"), abbrv = list(TRUE, FALSE), alpha = 0.8, size = 1.5, lwd = 0.01, height = unit(0.4, "line")) + 
  scale_fill_viridis("Thermal adaptation\ncapacity", option = "A", discrete = TRUE, labels = c(lab1, lab2, lab3, lab4, lab5, lab6, lab7, lab8)) + 
  scale_color_viridis("Thermal adaptation\ncapacity", option = "A", discrete = TRUE) + 
  border() + theme(legend.position="left") +
  theme(text = element_text(size = 6))  + 
  theme(legend.key.size = unit(0.35, "cm")) +
  scale_x_reverse("Age [Ma]",expand=c(0.2, 0)) + 
  theme(axis.title.x = element_text(vjust= 8)) + 
  theme(plot.margin = unit(c(0.25, 0.25, 0, 0), "cm")) + 
  theme(legend.position = "left", legend.box.spacing = unit(0, "pt"))
temp_imbalance


out_lowsilwclimresponse_df$adaptation_capacity <- interaction(out_lowsilwclimresponse_df$therm_adapt, out_lowsilwclimresponse_df$disp_par)
out_lowsilwclimresponse_boxplot <- subset(re_out_lowsilwclimresponse_df, disp_par == 650 & therm_adapt < 60e-6 | (disp_par == 4000 & therm_adapt == 60e-06))
boxplot_lowsilwclimresponse <- ggplot(out_lowsilwclimresponse_boxplot) + 
  geom_boxplot(aes(x = adaptation_capacity, y = median_dA, fill = adaptation_capacity), lwd = 0.25, outlier.size = 0.2) + 
  theme_bw() + 
  theme(legend.position = "none") + 
  geom_hline(yintercept = 0, linetype = 2) + 
  # ylim(c(-1.75e+13, 1.75e+13)) + 
  xlab("") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  scale_x_discrete(labels = c(lab1, lab2, lab3, lab4, lab5, lab6, lab7, lab8)) +
  scale_fill_viridis("Thermal adaptation\ncapacity", option = "A", discrete = TRUE) + 
  ylab("") + 
  theme(text = element_text(size = 6)) + 
  theme(plot.margin = unit(c(0.25, 0.25, 0, 0), "cm"))
boxplot_lowsilwclimresponse


# E) High silicate weathering activation energy 
out_A <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/high_silw_clim_respons/compiled_output_highresponse_A.rds")
out_B <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/high_silw_clim_respons/compiled_output_highresponse_B.rds")
out_C <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/high_silw_clim_respons/compiled_output_highresponse_C.rds")
out_D <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/high_silw_clim_respons/compiled_output_highresponse_D.rds")
out_E <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/high_silw_clim_respons/compiled_output_highresponse_E.rds")
out_F <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/high_silw_clim_respons/compiled_output_highresponse_F.rds")
out_G <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/high_silw_clim_respons/compiled_output_highresponse_G.rds")
out_H <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/high_silw_clim_respons/compiled_output_highresponse_H.rds")
out_highsilwclimresponse <- c(out_A, out_B, out_C, out_D, out_E, out_F, out_G, out_H)
out_highsilwclimresponse_df <- do.call(rbind.data.frame, out_highsilwclimresponse)


split_list <- split(out_highsilwclimresponse_df, f=list(out_highsilwclimresponse_df$therm_adapt, out_highsilwclimresponse_df$disp_par))
re_out_highsilwclimresponse <- lapply(split_list, FUN=function(input){
  min <- aggregate(input$dA, by = list(input$time), FUN = min)$x
  max <- aggregate(input$dA, by = list(input$time), FUN = max)$x
  median <- aggregate(input$dA, by = list(input$time), FUN = median)$x
  mean <- aggregate(input$dA, by = list(input$time), FUN = mean)$x
  time <- unique(input$time)
  therm_adapt <- input$therm_adapt[1]
  disp_par <- input$disp_par[1]
  outframe <- data.frame("time" = time, "min_dA" = min, "max_dA" = max, "median_dA" = median, "mean_dA" = mean,   "therm_adapt" = therm_adapt, "disp_par" = disp_par)
  return(outframe)
})

ylab <- expression(C ~ Cycle ~ Imabalance ~ "("*mol ~ C ~ yr^-1*")")
lab1 <- expression(0.01 ~ "°"*C ~ Myr^-1)
lab2 <- expression(0.025 ~ "°"*C ~ Myr^-1)
lab3 <- expression(0.05 ~ "°"*C ~ Myr^-1)
lab4 <- expression(0.1 ~ "°"*C ~ Myr^-1)
lab5 <- expression(0.25 ~ "°"*C ~ Myr^-1)
lab6 <- expression(0.5 ~ "°"*C ~ Myr^-1)
lab7 <- expression(1 ~ "°"*C ~ Myr^-1)
lab8 <- "Immediate"
re_out_highsilwclimresponse_df <- do.call(rbind.data.frame, re_out_highsilwclimresponse)
re_out_highsilwclimresponse_df$time_geol <- re_out_highsilwclimresponse_df$time/(-1e+6)
re_out_highsilwclimresponse_df$adaptation_capacity <- interaction(re_out_highsilwclimresponse_df$therm_adapt, re_out_highsilwclimresponse_df$disp_par)
temp_imbalance <- ggplot(subset(re_out_highsilwclimresponse_df, disp_par == 650 & therm_adapt < 60e-6 | (disp_par == 4000 & therm_adapt == 60e-06))) + 
  geom_ribbon(aes(x = time_geol, ymin = min_dA, ymax=max_dA, fill = as.factor(adaptation_capacity)), alpha = 1) +
  # geom_line(aes(x = time_geol, y = median_dA, col = as.factor(adaptation_capacity)), alpha = 1) +
  theme_bw() + 
  geom_hline(yintercept = 0, linetype = 2, col = "darkgrey") + 
  ylab(ylab)  +
  coord_geo(xlim=c(390, 0), ylim=c(-2.5e+13, 2.5e+13), dat = list("periods", "eras"), pos = list("b", "b"), abbrv = list(TRUE, FALSE), alpha = 0.8, size = 1.5, lwd = 0.01, height = unit(0.4, "line")) + 
  scale_fill_viridis("Thermal adaptation\ncapacity", option = "A", discrete = TRUE, labels = c(lab1, lab2, lab3, lab4, lab5, lab6, lab7, lab8)) + 
  scale_color_viridis("Thermal adaptation\ncapacity", option = "A", discrete = TRUE) + 
  border() + theme(legend.position="left") +
  theme(text = element_text(size = 6))  + 
  theme(legend.key.size = unit(0.35, "cm")) +
  scale_x_reverse("Age [Ma]",expand=c(0.2, 0)) + 
  theme(axis.title.x = element_text(vjust= 8)) + 
  theme(plot.margin = unit(c(0.25, 0.25, 0, 0), "cm")) + 
  theme(legend.position = "left", legend.box.spacing = unit(0, "pt"))
temp_imbalance


out_highsilwclimresponse_df$adaptation_capacity <- interaction(out_highsilwclimresponse_df$therm_adapt, out_highsilwclimresponse_df$disp_par)
out_highsilwclimresponse_boxplot <- subset(re_out_highsilwclimresponse_df, disp_par == 650 & therm_adapt < 60e-6 | (disp_par == 4000 & therm_adapt == 60e-06))
boxplot_highsilwclimresponse <- ggplot(out_highsilwclimresponse_boxplot) + 
  geom_boxplot(aes(x = adaptation_capacity, y = median_dA, fill = adaptation_capacity), lwd = 0.25, outlier.size = 0.2) + 
  theme_bw() + 
  theme(legend.position = "none") + 
  geom_hline(yintercept = 0, linetype = 2) + 
  # ylim(c(-1.75e+13, 1.75e+13)) + 
  xlab("") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  scale_x_discrete(labels = c(lab1, lab2, lab3, lab4, lab5, lab6, lab7, lab8)) +
  scale_fill_viridis("Thermal adaptation\ncapacity", option = "A", discrete = TRUE) + 
  ylab("") + 
  theme(text = element_text(size = 6)) + 
  theme(plot.margin = unit(c(0.25, 0.25, 0, 0), "cm"))
boxplot_highsilwclimresponse


# F) With CO2 fertilization effect 
out_A <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/with_co2_fert/compiled_output_co2_fert_A.rds")
out_B <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/with_co2_fert/compiled_output_co2_fert_B.rds")
out_C <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/with_co2_fert/compiled_output_co2_fert_C.rds")
out_D <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/with_co2_fert/compiled_output_co2_fert_D.rds")
out_E <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/with_co2_fert/compiled_output_co2_fert_E.rds")
out_F <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/with_co2_fert/compiled_output_co2_fert_F.rds")
out_G <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/with_co2_fert/compiled_output_co2_fert_G.rds")
out_H <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/with_co2_fert/compiled_output_co2_fert_H.rds")
out_withco2fert <- c(out_A, out_B, out_C, out_D, out_E, out_F, out_G, out_H)

out_withco2fert_df <- do.call(rbind.data.frame, out_withco2fert)

split_list <- split(out_withco2fert_df, f=list(out_withco2fert_df$therm_adapt, out_withco2fert_df$disp_par))
re_out_withco2fert <- lapply(split_list, FUN=function(input){
  min <- aggregate(input$dA, by = list(input$time), FUN = min)$x
  max <- aggregate(input$dA, by = list(input$time), FUN = max)$x
  median <- aggregate(input$dA, by = list(input$time), FUN = median)$x
  mean <- aggregate(input$dA, by = list(input$time), FUN = mean)$x
  time <- unique(input$time)
  therm_adapt <- input$therm_adapt[1]
  disp_par <- input$disp_par[1]
  outframe <- data.frame("time" = time, "min_dA" = min, "max_dA" = max, "median_dA" = median, "mean_dA" = mean,   "therm_adapt" = therm_adapt, "disp_par" = disp_par)
  return(outframe)
})

ylab <- expression(C ~ Cycle ~ Imabalance ~ "("*mol ~ C ~ yr^-1*")")
lab1 <- expression(0.01 ~ "°"*C ~ Myr^-1)
lab2 <- expression(0.025 ~ "°"*C ~ Myr^-1)
lab3 <- expression(0.05 ~ "°"*C ~ Myr^-1)
lab4 <- expression(0.1 ~ "°"*C ~ Myr^-1)
lab5 <- expression(0.25 ~ "°"*C ~ Myr^-1)
lab6 <- expression(0.5 ~ "°"*C ~ Myr^-1)
lab7 <- expression(1 ~ "°"*C ~ Myr^-1)
lab8 <- "Immediate"
re_out_withco2fert_df <- do.call(rbind.data.frame, re_out_withco2fert)
re_out_withco2fert_df$time_geol <- re_out_withco2fert_df$time/(-1e+6)
re_out_withco2fert_df$adaptation_capacity <- interaction(re_out_withco2fert_df$therm_adapt, re_out_withco2fert_df$disp_par)
temp_imbalance <- ggplot(subset(re_out_withco2fert_df, disp_par == 650 & therm_adapt < 60e-6 | (disp_par == 4000 & therm_adapt == 60e-06))) + 
  geom_ribbon(aes(x = time_geol, ymin = min_dA, ymax=max_dA, fill = as.factor(adaptation_capacity)), alpha = 1) +
  # geom_line(aes(x = time_geol, y = median_dA, col = as.factor(adaptation_capacity)), alpha = 1) +
  theme_bw() + 
  geom_hline(yintercept = 0, linetype = 2, col = "darkgrey") + 
  ylab(ylab)  +
  coord_geo(xlim=c(390, 0), ylim=c(-2.5e+13, 2.5e+13), dat = list("periods", "eras"), pos = list("b", "b"), abbrv = list(TRUE, FALSE), alpha = 0.8, size = 1.5, lwd = 0.01, height = unit(0.4, "line")) + 
  scale_fill_viridis("Thermal adaptation\ncapacity", option = "A", discrete = TRUE, labels = c(lab1, lab2, lab3, lab4, lab5, lab6, lab7, lab8)) + 
  scale_color_viridis("Thermal adaptation\ncapacity", option = "A", discrete = TRUE) + 
  border() + theme(legend.position="left") +
  theme(text = element_text(size = 6))  + 
  theme(legend.key.size = unit(0.35, "cm")) +
  scale_x_reverse("Age [Ma]",expand=c(0.2, 0)) + 
  theme(axis.title.x = element_text(vjust= 8)) + 
  theme(plot.margin = unit(c(0.25, 0.25, 0, 0), "cm")) + 
  theme(legend.position = "left", legend.box.spacing = unit(0, "pt"))
temp_imbalance


out_withco2fert_df$adaptation_capacity <- interaction(out_withco2fert_df$therm_adapt, out_withco2fert_df$disp_par)
out_withco2fert_boxplot <- subset(re_out_withco2fert_df, disp_par == 650 & therm_adapt < 60e-6 | (disp_par == 4000 & therm_adapt == 60e-06))
boxplot_withco2fert <- ggplot(out_withco2fert_boxplot) + 
  geom_boxplot(aes(x = adaptation_capacity, y = median_dA, fill = adaptation_capacity), lwd = 0.25, outlier.size = 0.2) + 
  theme_bw() + 
  theme(legend.position = "none") + 
  geom_hline(yintercept = 0, linetype = 2) + 
  xlab("") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  scale_x_discrete(labels = c(lab1, lab2, lab3, lab4, lab5, lab6, lab7, lab8)) +
  scale_fill_viridis("Thermal adaptation\ncapacity", option = "A", discrete = TRUE) + 
  ylab("") + 
  theme(text = element_text(size = 6)) + 
  theme(plot.margin = unit(c(0.25, 0.25, 0, 0), "cm"))
boxplot_withco2fert




# G) PREPLANT of 0.25
out_A <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/low_weathering_enhancement/compiled_output_default_A.rds")
out_B <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/low_weathering_enhancement/compiled_output_default_B.rds")
out_C <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/low_weathering_enhancement/compiled_output_default_C.rds")
out_D <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/low_weathering_enhancement/compiled_output_default_D.rds")
out_E <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/low_weathering_enhancement/compiled_output_default_E.rds")
out_F <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/low_weathering_enhancement/compiled_output_default_F.rds")
out_G <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/low_weathering_enhancement/compiled_output_default_G.rds")
out_H <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/low_weathering_enhancement/compiled_output_default_H.rds")
out_lowenh<- c(out_A, out_B, out_C, out_D, out_E, out_F, out_G, out_H)
out_lowenh_df <- do.call(rbind.data.frame, out_lowenh)


split_list <- split(out_lowenh_df, f=list(out_lowenh_df$therm_adapt, out_lowenh_df$disp_par))
re_out_lowenh <- lapply(split_list, FUN=function(input){
  min <- aggregate(input$dA, by = list(input$time), FUN = min)$x
  max <- aggregate(input$dA, by = list(input$time), FUN = max)$x
  median <- aggregate(input$dA, by = list(input$time), FUN = median)$x
  mean <- aggregate(input$dA, by = list(input$time), FUN = mean)$x
  time <- unique(input$time)
  therm_adapt <- input$therm_adapt[1]
  disp_par <- input$disp_par[1]
  outframe <- data.frame("time" = time, "min_dA" = min, "max_dA" = max, "median_dA" = median, "mean_dA" = mean,   "therm_adapt" = therm_adapt, "disp_par" = disp_par)
  return(outframe)
})

ylab <- expression(C ~ Cycle ~ Imabalance ~ "("*mol ~ C ~ yr^-1*")")
lab1 <- expression(0.01 ~ "°"*C ~ Myr^-1)
lab2 <- expression(0.025 ~ "°"*C ~ Myr^-1)
lab3 <- expression(0.05 ~ "°"*C ~ Myr^-1)
lab4 <- expression(0.1 ~ "°"*C ~ Myr^-1)
lab5 <- expression(0.25 ~ "°"*C ~ Myr^-1)
lab6 <- expression(0.5 ~ "°"*C ~ Myr^-1)
lab7 <- expression(1 ~ "°"*C ~ Myr^-1)
lab8 <- "Immediate"
re_out_lowenh_df <- do.call(rbind.data.frame, re_out_lowenh)
re_out_lowenh_df$time_geol <- re_out_lowenh_df$time/(-1e+6)
re_out_lowenh_df$adaptation_capacity <- interaction(re_out_lowenh_df$therm_adapt, re_out_lowenh_df$disp_par)
temp_imbalance <- ggplot(subset(re_out_lowenh_df, disp_par == 650 & therm_adapt < 60e-6 | (disp_par == 4000 & therm_adapt == 60e-06))) + 
  geom_ribbon(aes(x = time_geol, ymin = min_dA, ymax=max_dA, fill = as.factor(adaptation_capacity)), alpha = 1) +
  # geom_line(aes(x = time_geol, y = median_dA, col = as.factor(adaptation_capacity)), alpha = 1) +
  theme_bw() + 
  geom_hline(yintercept = 0, linetype = 2, col = "darkgrey") + 
  ylab(ylab)  +
  coord_geo(xlim=c(390, 0), ylim=c(-2.5e+13, 2.5e+13), dat = list("periods", "eras"), pos = list("b", "b"), abbrv = list(TRUE, FALSE), alpha = 0.8, size = 1.5, lwd = 0.01, height = unit(0.4, "line")) + 
  scale_fill_viridis("Thermal adaptation\ncapacity", option = "A", discrete = TRUE, labels = c(lab1, lab2, lab3, lab4, lab5, lab6, lab7, lab8)) + 
  scale_color_viridis("Thermal adaptation\ncapacity", option = "A", discrete = TRUE) + 
  border() + theme(legend.position="left") +
  theme(text = element_text(size = 6))  + 
  theme(legend.key.size = unit(0.35, "cm")) +
  scale_x_reverse("Age [Ma]",expand=c(0.2, 0)) + 
  theme(axis.title.x = element_text(vjust= 8)) + 
  theme(plot.margin = unit(c(0.25, 0.25, 0, 0), "cm")) + 
  theme(legend.position = "left", legend.box.spacing = unit(0, "pt"))
temp_imbalance


out_lowenh_df$adaptation_capacity <- interaction(out_lowenh_df$therm_adapt, out_lowenh_df$disp_par)
out_lowenh_boxplot <- subset(re_out_lowenh_df, disp_par == 650 & therm_adapt < 60e-6 | (disp_par == 4000 & therm_adapt == 60e-06))
boxplot_lowenh <- ggplot(out_lowenh_boxplot) + 
  geom_boxplot(aes(x = adaptation_capacity, y = median_dA, fill = adaptation_capacity), lwd = 0.25, outlier.size = 0.2) + 
  theme_bw() + 
  theme(legend.position = "none") + 
  geom_hline(yintercept = 0, linetype = 2) + 
  xlab("") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  scale_x_discrete(labels = c(lab1, lab2, lab3, lab4, lab5, lab6, lab7, lab8)) +
  scale_fill_viridis("Thermal adaptation\ncapacity", option = "A", discrete = TRUE) + 
  ylab("") + 
  theme(text = element_text(size = 6)) + 
  theme(plot.margin = unit(c(0.25, 0.25, 0, 0), "cm"))
boxplot_lowenh


# I) silicate weathering CO2 enhancement in the absence of plants
out_A <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/with_silwco2_fert/compiled_output_co2_fert_A.rds")
out_B <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/with_silwco2_fert/compiled_output_co2_fert_B.rds")
out_C <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/with_silwco2_fert/compiled_output_co2_fert_C.rds")
out_D <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/with_silwco2_fert/compiled_output_co2_fert_D.rds")
out_E <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/with_silwco2_fert/compiled_output_co2_fert_E.rds")
out_F <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/with_silwco2_fert/compiled_output_co2_fert_F.rds")
out_G <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/with_silwco2_fert/compiled_output_co2_fert_G.rds")
out_H <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/with_silwco2_fert/compiled_output_co2_fert_H.rds")
out_silwco2fert<- c(out_A, out_B, out_C, out_D, out_E, out_F, out_G, out_H)

out_silwco2fert_df <- do.call(rbind.data.frame, out_silwco2fert)


split_list <- split(out_silwco2fert_df, f=list(out_silwco2fert_df$therm_adapt, out_silwco2fert_df$disp_par))
re_out_silwco2fert <- lapply(split_list, FUN=function(input){
  min <- aggregate(input$dA, by = list(input$time), FUN = min)$x
  max <- aggregate(input$dA, by = list(input$time), FUN = max)$x
  median <- aggregate(input$dA, by = list(input$time), FUN = median)$x
  mean <- aggregate(input$dA, by = list(input$time), FUN = mean)$x
  time <- unique(input$time)
  therm_adapt <- input$therm_adapt[1]
  disp_par <- input$disp_par[1]
  outframe <- data.frame("time" = time, "min_dA" = min, "max_dA" = max, "median_dA" = median, "mean_dA" = mean,   "therm_adapt" = therm_adapt, "disp_par" = disp_par)
  return(outframe)
})

ylab <- expression(C ~ Cycle ~ Imabalance ~ "("*mol ~ C ~ yr^-1*")")
lab1 <- expression(0.01 ~ "°"*C ~ Myr^-1)
lab2 <- expression(0.025 ~ "°"*C ~ Myr^-1)
lab3 <- expression(0.05 ~ "°"*C ~ Myr^-1)
lab4 <- expression(0.1 ~ "°"*C ~ Myr^-1)
lab5 <- expression(0.25 ~ "°"*C ~ Myr^-1)
lab6 <- expression(0.5 ~ "°"*C ~ Myr^-1)
lab7 <- expression(1 ~ "°"*C ~ Myr^-1)
lab8 <- "Immediate"
re_out_silwco2fert_df <- do.call(rbind.data.frame, re_out_silwco2fert)
re_out_silwco2fert_df$time_geol <- re_out_silwco2fert_df$time/(-1e+6)
re_out_silwco2fert_df$adaptation_capacity <- interaction(re_out_silwco2fert_df$therm_adapt, re_out_silwco2fert_df$disp_par)
temp_imbalance <- ggplot(subset(re_out_silwco2fert_df, disp_par == 650 & therm_adapt < 60e-6 | (disp_par == 4000 & therm_adapt == 60e-06))) + 
  geom_ribbon(aes(x = time_geol, ymin = min_dA, ymax=max_dA, fill = as.factor(adaptation_capacity)), alpha = 1) +
  # geom_line(aes(x = time_geol, y = median_dA, col = as.factor(adaptation_capacity)), alpha = 1) +
  theme_bw() + 
  geom_hline(yintercept = 0, linetype = 2, col = "darkgrey") + 
  ylab(ylab)  +
  coord_geo(xlim=c(390, 0), ylim=c(-2.5e+13, 2.5e+13), dat = list("periods", "eras"), pos = list("b", "b"), abbrv = list(TRUE, FALSE), alpha = 0.8, size = 1.5, lwd = 0.01, height = unit(0.4, "line")) + 
  scale_fill_viridis("Thermal adaptation\ncapacity", option = "A", discrete = TRUE, labels = c(lab1, lab2, lab3, lab4, lab5, lab6, lab7, lab8)) + 
  scale_color_viridis("Thermal adaptation\ncapacity", option = "A", discrete = TRUE) + 
  border() + theme(legend.position="left") +
  theme(text = element_text(size = 6))  + 
  theme(legend.key.size = unit(0.35, "cm")) +
  scale_x_reverse("Age [Ma]",expand=c(0.2, 0)) + 
  theme(axis.title.x = element_text(vjust= 8)) + 
  theme(plot.margin = unit(c(0.25, 0.25, 0, 0), "cm")) + 
  theme(legend.position = "left", legend.box.spacing = unit(0, "pt"))

temp_imbalance


out_silwco2fert_df$adaptation_capacity <- interaction(out_silwco2fert_df$therm_adapt, out_silwco2fert_df$disp_par)
out_silwco2fert_boxplot <- subset(re_out_silwco2fert_df, disp_par == 650 & therm_adapt < 60e-6 | (disp_par == 4000 & therm_adapt == 60e-06))
boxplot_silwco2fert <- ggplot(out_silwco2fert_boxplot) + 
  geom_boxplot(aes(x = adaptation_capacity, y = median_dA, fill = adaptation_capacity), lwd = 0.25, outlier.size = 0.2) + 
  theme_bw() + 
  theme(legend.position = "none") + 
  geom_hline(yintercept = 0, linetype = 2) + 
  xlab("") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  scale_x_discrete(labels = c(lab1, lab2, lab3, lab4, lab5, lab6, lab7, lab8)) +
  scale_fill_viridis("Thermal adaptation\ncapacity", option = "A", discrete = TRUE) + 
  ylab("") + 
  theme(text = element_text(size = 6)) + 
  theme(plot.margin = unit(c(0.25, 0.25, 0, 0), "cm"))
boxplot_silwco2fert


# H) k_sensitivity
# Import results 
out_A <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/k_sensitivity/compiled_output_A.rds")
out_A <- lapply(out_A, FUN = function(input){
  input$k_sensitivity <- 0.25
  input$deltaT <- NA
  input$delta_locb <- NA
  input$delta_locb_rel <- NA
  for(i in c(2:length(input[,"GAST"]))){
    input$deltaT[i] <- abs(input$GAST[i] - input$GAST[i - 1])
    input$delta_locb[i] <- (input$locb[i] - input$locb[i - 1])
    input$delta_locb_rel[i] <- (input$locb[i] - input$locb[i - 1])/input$locb[i]
  }
  return(input)
})
out_B <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/k_sensitivity/compiled_output_B.rds")
out_B <- lapply(out_B, FUN = function(input){
  input$k_sensitivity <- 1
  input$deltaT <- NA
  input$delta_locb <- NA
  input$delta_locb_rel <- NA
  for(i in c(2:length(input[,"GAST"]))){
    input$deltaT[i] <- abs(input$GAST[i] - input$GAST[i - 1])
    input$delta_locb[i] <- (input$locb[i] - input$locb[i - 1])
    input$delta_locb_rel[i] <- (input$locb[i] - input$locb[i - 1])/input$locb[i]
  }
  return(input)
})
out_C <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/k_sensitivity/compiled_output_C.rds")
out_C <- lapply(out_C, FUN = function(input){
  input$k_sensitivity <- 0.1
  input$deltaT <- NA
  input$delta_locb <- NA
  input$delta_locb_rel <- NA
  for(i in c(2:length(input[,"GAST"]))){
    input$deltaT[i] <- abs(input$GAST[i] - input$GAST[i - 1])
    input$delta_locb[i] <- (input$locb[i] - input$locb[i - 1])
    input$delta_locb_rel[i] <- (input$locb[i] - input$locb[i - 1])/input$locb[i]
  }
  return(input)
})
out_D <- readRDS("../Version_14_Aug23/OUTPUT_collection/supp_models/k_sensitivity/compiled_output_D.rds")
out_D <- lapply(out_D, FUN = function(input){
  input$k_sensitivity <- 0.02
  input$deltaT <- NA
  input$delta_locb <- NA
  input$delta_locb_rel <- NA
  for(i in c(2:length(input[,"GAST"]))){
    input$deltaT[i] <- abs(input$GAST[i] - input$GAST[i - 1])
    input$delta_locb[i] <- (input$locb[i] - input$locb[i - 1])
    input$delta_locb_rel[i] <- (input$locb[i] - input$locb[i - 1])/input$locb[i]
  }
  return(input)
})
out_ksens<- c(out_A, out_B, out_C, out_D)
out_ksens_df <- do.call(rbind.data.frame, out_ksens)


split_list <- split(out_ksens_df, f=list(out_ksens_df$k_sensitivity))
re_out_ksens <- lapply(split_list, FUN=function(input){
  min <- aggregate(input$dA, by = list(input$time), FUN = min)$x
  max <- aggregate(input$dA, by = list(input$time), FUN = max)$x
  median <- aggregate(input$dA, by = list(input$time), FUN = median)$x
  median_deltaT <- aggregate(input$deltaT, by = list(input$time), FUN = median)$x
  median_delta_locb <- aggregate(input$delta_locb, by = list(input$time), FUN = median)$x
  median_delta_locb_rel <- aggregate(input$delta_locb_rel, by = list(input$time), FUN = median)$x
  mean <- aggregate(input$dA, by = list(input$time), FUN = mean)$x
  time <- unique(input$time)
  k_sensitivity <- input$k_sensitivity[1]
  outframe <- data.frame("time" = time, "min_dA" = min, "max_dA" = max, "median_dA" = median, "mean_dA" = mean,
                         "k_sensitivity" = k_sensitivity, 
                         "median_deltaT" = median_deltaT, "median_delta_locb" = median_delta_locb, 
                         "median_delta_locb_rel" = median_delta_locb_rel)
  return(outframe)
})
re_out_ksens_df <- do.call(rbind.data.frame,re_out_ksens)

boxplot_ksens <- ggplot(re_out_ksens_df) + 
  geom_boxplot(aes(x = as.factor(k_sensitivity), y = median_dA, fill = as.factor(k_sensitivity)), lwd = 0.25, outlier.size = 0.2) + 
  theme_bw() + 
  theme(legend.position = "none") + 
  geom_hline(yintercept = 0, linetype = 2) + 
  xlab("") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  # scale_x_discrete("Temperature sensitivity\n penalty") +
  scale_fill_viridis("", option = "B", discrete = TRUE) +
  ylab(expression(C ~ Cycle ~ Imabalance ~ "("*mol ~ C ~ yr^-1*")")) + 
  theme(text = element_text(size = 6)) + 
  theme(plot.margin = unit(c(0.25, 0.25, 0, 0), "cm"))
boxplot_ksens


unique(re_out_ksens_df$k_sensitivity)

lm_ksens <- ggplot(subset(re_out_ksens_df)) + 
  geom_point(aes(x = median_deltaT, y = median_delta_locb_rel, col = as.factor(k_sensitivity)), alpha = 0.25, size = 0.25) + 
  geom_smooth(aes(x = median_deltaT, y = median_delta_locb_rel, col = as.factor(k_sensitivity)) , method = lm) + 
  scale_color_viridis("Penalty\nparameter", option = "B", discrete = TRUE) + 
  theme_bw() + 
  ylab(expression("Relative change "~F[locb]~"(-)")) + 
  xlab(expression(Delta*"GAST"~"(°C)")) +
  theme(text = element_text(size = 6)) + 
  theme(legend.position = "right") +
  theme(legend.key.size = unit(0.35, "cm"))
lm_ksens


model1 <- lm(data = subset(re_out_ksens_df, k_sensitivity == 1), median_delta_locb_rel ~ median_deltaT)
summary(model1)
model2 <- lm(data = subset(re_out_ksens_df, k_sensitivity == 0.25), median_delta_locb_rel ~ median_deltaT)
summary(model2)
model3 <- lm(data = subset(re_out_ksens_df, k_sensitivity == 0.1), median_delta_locb_rel ~ median_deltaT)
summary(model3)
model4 <- lm(data = subset(re_out_ksens_df, k_sensitivity == 0.02), median_delta_locb_rel ~ median_deltaT)
summary(model4)

