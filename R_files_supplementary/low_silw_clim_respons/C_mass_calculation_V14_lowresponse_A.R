###############################################################################################################################################
################### Explicit calculation of carbon fluxes - re-implementation of chapter I (done in MATLAB)
###################
# Authoer: Julian Rogger
# Start date: 26.10.2022
# Log: 26.10.22 Start coding, no implementation of uncertainty, there are some differences in focal smoothing compared to MATLAB version
# 15.11.22: Restructure code - experimental design: 
#             1.) Reconstruct CENOZOIC climate with PLASIM Simulations and CO2 values of RAE: 2 Version - smooth CO2 curve and one with abrupt changes implemented 
#             2.) Run the model with the smooth curve - once with life and once without life: Gradients introduced by plants are important to obtain a balance
#             3.) Introduce eco-evolutionary adaptation - should be still in balance + explain carbon istope excursion for the cenozoic, which coincide with climate changes (Katz or so et al.) - check out caves paper et al. 
# 01.12.22: To the model setup for the entire last 400 Myrs
###############################################################################################################################################
# 

setwd(getwd())

# Define parameter space 
thermal_adapt_par_space <- c(0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 1, 60)/1e+6
dispersal_par_space <- c(650, 4000)
PREPLANT_par_space <- c(1/10)
CO2_range_space <- c("low", "mid", "high")
degass_range <- c("low", "mid", "high")

pars <- expand.grid(thermal_adapt_par_space, dispersal_par_space, PREPLANT_par_space, CO2_range_space, degass_range)
colnames(pars) <- c("thermal_adapt_par", "dispersal_par", "PREPLANT_par", "CO2_range", "degass_range")
pars$ID <- c(1:nrow(pars))
pars_list <- split(pars, seq(nrow(pars)))

pars_list_A <- pars_list[1:18]
pars_list_B <- pars_list[19:36]
pars_list_C <- pars_list[37:54]
pars_list_D <- pars_list[55:72]
pars_list_E <- pars_list[73:90]
pars_list_F <- pars_list[91:108]
pars_list_G <- pars_list[109:126]
pars_list_H <- pars_list[127:144]


# Intialise cluster
library(parallel)
cl <- makeCluster(18, outfile = "./out_report_lowresponse_A.txt", overwrite = F)

out <- parLapply(cl = cl, pars_list_A, fun = function(list_entry){
  
  if(file.exists(paste("experiment_run_",list_entry$ID,"/summary_out.rds",sep=""))){
    print(paste("This experiment has already been run:", list_entry$ID, sep = " "))
    out <- readRDS(paste("experiment_run_",list_entry$ID,"/summary_out.rds",sep=""))
  } else {
    
    print(paste("Starting with:", list_entry$ID, sep = "#"))
    # Clear every trace of potentially already initialized run 
    if(dir.exists(paste("experiment_run_",list_entry$ID,"/", sep = ""))){
      unlink(paste("experiment_run_",list_entry$ID,"/", sep = ""), recursive = TRUE)
    }

  ######## Import model parameters
  thermal_adapt_par <- list_entry$thermal_adapt_par
  dispersal_par <- list_entry$dispersal_par
  PREPLANT_par <- list_entry$PREPLANT_par
  CO2_range <- list_entry$CO2_range
  degass_range <- list_entry$degass_range
  ID <- list_entry$ID
  
  # Prefixes - Euler vs. Local 
  prefix_topo <- "../../main_simulation/scotese_maps/"
  prefix_climate <- "../../main_simulation/climate_data/"
  prefix_geochemical_input <- "../../main_simulation/geochemical_input_data/"
  
  
  ####### Activate libraries on worker
  library(raster)
  library(gen3sis)
  library(rmatio)
  library(dplyr)
  library(matrixStats)

  ####### Create experimental folder 
  experiment_directory <- paste("experiment_run_",ID,sep="")
  dir.create(experiment_directory, showWarnings = FALSE)
  dir.create(paste(experiment_directory, "/gen3sis/", sep = ""), showWarnings = FALSE)
  dir.create(paste(experiment_directory, "/gen3sis/model_output/", sep = ""), showWarnings = FALSE)
  dir.create(paste(experiment_directory, "/gen3sis/model_output/biology_configuration_file/", sep = ""), showWarnings = FALSE)
  dir.create(paste(experiment_directory, "/gen3sis/model_output/biology_configuration_file/world_trait_map/", sep = ""), showWarnings = FALSE)
  file.copy("./biology_configuration_file.R", paste(experiment_directory, "/gen3sis/", sep = ""))
  
  
  ######## Resolution 
  # Resolution has a large impact on model speed 
  template_topo <- raster(nrows=180, ncols=360, crs="+proj=longlat +datum=WGS84 +no_defs", resolution = c(1, 1))
  template <- raster(nrows=48, ncols=96, crs="+proj=longlat +datum=WGS84 +no_defs", resolution = c(3.75, 3.75))
  
  ######## Initialization of fluxes
  # Present-day fluxes [mol C yr-1]
  k_ocdeg <- 1.25e+12
  k_ccdeg <- 1.5e+13
  k_sfw <- 1.75e+12
  k_mocb <- 3.5e+12
  k_locb <- 3.5e+12
  k_oxidw <- k_mocb + k_locb - k_ocdeg # Assuming steady state at present day
  k_silw <- 1.325e+13
  k_carbw <- 8e+12 # not relevant for atmosphere-ocean carbon mass balance
  
  ######## Compile geochemical data
  degass_read_from_SCION <- read.mat(paste(prefix_geochemical_input,"combined_D_force_revised.mat", sep = ""))
  degass_data_mills <- data.frame("time" = degass_read_from_SCION$D_force_x, "D_min" = degass_read_from_SCION$D_force_min, "D_mid" = degass_read_from_SCION$D_force_mid, "D_max" = degass_read_from_SCION$D_force_max)
  degass_data_mills$time <- degass_data_mills$time * 1e+6
  degass_data_marcilly <- (read.csv(paste(prefix_geochemical_input,"spreading_rates_marcilly_etal_21.csv", sep = ""), header=TRUE, sep = ";"))
  degass_data_marcilly <- degass_data_marcilly[order(nrow(degass_data_marcilly):1),]
  degass_data_marcilly$time <- degass_data_marcilly$time * (-1e+6)

  
  O2_x <- read.mat(paste(prefix_geochemical_input, "geochem_data_2020.mat", sep = ""))$O2_x
  O2_y <- read.mat(paste(prefix_geochemical_input, "geochem_data_2020.mat", sep = ""))$O2_y
  O2_data <- data.frame("time" = NA, "min_O2" = NA, "max_O2" = NA)
  i <- 1
  for(u in seq(1, length(O2_x)-1, 2)){
    O2_data[i,"time"] <- O2_x[u]
    O2_data[i,"min_O2"] <- min(c(O2_y[u], O2_y[u+1]))
    O2_data[i,"max_O2"] <- max(c(O2_y[u], O2_y[u+1]))
    i <- i + 1
  }
  O2_data[i, 1] <- 0; O2_data[i,2] <- 21.64; O2_data[i,3] <- 21.64
  O2_data$time <- O2_data$time * (1e+6)
  O2_spline_min <- smooth.spline(O2_data$time, O2_data$min_O2, spar = 0.1)
  O2_spline_max <- smooth.spline(O2_data$time, O2_data$max_O2, spar = 0.1)
  
  
  
  ######### Produce CO2 curve 
  CO2_data <- read.csv("../../main_simulation/co2_loess_fit.csv", header=TRUE, sep = ";")
  colnames(CO2_data)
  colnames(CO2_data) <- c("Age", "pCO2_maxprob", "lw95", "lw68", "up68", "up95")
  CO2_data$time <- CO2_data$Age * (-1e+6)
  CO2_smooth_spline <- smooth.spline(CO2_data$time, CO2_data$pCO2_maxprob, spar = NULL)
  CO2_smooth_spline_high <- smooth.spline(CO2_data$time, CO2_data$up68, spar = NULL)
  CO2_smooth_spline_low <- smooth.spline(CO2_data$time, CO2_data$lw68, spar = NULL)
  
  # Initialize
  climate_model_timesteps <- seq(-400, 0, 5)*1e+6
  climate_model_co2_levels <- c(50, seq(500, 3000, 500))
  
  out <- data.frame("time" = rep(NA, length(climate_model_timesteps)), 
                    "degass" = rep(NA,  length(climate_model_timesteps)),
                    "oxidw" = rep(NA,  length(climate_model_timesteps)),
                    "silw" = rep(NA,  length(climate_model_timesteps)), 
                    "locb" = rep(NA,  length(climate_model_timesteps)), 
                    "mocb" = rep(NA,  length(climate_model_timesteps)),
                    "sfw" = rep(NA,  length(climate_model_timesteps)), 
                    "dA" = rep(NA,  length(climate_model_timesteps)), 
                    "DoA" = rep(NA,  length(climate_model_timesteps)),
                    "GAST" = rep(NA, length(climate_model_timesteps)),
                    "carbw" = rep(NA, length(climate_model_timesteps)), 
                    "erosion" = rep(NA, length(climate_model_timesteps)), 
                    "CO2" = rep(NA, length(climate_model_timesteps)))
  output_collection <- list()
  
  ######## Start of the loop
  DT <- 5e+5
  step_index <- 1
  timesteps <- seq(-390e+6, 0, DT)
  for(t in timesteps){

    
    # CO2_values at t derived from smooth spline 
    if(CO2_range == "mid"){
      CO2_current <- predict(CO2_smooth_spline, t)$y
    } else if(CO2_range == "high"){
      CO2_current <- predict(CO2_smooth_spline_high, t)$y
    } else if(CO2_range == "low"){
      CO2_current <- predict(CO2_smooth_spline_low, t)$y
    }
  
    # Time contribution 
    if(t %in% climate_model_timesteps){
      earlier_timestep <- t
      later_timestep <- t
      earlier_contribution <- 1
      later_contribution <- 0
    } else {
      later_timestep <- min(climate_model_timesteps[climate_model_timesteps > t])
      earlier_timestep <- max(climate_model_timesteps[climate_model_timesteps < t])
      earlier_contribution <- 1 - abs((earlier_timestep - t))/(later_timestep - earlier_timestep)
      later_contribution <- 1 - earlier_contribution
    }
    
    if(CO2_current %in% climate_model_co2_levels){
      lower_CO2 <- CO2_current
      upper_CO2 <- CO2_current
      CO2_contribution_lower <- 1
      CO2_contribution_upper <- 0
    } else {
      lower_CO2 <- climate_model_co2_levels[which(CO2_current - climate_model_co2_levels == min(CO2_current - climate_model_co2_levels[CO2_current - climate_model_co2_levels > 0]))]
      upper_CO2 <- climate_model_co2_levels[which(CO2_current - climate_model_co2_levels == max(CO2_current - climate_model_co2_levels[CO2_current - climate_model_co2_levels <= 0]))]
      CO2_contribution_lower <- 1 - (CO2_current - lower_CO2)/(upper_CO2 - lower_CO2)
      CO2_contribution_upper <- 1 - (upper_CO2 - CO2_current)/(upper_CO2 - lower_CO2)
    }
    
    
    # topography
    template <- raster(nrows=48, ncols=96, crs="+proj=longlat +datum=WGS84 +no_defs", resolution = c(3.75, 3.75))
    file_list <- list.files(prefix_topo, pattern=".nc$")
    index <- grep(paste("_",earlier_timestep/(-1e+6), "Ma.nc", sep = ""), file_list)
    earlier_topo <- raster(paste(prefix_topo,file_list[index], sep=""))
    index <- grep(paste("_",later_timestep/(-1e+6), "Ma.nc", sep = ""), file_list)
    later_topo <- raster(paste(prefix_topo,file_list[index], sep=""))
    topo_current <- earlier_contribution * earlier_topo + later_contribution * later_topo
    topo_current <- resample(topo_current, template, method="bilinear")
    topo_current_land <- topo_current
    topo_current_land[topo_current <= 0] <- NA
    
    
    
    # temperature 
    # Interpolation: 1.) To CO2 value 2.) Temporal 
    temp_early_lower_co2 <- raster(paste(prefix_climate,"run_",earlier_timestep/(-1e+6),"Ma/processed_",lower_CO2,"_ppm/timmean_",lower_CO2,"_ppm.nc", sep=""), var = "ts", level = 4)
    temp_early_lower_co2 <- rotate(temp_early_lower_co2)
    temp_early_higher_co2 <- raster(paste(prefix_climate,"run_",earlier_timestep/(-1e+6),"Ma/processed_",upper_CO2,"_ppm/timmean_",upper_CO2,"_ppm.nc", sep=""), var = "ts", level = 4)
    temp_early_higher_co2 <- rotate(temp_early_higher_co2)
    temp_early_current_co2 <- (CO2_contribution_lower * temp_early_lower_co2 + CO2_contribution_upper * temp_early_higher_co2) - 273.15
    
    temp_later_lower_co2 <- raster(paste(prefix_climate,"run_",later_timestep/(-1e+6),"Ma/processed_",lower_CO2,"_ppm/timmean_",lower_CO2,"_ppm.nc", sep=""), var = "ts", level = 4)
    temp_later_lower_co2 <- rotate(temp_later_lower_co2)
    temp_later_higher_co2 <- raster(paste(prefix_climate,"run_",later_timestep/(-1e+6),"Ma/processed_",upper_CO2,"_ppm/timmean_",upper_CO2,"_ppm.nc", sep=""), var = "ts", level = 4)
    temp_later_higher_co2 <- rotate(temp_later_higher_co2)
    temp_later_current_co2 <- (CO2_contribution_lower * temp_later_lower_co2 + CO2_contribution_upper * temp_later_higher_co2) - 273.15
    
    temp_current <- earlier_contribution * temp_early_current_co2 + later_contribution * temp_later_current_co2
    GAST <- weighted.mean(values(temp_current), values(area(temp_current)))
    
    ###### runoff
    precip_early_lower_co2 <- raster(paste(prefix_climate,"run_",earlier_timestep/(-1e+6),"Ma/processed_",lower_CO2,"_ppm/timmean_",lower_CO2,"_ppm.nc", sep=""), var = "pr", level = 1)
    evap_early_lower_co2 <- raster(paste(prefix_climate,"run_",earlier_timestep/(-1e+6),"Ma/processed_",lower_CO2,"_ppm/timmean_",lower_CO2,"_ppm.nc", sep=""), var = "evap", level = 1)
    runoff_early_lower_co2 <- precip_early_lower_co2 + evap_early_lower_co2
    runoff_early_lower_co2 <- rotate(runoff_early_lower_co2)
    precip_early_higher_co2 <- raster(paste(prefix_climate,"run_",earlier_timestep/(-1e+6),"Ma/processed_",upper_CO2,"_ppm/timmean_",upper_CO2,"_ppm.nc", sep=""), var = "pr", level = 1)
    evap_early_higher_co2 <- raster(paste(prefix_climate,"run_",earlier_timestep/(-1e+6),"Ma/processed_",upper_CO2,"_ppm/timmean_",upper_CO2,"_ppm.nc", sep=""), var = "evap", level = 1)
    runoff_early_higher_co2 <- precip_early_higher_co2 + evap_early_higher_co2
    runoff_early_higher_co2 <- rotate(runoff_early_higher_co2)
    runoff_early_current_co2 <- CO2_contribution_lower * runoff_early_lower_co2 + CO2_contribution_upper * runoff_early_higher_co2
    
    precip_later_lower_co2 <- raster(paste(prefix_climate,"run_",later_timestep/(-1e+6),"Ma/processed_",lower_CO2,"_ppm/timmean_",lower_CO2,"_ppm.nc", sep=""), var = "pr", level = 1)
    evap_later_lower_co2 <- raster(paste(prefix_climate,"run_",later_timestep/(-1e+6),"Ma/processed_",lower_CO2,"_ppm/timmean_",lower_CO2,"_ppm.nc", sep=""), var = "evap", level = 1)
    runoff_later_lower_co2 <- precip_later_lower_co2 + evap_later_lower_co2
    runoff_later_lower_co2 <- rotate(runoff_later_lower_co2)
    precip_later_higher_co2 <- raster(paste(prefix_climate,"run_",later_timestep/(-1e+6),"Ma/processed_",upper_CO2,"_ppm/timmean_",upper_CO2,"_ppm.nc", sep=""), var = "pr", level = 1)
    evap_later_higher_co2 <- raster(paste(prefix_climate,"run_",later_timestep/(-1e+6),"Ma/processed_",upper_CO2,"_ppm/timmean_",upper_CO2,"_ppm.nc", sep=""), var = "evap", level = 1)
    runoff_later_higher_co2 <- precip_later_higher_co2 + evap_later_higher_co2
    runoff_later_higher_co2 <- rotate(runoff_later_higher_co2)
    runoff_later_current_co2 <- CO2_contribution_lower * runoff_later_lower_co2 + CO2_contribution_upper * runoff_later_higher_co2
    
    runoff_current <- earlier_contribution * runoff_early_current_co2 + later_contribution * runoff_later_current_co2
    runoff_current[topo_current <= 0] <- NA
    runoff_current[runoff_current <= 0] <- 0
    runoff_current <- runoff_current*(365*24*60*60*1000)
    
    
    ###### Aridity --> Three raster needed: 1) Temperature, 2) Net Shortwave radiation, 3) Precipitation 
    # 1) Temperature 
    # plot(temp_raster_present)
    temp_current_K <- temp_current + 273.15 # [K]
    #lambda is latent heat of evaporation in mJ m-3 and is a function of temperature lambda = (3.146 - 0.002361 Td) * 10^3 in MJ m-3
    latent_heat <- (3.146 - 0.002361*temp_current_K) * 10^3
    
    # 2) Net Shortwave radiation [W m-2] = [J s-1 m-2]
    rss_early_lower_co2 <- raster(paste(prefix_climate,"run_",earlier_timestep/(-1e+6),"Ma/processed_",lower_CO2,"_ppm/timmean_",lower_CO2,"_ppm.nc", sep=""), var = "rss", level = 4)
    rss_early_lower_co2 <- rotate(rss_early_lower_co2)
    rss_early_higher_co2 <- raster(paste(prefix_climate,"run_",earlier_timestep/(-1e+6),"Ma/processed_",upper_CO2,"_ppm/timmean_",upper_CO2,"_ppm.nc", sep=""), var = "rss", level = 4)
    rss_early_higher_co2 <- rotate(rss_early_higher_co2)
    rss_early_current_co2 <- (CO2_contribution_lower * rss_early_lower_co2 + CO2_contribution_upper * rss_early_higher_co2)
    
    rss_later_lower_co2 <- raster(paste(prefix_climate,"run_",later_timestep/(-1e+6),"Ma/processed_",lower_CO2,"_ppm/timmean_",lower_CO2,"_ppm.nc", sep=""), var = "rss", level = 4)
    rss_later_lower_co2 <- rotate(rss_later_lower_co2)
    rss_later_higher_co2 <- raster(paste(prefix_climate,"run_",later_timestep/(-1e+6),"Ma/processed_",upper_CO2,"_ppm/timmean_",upper_CO2,"_ppm.nc", sep=""), var = "rss", level = 4)
    rss_later_higher_co2 <- rotate(rss_later_higher_co2)
    rss_later_current_co2 <- (CO2_contribution_lower * rss_later_lower_co2 + CO2_contribution_upper * rss_later_higher_co2)
    
    rss_current <- earlier_contribution * rss_early_current_co2 + later_contribution * rss_later_current_co2
    rss_current[topo_current <= 0] <- NA
    rss_current <- rss_current * (365*24*60*60*1e-6) # from [J s-1 m-2] to [MJ year-1 m-2]
    
    
    
    # 3) Total precipitation [m s-1]
    precip_early_lower_co2 <- raster(paste(prefix_climate,"run_",earlier_timestep/(-1e+6),"Ma/processed_",lower_CO2,"_ppm/timmean_",lower_CO2,"_ppm.nc", sep=""), var = "pr", level = 1)
    precip_early_lower_co2 <- rotate(precip_early_lower_co2)
    precip_early_higher_co2 <- raster(paste(prefix_climate,"run_",earlier_timestep/(-1e+6),"Ma/processed_",upper_CO2,"_ppm/timmean_",upper_CO2,"_ppm.nc", sep=""), var = "pr", level = 1)
    precip_early_higher_co2 <- rotate(precip_early_higher_co2)
    precip_early_current_co2 <- (CO2_contribution_lower * precip_early_lower_co2 + CO2_contribution_upper * precip_early_higher_co2)
    
    precip_later_lower_co2 <- raster(paste(prefix_climate,"run_",later_timestep/(-1e+6),"Ma/processed_",lower_CO2,"_ppm/timmean_",lower_CO2,"_ppm.nc", sep=""), var = "pr", level = 1)
    precip_later_lower_co2 <- rotate(precip_later_lower_co2)
    precip_later_higher_co2 <- raster(paste(prefix_climate,"run_",later_timestep/(-1e+6),"Ma/processed_",upper_CO2,"_ppm/timmean_",upper_CO2,"_ppm.nc", sep=""), var = "pr", level = 1)
    precip_later_higher_co2 <- rotate(precip_later_higher_co2)
    precip_later_current_co2 <- (CO2_contribution_lower * precip_later_lower_co2 + CO2_contribution_upper * precip_later_higher_co2)
    
    precip_current <- earlier_contribution * precip_early_current_co2 + later_contribution * precip_later_current_co2
    precip_current <- precip_current * 365*24*60*60 # conversion from [m s-1] to [m year-1]
    precip_current[topo_current <= 0] <- NA
    
    
    #Budyko aridity index -> AI = Rn/(Lamda * P), Rn->average daily radiation MJ m-2 day-1, P is average daily rainfall in m day-1, lambda is latent heat of evaporation in mJ m-3 and is a function of temperature lambda = (3.146 - 0.002361 Td) * 10^3 in MJ m-3
    # currently rss [W m-2] = [J s-1 m-2]: conversion factor *24*60*60*1e-6 [MJ day-1 m-2]
    budyko_aridity <- rss_current/(latent_heat*precip_current)
    # declare everyting above 7 as super arid 
    budyko_aridity[budyko_aridity >= 7] <- 7
    # Normalize and inverse 
    aridity <- 1-(budyko_aridity/7) # 0 --> full aridity, 1 --> no aridity
    
    
    
  
    #################
    #################
    # # Degassing approximation: Mills 
    
    Bforcing <- approx(c(-1000e+6, -150e+6, -100e+6, 0), y = c(0.75, 0.75, 1, 1), xout = t)$y
    if(degass_range == "mid"){
      degass_current_cc <- approx(degass_data_mills$time, degass_data_mills$D_mid, xout = t)$y * (k_ccdeg) * Bforcing
      degass_current_oc <- approx(degass_data_mills$time, degass_data_mills$D_mid, xout = t)$y * (k_ocdeg)
      degass_current <- degass_current_cc + degass_current_oc
    } else if(degass_range == "low"){
      degass_current_cc_mills <- approx(degass_data_mills$time, degass_data_mills$D_min, xout = t)$y * (k_ccdeg)
      degass_current_oc_mills <- approx(degass_data_mills$time, degass_data_mills$D_min, xout = t)$y * (k_ocdeg)
      # Take Marcilly as the low estimate 
      degass_current_cc_marcilly <- approx(degass_data_marcilly$time, degass_data_marcilly$fSR_zircon, xout = t)$y * (k_ccdeg)
      degass_current_oc_marcilly <- approx(degass_data_marcilly$time, degass_data_marcilly$fSR_zircon, xout = t)$y * (k_ocdeg)
      degass_current_cc <- min(c(degass_current_cc_mills, degass_current_cc_marcilly)) * Bforcing
      degass_current_oc <- min(c(degass_current_oc_mills, degass_current_oc_marcilly))
      degass_current <- degass_current_cc + degass_current_oc
    } else if(degass_range == "high"){
      degass_current_cc <- approx(degass_data_mills$time, degass_data_mills$D_max, xout = t)$y * (k_ccdeg) * Bforcing
      degass_current_oc <- approx(degass_data_mills$time, degass_data_mills$D_max, xout = t)$y * (k_ocdeg)
      degass_current <- degass_current_cc + degass_current_oc 
    }
      # degass_current <- sample(c(degass_current_low, degass_current_mid, degass_current_high), 1)
        
    #################
    #################
    # Land organic carbon burial and weathering enhancement 
    # Spatially explicit with evolution 
      
      # Transfer parameters to gen3sis
      parameter_list <- list()
      parameter_list[["thermal_adapt_par"]] <- thermal_adapt_par*DT # Scale parameter with the timestep
      parameter_list[["dispersal_par"]] <- dispersal_par # Scale dispersal parameter with the timestep, dispersal parameter is given in [m year-1]
      parameter_list[["timestep"]] <- length(timesteps)-step_index
      saveRDS(parameter_list, paste(experiment_directory,"/gen3sis/parameter_list.rds", sep = ""))
      setwd(paste(experiment_directory,"/gen3sis/", sep = ""))
      
      # Landscape depends on CO2 trajectory
      if(CO2_range == "mid"){
        landscape_directory <- "../../../../main_simulation/model_landscape/landscape_smooth_midCO2"
      } else if(CO2_range == "high"){
        landscape_directory <- "../../../../main_simulation/model_landscape/landscape_smooth_highCO2"
      } else if(CO2_range == "low"){
        landscape_directory <- "../../../../main_simulation/model_landscape/landscape_smooth_lowCO2"
      }

      
      
      biosphere_simulation <- run_simulation(config = "./biology_configuration_file.R",
                                             landscape = landscape_directory,
                                             verbose=T, #  progress printed
                                             output_directory="./model_output")
      #unlink("./biology/parameter_list.rds")
      setwd("../../")
      
      
      # Template
      template$area <- area(template) * 1e+6 # in m2
      template_df <- as.data.frame(template, xy=T)
      
      # Import landscape 
      current_gen3sis_time_step <- length(timesteps)-step_index  
      landscape <- readRDS(paste(experiment_directory,"/gen3sis/model_output/biology_configuration_file/landscapes/landscape_t_",current_gen3sis_time_step,".rds", sep=""))
      landscape_df <- as.data.frame(landscape$coordinates, xy=TRUE)
      landscape_df$temp <- landscape$environment[, "temp"]
      landscape_df$aridity <- landscape$environment[, "aridity"]
      landscape_df$topo <- landscape$environment[, "topo"]
      landscape_df$radiation <- landscape$environment[, "radiation"]
      output <- left_join(template_df, landscape_df, by=c("x", "y"))
      
      
      # Import biosphere state 
      biome_traits <- readRDS(paste(experiment_directory,"/gen3sis/model_output/biology_configuration_file/world_trait_map/map_t_",current_gen3sis_time_step,".rds", sep=""))
      biome_traits_df <- as.data.frame(biome_traits, xy=TRUE)
      output <- left_join(output, biome_traits_df, by=c("x", "y"))
      
       
      ######################################################################
      ################ Productivity calculations 
      #####
      NPP_init_adapt <- 0.01834285  # Energy conversion factor 
      
      # New vegetation module
      # Temperature step function 
      output$temp_limit <- approx(c(-30, -20, 35, 45), y = c(0, 1, 1, 0), xout = output$temp, method = "linear", n = 50,
                                  yleft = 0, yright = 0, rule = 1, f = 0, ties = mean, na.rm = TRUE)$y
      
      
      # Aridity response function 
      output$aridity_limit <- output$aridity
      
      # Adaptation scaling
      output$temp_adapt_limit <-  1 * exp(-0.25 * (output$Tdiff^2))
      output$temp_adapt_limit[which(output$temp_adapt_limit<0)] <- 0
      
      
      # Set location with no/extinct species adaptation to 0 
      output$temp_adapt_limit[is.na(output$Tdiff) & output$topo > 0] <- 0
  
      
      # Productivity 
      burial_rate <- 0.0007
      output$productivity_adapt <- NPP_init_adapt * output$radiation * output$temp_limit * output$aridity_limit * output$temp_adapt_limit * output$area * burial_rate 
      output$productivity_adapt_sum_global <- sum(output$productivity_adapt, na.rm=TRUE)
      output$weathering_limit_adapt <-  (output$radiation/5000) * output$temp_limit * output$aridity_limit * output$temp_adapt_limit 
      output$weathering_limit_adapt[output$weathering_limit_adapt > 1] = 1
      output$adaptation_limit <-  output$temp_adapt_limit 
      output$physiological_limit <-  (output$radiation/5000) * output$temp_limit * output$aridity_limit 
      output$physiological_limit[output$physiological_limit > 1] = 1
      # Assuming that after 5000 [MJ year-1 m-2], rss is not limiting anymore
      
      # Mean degree of adaptation 
      DoA <- mean(output$temp_adapt_limit, na.rm=T)
      
      # For calibration
      output$productivity_adapt_sum_global[1]
    
  
    
    
      ######################################################################
      ################ Silicate weathering calculation 
      #####
      
      # Calculate slope 
      output_raster <- rasterFromXYZ(output, crs="+proj=longlat +datum=WGS84 +no_defs")
      topo_slope <- output_raster$topo
      topo_slope[is.na(topo_slope)] <- 0 # assume 0 elevation for oceans, else no slope is calculated for coast
      slope <- terrain(topo_slope, opt="slope", unit="tangent", neighbors=8)
      slope[is.na(output$topo)] <- NA # Crop to continents
      output_raster$slope <- slope
      output <- as.data.frame(output_raster, xy=TRUE)
      
      
      # Runoff
      output_raster$runoff <- runoff_current

      # Add topography for completion
      output_raster$topo_complete <- topo_current
      output_raster$rss_complete <- rss_current
      output_raster$temp_complete <- temp_current

      output <- as.data.frame(output_raster, xy=TRUE)
      
      
      # Erosion without temperature dependency 
      # Scaled for mid pre-industrial CO2 estimate
      k_e <- 0.95e-2
      output$erosion_rate <- k_e * (output$runoff)^0.5 * output$slope # No temperature dependency here
      output$erosion_rate_area_scaled <- output$erosion_rate * output$area
      # Calibration of k_e --> 1.6e+10
      tot_erosion <- sum(output$erosion_rate_area_scaled, na.rm=TRUE)/1e+10
      # tot_erosion
      
      # Runoff dependency 
      k_w <- 1.5e-6 # Silicate weathering flow dependence
      f_Q <- (1-exp(-k_w*output$runoff))
      
      # Temperature depdendence
      E_a <- 10 
      R <- 8.314e-3
      T_0 <- 286
      f_T <- exp((E_a/(R*T_0))-(E_a/(R*(output$temp+273.15))))
      
      # Erosion dependency 
      sigma_factor <- 0.9 # Reaction time parameter
      Z <- 10 # Silicate weathering zone depth [m]
      f_E <- ((Z/output$erosion_rate)^sigma_factor)/(sigma_factor)
      f_E[is.infinite(f_E)] = 2e+7 
      
      f_kinetic = f_Q * f_T * f_E
      
      # Silicate weathering 
      X_m <- 0.1 # Silicate cation weight fraction
      K <- 6e-5 # Silicate weathering grain size dependence
      W_sil = output$erosion_rate * X_m * (1-exp(-K*f_kinetic))
      PREPLANT = PREPLANT_par
      bio_enhancement = ( 1 - pmin(output$weathering_limit_adapt , 1, na.rm=FALSE) ) * PREPLANT + (output$weathering_limit_adapt) # * ((CO2_current/278)^0.5)
      output$bio_enhancement <- bio_enhancement
      W_sil = W_sil  * bio_enhancement
      output$silw <- W_sil
      W_sil_area_scaled = W_sil * output$area 
      W_sil_tot = sum(W_sil_area_scaled, na.rm=T)
      # Calibration
      if(as.factor(PREPLANT) == "0.25"){
        output$silw_total_mol <- k_silw * (W_sil_tot/285251212)
      }
      if(as.factor(PREPLANT) != "0.25"){
        output$silw_total_mol <- k_silw * (W_sil_tot/274123604)
      }
      
      
      output_raster <- rasterFromXYZ(output, crs="+proj=longlat +datum=WGS84 +no_defs")
    
    #################
    #################
    # Marine organic carbon burial 
      
    NPPmarine_initial <- 3.157142 # Placeholder
      
    temp_sea <- temp_current
    temp_sea[topo_current > 0] <- NA
    Tlim_marine <- -3.27*(10^-8)*(temp_sea^7) + 3.4132*(10^-6)*(temp_sea^6) - 1.348*(10^-4)*(temp_sea^5) + 2.462*(10^-3)*(temp_sea^4) - 0.0205*(temp_sea^3) + 0.0617*(temp_sea^2) + 0.2749*temp_sea + 1.2956
    Tlim_marine[temp_sea > 33 | temp_sea <= 0] <- 0
      
    burial_rate <- 0.0007  
    NPPmarine <- NPPmarine_initial * Tlim_marine  * output_raster$area * burial_rate
    output_raster$NPPmarine <- NPPmarine
    mocb <- sum(values(NPPmarine), na.rm=TRUE)  
    # mocb
    
    #################
    #################
    # seafloor weathering
    # Global approximation - normalized to medium estimate of GAST pre industrial
    GAST_presentday <- 15.4775
    f_T_sfw = exp(0.0608*(GAST-GAST_presentday))  
    sfw = k_sfw * f_T_sfw * (degass_current/ (k_ocdeg + k_ccdeg))
    
    #################
    #################
    # carbw
    k_carbw_scale <- 0.00031
    carbw <- k_carbw_scale * output$runoff * output$bio_enhancement
    output$carbw <- carbw
    carbw_area_scaled <- carbw * output$area
    carbw_tot <- sum(carbw_area_scaled, na.rm=T)
    output$carbw_total_mol <- carbw_tot
    
    # caliration  - 8e+12
    # carbw_tot
    
    #################
    #################
    # oxidw
    # Global approximation - normalized to medium estimate of GAST pre industrial
    a = 0.5
    O2_present_day <- 21.62824
    O2_low <- predict(O2_spline_min, x = t)$y
    O2_high <- predict(O2_spline_max, x = t)$y
    O2_current <- mean(c(O2_low, O2_high))
    runoff_presentday <- 252.4714
    runoff_current_ave <- weighted.mean(values(runoff_current), values(area(runoff_current)), na.rm=T)
    oxidw <- k_oxidw * (tot_erosion/1.625522) * ((O2_current/O2_present_day)^a)    # (carbw_tot/k_carbw)
    
    #################
    #################
    # Plotting during code development
    # par(mfrow = c(2, 1))
    # # plot(topo_current_land, main="topo")
    # # plot(output_raster$bio_enhancement, main="bio_enhancement")
    # plot(output_raster$productivity_adapt,  main = "Productivity")
    # plot(output_raster$adaptation_limit, main = "State of Adaptation", zlim=c(0,1))
    
  
    #################
    #################
    # Record the fluxes 
    out[step_index, "time"] <- t
    out[step_index, "degass"] <- degass_current
    out[step_index, "silw"] <- output$silw_total_mol[1]
    out[step_index, "locb"] <- output$productivity_adapt_sum_global[1]
    out[step_index, "mocb"] <- mocb
    out[step_index, "sfw"] <- sfw
    out[step_index, "oxidw"] <- oxidw
    out[step_index, "dA"] <- degass_current + oxidw - output$silw_total_mol[1] - output$productivity_adapt_sum_global[1] - mocb - sfw
    out[step_index, "DoA"] <- DoA
    out[step_index, "GAST"] <- GAST
    out[step_index, "carbw"] <- carbw_tot
    out[step_index, "erosion"] <- tot_erosion
    out[step_index, "CO2"] <- CO2_current
    
    output_collection[[paste(step_index)]] <- output
    
    
    
    step_index <- step_index + 1
    print(t)
    print(GAST)
    print(DoA)
    
    
    if(t == 0){
      # Final calibration to get a steady state at present 
      out$silw <- out$silw *  (k_silw/tail(out$silw, n = 1))
      out$locb <- out$locb * (k_locb/tail(out$locb, n = 1))
      out$mocb <- out$mocb * (k_mocb/tail(out$mocb, n = 1))
      out$carbw <- out$carbw * (k_carbw/tail(out$carbw, n = 1))
      out$sfw <- out$sfw * (k_sfw/tail(out$sfw, n = 1))
      out$oxidw <- out$oxidw * (k_oxidw/tail(out$oxidw, n = 1))
      out$dA <- out$degass + out$oxidw - out$silw - out$locb - out$mocb - out$sfw
      
      
      # Add experiment info 
      out$therm_adapt <- thermal_adapt_par
      out$disp_par <- dispersal_par
      out$PREPLANT_par <- PREPLANT_par
      out$CO2_range <- CO2_range
      out$degass_range <- degass_range 
      out$ID <- ID
      
      
      # Write output 
      saveRDS(out, paste(experiment_directory,"/summary_out.rds", sep = ""))
      saveRDS(output_collection, paste(experiment_directory,"/detailed_out.rds", sep = ""))
      
    }
    
  } # End of timeloop
  
  } # Run only if output doesn't exist yet 
  return(out)
  
}) # End of lapply function
saveRDS(out, "compiled_output_lowresponse_A.rds")

# # 
