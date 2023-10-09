###############################################################################################################################################
############### Compile gen3sis landscape for C mass balance calculation model 

###############################################################################################################################################

# Working directory 
setwd(getwd())


library(parallel)
function_list <- list("low", "mid", "high")
cl <- makeCluster(3, outfile = "./landscape_creation_out.txt", overwrite = F)

parLapply(cl, function_list, fun = function(list_entry){
 
climate_model_timesteps <- seq(-400, 0, 5)*1e+6
climate_model_co2_levels <- c(50, seq(500, 3000, 500))
landscapes_list <- list() # , aridity = NULL temp = NULL,
CO2_range <- list_entry

# Libraries in worker
library(raster)
library(gen3sis)
library(rmatio)
library(dplyr)
library(geoscale)

# Prefixes - Euler vs. Local 
prefix_topo <- "./scotese_maps/"
prefix_climate <- "./climate_data/"

######### Produce CO2 curve 
CO2_data <- read.csv("./co2_loess_fit.csv", header=TRUE, sep = ";")
colnames(CO2_data)
colnames(CO2_data) <- c("Age", "pCO2_maxprob", "lw95", "lw68", "up68", "up95")
CO2_data$time <- CO2_data$Age * (-1e+6)
CO2_smooth_spline <- smooth.spline(CO2_data$time, CO2_data$pCO2_maxprob, spar = NULL)
CO2_smooth_spline_high <- smooth.spline(CO2_data$time, CO2_data$up68, spar = NULL)
CO2_smooth_spline_low <- smooth.spline(CO2_data$time, CO2_data$lw68, spar = NULL)

for (t in seq(0, -400e+6, -5e+5)){
  
    # CO2_values at t derived from smooth spline 
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
    #runoff_current <- runoff_current*(365*24*60*60*1000)


    ###### Aridity --> Three raster needed: 1) Temperature, 2) Net Shortwave radiation, 3) Precipitation
    # 1) Temperature
    # plot(Tair_raster_present)
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


    landscapes_list$temp <- c(landscapes_list$temp, temp_current)
    landscapes_list$topo <- c(landscapes_list$topo, topo_current_land)
    landscapes_list$aridity <- c(landscapes_list$aridity, aridity)
    landscapes_list$radiation <- c(landscapes_list$radiation, rss_current)

    # par(mfrow = c(2,2))
    # plot(topo_current_land, zlim=c(0, 5000))
    # plot(temp_current, zlim=c(-50, 40), main = paste(t))
    # plot(aridity, zlim=c(0, 1))
    # plot(rss_current, zlim=c(0, 9000))
    print(t)
    print(CO2_current)
    print(mean(values(temp_current), na.rm=T))
    
}
filename <- paste("landscapes_list_smooth_", CO2_range, "CO2.rds", sep = "")
saveRDS(landscapes_list, filename)
landscapes_list <- readRDS(filename)

cost_function_water_topo <- function(source, habitable_src, dest, habitable_dest) {
  if (!all(habitable_src, habitable_dest)) {
    return(4/1000) # Costs 4 times more to cross water cells
  } else {
    return(1/1000)
  }
}

directory_name <- paste("./model_landscape/landscape_smooth_",CO2_range,"CO2", sep = "")
create_input_landscape(landscapes = landscapes_list,
                       timesteps = as.character(seq(0,length(seq(-400e+6, 0, 5e+5))-1, 1)),
                       cost_function = cost_function_water_topo,
                       directions = 8,
                       output_directory = directory_name,
                       overwrite = T,
                       crs = "+proj=longlat +datum=WGS84 +no_defs",
                       calculate_full_distance_matrices = T,
                       verbose = T)

})


