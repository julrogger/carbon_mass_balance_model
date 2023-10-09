
######################################
###            METADATA            ###
######################################
# gen3sis configuration
#
# Version: 1.0
#
# Author: Julian Rogger
#
# Date: 28.03.2022
#
# Landscape: global_Phanerozoic_v1
#
# Publications:
# Description: 
#
######################################


######################################
###         General settings       ###
######################################

# set the random seed for the simulation.
random_seed = 666

# Transfer parameters from frontend needed in the configuration file 
# paste("reading in the parameters")
parameter_list <- readRDS("parameter_list.rds")
thermal_adapt_par <- parameter_list[["thermal_adapt_par"]]
dispersal_par <- parameter_list[["dispersal_par"]]
timestep <- parameter_list[["timestep"]]


# set the starting time step or leave NA to use the earliest/highest time-step.
start_time = timestep

# set the end time step or leave as NA to use the latest/lowest time-step (0).
end_time = start_time

# maximum total number of species in the simulation before it is aborted.
max_number_of_species = 100000

# maximum number of species within one cell before the simulation is aborted.
max_number_of_coexisting_species = 100000

# a list of traits to include with each species
# a "dispersal" trait is implicitly added in any case
trait_names = c("Topt", "Tdiff")

# ranges to scale the input environments with:
# not listed variable:         no scaling takes place
# listed, set to NA:           the environmental variable will be scaled from [min, max] to [0, 1]
# listed with a given range r: the environmental variable will be scaled from [r1, r2] to [0, 1]
environmental_ranges = NA



######################################
###            Observer            ###
######################################

# a place to inspect the internal state of the simulation and collect additional information if desired.
end_of_timestep_observer = function(data, vars, config){
  # browser()
  save_landscape()
  save_species()
  landscape_coordinates <- as.data.frame(data$landscape$coordinates)
  landscape_coordinates$Topt <- NA
  landscape_coordinates$Tdiff <- NA
  landscape_coordinates$species_id <- NA
  for(i in 1:length(data$all_species)){
    traits <- data$all_species[[i]]$traits
    landscape_coordinates$Topt[which(rownames(landscape_coordinates) %in% rownames(traits))] <- traits[, "Topt"]
    landscape_coordinates$Tdiff[which(rownames(landscape_coordinates) %in% rownames(traits))] <- traits[, "Tdiff"]
    landscape_coordinates$species_id[which(rownames(landscape_coordinates) %in% rownames(traits))] <- i
  }
  traits_world <- raster::rasterFromXYZ(landscape_coordinates, res=c(3.75, 3.75) , crs="+proj=longlat +datum=WGS84 +no_defs")
  saveRDS(traits_world, file=paste(config$directories$output,"/world_trait_map/map_t_",vars$ti,".rds", sep=""))
}


######################################
###         Initialization         ###
######################################

# the initial abundance of a newly colonized cell, both during setup and later when 
# colonizing a cell during the dispersal.
initial_abundance = 1

# place species in the landscape:
create_ancestor_species <- function(landscape, config) {
  # browser()
  if(timestep == 780){
    
    landcells <- length(landscape$environment[, 1]) # as many ancestors as there are land cells - each optimally adapted
    all_species <- list()
    for(i in 1:landcells){
      new_species <- create_species(rownames(landscape$environment)[i], config)
      new_species$traits[ , "Topt"] <- landscape$environment[rownames(landscape$environment)[i], "temp"] + rnorm(1, mean = 0, sd = 1.0)
      new_species$traits[ , "Tdiff"] <- abs(new_species$traits[, "Topt"] - landscape$environment[rownames(landscape$environment)[i], "temp"])
      
      all_species <- append(all_species, list(new_species))
    }
    return(all_species)
    
  } else {

    # Import species from last time step
    filename <- paste('/species/species_t_', timestep + 1, '.rds', sep="")
    species <- readRDS(paste(config$directories$output, filename, sep=""))
    
    # habitable cells
    # browser()
    habitable <- rownames(landscape[["coordinates"]])
    species <- lapply(species, FUN = function(input){
      limited_cells <- names(input[["abundance"]])
      limited_cells <- limited_cells[which(limited_cells %in% habitable)]
      
      input[["abundance"]] <- input[["abundance"]][limited_cells]
      input[["traits"]] <- input[["traits"]][limited_cells, ,drop = FALSE]
      input[["divergence"]] <- limit_divergence_to_cells(input[["divergence"]], limited_cells)
      
      return(input)
    })
    
    # browser()
    index_list <- lapply(species, FUN = function(input){
      index <- NA
      if(length(input$abundance) >= 1){index <- 1}
      if(length(input$abundance) == 0){index <- 0}
      return(index)
    })
    index_vec <- do.call(rbind, index_list)
    if(any(index_vec == 0)){
    index <- which(index_vec == 0)
    species <- species[-c(index)]
    }
    


    return(species) 
    
    
  }
}



######################################
###             Dispersal          ###
######################################

# the maximum range to consider when calculating the distances from local distance inputs.
max_dispersal <- Inf

# returns n dispersal values.
get_dispersal_values <- function(n, species, landscape, config) {
  
  # Dispersal values only depends on aridity 
  values <- rep(dispersal_par, n)
  occurrence <- names(species$abundance)
  aridity_occurrence <- landscape$environment[occurrence, "aridity"]
  ariditylim_occurrence <- aridity_occurrence
  ariditylim_occurrence[ariditylim_occurrence < 0.25] <- 0.25
  return(values * ariditylim_occurrence)
  
  }


######################################
###          Speciation            ###
######################################

# threshold for genetic distance after which a speciation event takes place.
divergence_threshold = 1
# factor by which the divergence is increased between geographically isolated population.
# can also be a matrix between the different population clusters.
get_divergence_factor <- function(species, cluster_indices, landscape, config) {
  
  return(1) # as soon as "species"/biomes et separated they evolve separately the next timestep
}


######################################
###      Trait Evolution           ###
######################################
adaptation_increment <- thermal_adapt_par # flags - imported from frontend

# mutate the traits of a species and return the new traits matrix.
apply_evolution <- function(species, cluster_indices, landscape, config){
  
  traits <- species[["traits"]]
  cells <- rownames(traits)
  
  # evolve Topt
  temperature_difference <- landscape$environment[cells, "temp"] - traits[cells, "Topt"]
  temperature_adaptation_change <- temperature_difference
  temperature_adaptation_change[temperature_adaptation_change >= adaptation_increment] <-  adaptation_increment
  temperature_adaptation_change[temperature_adaptation_change <= (-adaptation_increment)] <- (-adaptation_increment)
  traits[cells, "Topt"] <- traits[cells, "Topt"] + temperature_adaptation_change 
  traits[cells, "Tdiff"] <- abs(traits[cells, "Topt"] - landscape$environment[cells,"temp"])
  
  return(traits)
}


######################################
###             Ecology            ###
######################################

# called for every cell with all occurring species, this function calculates abundances and/or 
# who survives for each sites.
# returns a vector of abundances.
# set the abundance to 0 for every species supposed to die.
# 


apply_ecology <- function(abundance, traits, landscape, config) {
  K <- 1
  # Competition - best adapted biome survives 
  Tair_diff <- abs( traits[, "Topt"] - landscape[, "temp"])
  
  Tair_limit <- 1 * exp(-0.02 * (Tair_diff^2))
  performance <- Tair_limit
  
  abundance[names(abundance[performance != max(performance)])] = 0
  abundance[names(abundance[performance == max(performance)])] = 1
  
  survivor_name <- sample(names(abundance[abundance == 1]), 1)  
  abundance[survivor_name] = 1
  abundance[names(abundance) != survivor_name] = 0
  
  # Hard limits on temperature difference and aridity 
  # This will kill biomes, in the next step a new one will take over or space is left empty
  # abundance[which(traits[, "Tdiff"] >= 15)] <- 0 
  
  
  return(abundance)
}
