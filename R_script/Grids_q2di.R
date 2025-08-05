################################################################################
# DGGS Hexagonal Grids Generation Script
# 
# This script generates hexagonal grids within areas of interest for 
# Discrete Global Grid Systems (DGGS) using the dggridR package.
#
# Author: Mingke Li
# 
# Usage: Rscript Grids_q2di.R <resolution> <study_area>
# Example: Rscript Grids_q2di.R 20 Calgary
################################################################################

# Load required libraries
suppressPackageStartupMessages({
  library(dggridR)
  library(rgdal)
  library(rgeos)
  library(dplyr)
  library(geosphere)
})

# Configure logging
log_message <- function(message, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] %s: %s\n", timestamp, level, message))
}

# Error handling function
handle_error <- function(e, context = "") {
  log_message(sprintf("Error in %s: %s", context, e$message), "ERROR")
  quit(status = 1)
}

# Input validation function
validate_inputs <- function(resolution, study_area) {
  # Validate resolution
  valid_resolutions <- c(19, 20, 21, 22, 23, 24)
  if (!resolution %in% valid_resolutions) {
    stop(sprintf("Invalid resolution: %d. Must be one of: %s", 
                 resolution, paste(valid_resolutions, collapse = ", ")))
  }
  
  # Validate study area
  valid_areas <- c("Calgary", "Canmore", "BuffaloLake")
  if (!study_area %in% valid_areas) {
    stop(sprintf("Invalid study area: %s. Must be one of: %s", 
                 study_area, paste(valid_areas, collapse = ", ")))
  }
  
  # Check if input file exists
  area_file <- sprintf("Data/Area_%s.shp", study_area)
  if (!file.exists(area_file)) {
    stop(sprintf("Study area file not found: %s", area_file))
  }
  
  log_message(sprintf("Input validation passed - Resolution: %d, Area: %s", 
                      resolution, study_area))
}

# DGGS configuration constants
DGGS_CONFIG <- list(
  vlat = 31,
  vlon = -56,
  vazimuth = 0,
  projection = "ISEA",
  aperture = 3,
  topology = "HEXAGON"
)

# Resolution lookup table
RESOLUTION_LOOKUP <- data.frame(
  res_list = c(19, 20, 21, 22, 23, 24),
  cell_size_list = c(0.0009, 0.0008, 0.0003, 0.0003, 0.0001, 0.0001)
)

#' Generate hexagonal grids within area of interest
#' 
#' This function generates hexagonal grids for DGGS cells within a specified
#' study area. It creates a regular grid of hexagonal polygons and converts
#' them to DGGS coordinates using the dggridR package.
#' 
#' @param area Character string specifying the study area name
#' @param resolution Integer specifying the DGGS resolution level
#' @return NULL (writes output to shapefile)
#' @examples
#' generate_grids("Calgary", 20)
generate_grids <- function(area, resolution) {
  tryCatch({
    log_message(sprintf("Starting grid generation for area: %s, resolution: %d", 
                        area, resolution))
    
    # Read study area shapefile
    log_message("Reading study area shapefile...")
    study_area <- readOGR(dsn = "Data", layer = area)
    study_area_bbox <- bbox(study_area[1,])
    
    # Extract bounding box coordinates
    minx <- study_area_bbox[1, 1]
    miny <- study_area_bbox[2, 1]
    maxx <- study_area_bbox[1, 2]
    maxy <- study_area_bbox[2, 2]
    
    log_message(sprintf("Study area bounds: lon [%.4f, %.4f], lat [%.4f, %.4f]", 
                        minx, maxx, miny, maxy))
    
    # Get cell size from lookup table
    dggs_cellsize <- RESOLUTION_LOOKUP$cell_size_list[RESOLUTION_LOOKUP$res_list == resolution]
    log_message(sprintf("Cell size: %.6f", dggs_cellsize))
    
    # Construct DGGS
    log_message("Constructing DGGS...")
    dgg <- dgconstruct(
      projection = DGGS_CONFIG$projection,
      aperture = DGGS_CONFIG$aperture,
      topology = DGGS_CONFIG$topology,
      res = resolution,
      azimuth_deg = DGGS_CONFIG$vazimuth,
      pole_lat_deg = DGGS_CONFIG$vlat,
      pole_lon_deg = DGGS_CONFIG$vlon
    )
    
    # Generate hexagonal grid
    log_message("Generating hexagonal grid...")
    grid <- dgrectgrid(
      dgg, 
      minlat = miny, 
      minlon = minx, 
      maxlat = maxy, 
      maxlon = maxx, 
      cellsize = dggs_cellsize, 
      frame = FALSE
    )
    
    log_message(sprintf("Generated %d grid cells", length(grid)))
    
    # Calculate centroids for coordinate conversion
    log_message("Calculating grid centroids...")
    grid_centroids <- centroid(grid)
    
    # Convert to SpatialPolygonsDataFrame
    grid <- as(grid, "SpatialPolygonsDataFrame")
    
    # Convert centroids to DGGS coordinates
    log_message("Converting centroids to DGGS coordinates...")
    dggs_coords <- dgGEO_to_Q2DI(dgg, grid_centroids[, 1], grid_centroids[, 2])
    grid@data$i <- dggs_coords$i
    grid@data$j <- dggs_coords$j
    
    # Create output directory if it doesn't exist
    if (!dir.exists("Result")) {
      dir.create("Result", recursive = TRUE)
      log_message("Created Result directory")
    }
    
    # Save to shapefile
    output_file <- sprintf("Result/Area_%s_hex_%d.shp", area, resolution)
    writeOGR(grid, output_file, layer = "cell", driver = 'ESRI Shapefile')
    
    log_message(sprintf("Grid saved to: %s", output_file))
    log_message(sprintf("Grid generation completed successfully. Total cells: %d", 
                        nrow(grid@data)))
    
  }, error = function(e) {
    handle_error(e, "generate_grids")
  })
}

#' Main execution function
#' 
#' Parses command line arguments and executes the grid generation process
main <- function() {
  tryCatch({
    # Parse command line arguments
    args <- commandArgs(trailingOnly = TRUE)
    
    if (length(args) != 2) {
      stop("Usage: Rscript Grids_q2di.R <resolution> <study_area>")
    }
    
    dggs_res <- as.numeric(args[1])
    study_area <- args[2]
    
    # Validate inputs
    validate_inputs(dggs_res, study_area)
    
    # Execute grid generation
    generate_grids(sprintf("Area_%s", study_area), dggs_res)
    
    log_message("Script execution completed successfully")
    
  }, error = function(e) {
    log_message(sprintf("Script failed: %s", e$message), "ERROR")
    quit(status = 1)
  })
}

# Execute main function if script is run directly
if (!interactive()) {
  main()
}

