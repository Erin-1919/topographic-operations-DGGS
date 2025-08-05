# DGGS Topographic Analysis - R Scripts

This directory contains professional R scripts for generating DGGS (Discrete Global Grid System) components for topographic analysis. The scripts have been improved with modern R practices, comprehensive documentation, error handling, and logging.

## Overview

The R scripts provide functionality for:
- Generating centroids within areas of interest
- Creating hexagonal grids with DGGS coordinates
- Converting between geographic and DGGS coordinate systems
- Supporting the Python analysis pipeline

## Scripts

### 1. Centroids_q2di.R
**Purpose**: Generate centroids within areas of interest for DGGS analysis

**Key Features**:
- Creates regular grid of points within study area bounds
- Converts geographic coordinates to DGGS Q2DI coordinates
- Removes duplicate centroids
- Outputs CSV file with centroid coordinates and DGGS indices

**Usage**:
```bash
Rscript Centroids_q2di.R <resolution> <study_area>
```

**Example**:
```bash
Rscript Centroids_q2di.R 20 Calgary
```

**Output**:
- `Result/Area_{study_area}_centroids_{resolution}.csv`
- Contains columns: `i`, `j`, `lon_c`, `lat_c`

### 2. Grids_q2di.R
**Purpose**: Generate hexagonal grids within areas of interest

**Key Features**:
- Creates hexagonal grid polygons within study area bounds
- Converts grid centroids to DGGS Q2DI coordinates
- Outputs shapefile with grid polygons and DGGS indices

**Usage**:
```bash
Rscript Grids_q2di.R <resolution> <study_area>
```

**Example**:
```bash
Rscript Grids_q2di.R 20 Calgary
```

**Output**:
- `Result/Area_{study_area}_hex_{resolution}.shp`
- Shapefile with grid polygons and `i`, `j` columns

## Requirements

### R Dependencies
```r
install.packages(c("dggridR", "rgdal", "rgeos", "dplyr", "geosphere"))
```

### Required Data Structure
```
Data/
├── Area_Calgary.shp
├── Area_Canmore.shp
└── Area_BuffaloLake.shp

Result/
└── (output files will be created here)
```

## Input Parameters

### Resolution Levels
- Valid range: 19-24
- Recommended: 20-22 for analysis
- Higher resolution = smaller cells, more detail

### Study Areas
- Calgary
- Canmore  
- BuffaloLake

## DGGS Configuration

The scripts use the following DGGS configuration:
- **Projection**: ISEA (Icosahedral Snyder Equal Area)
- **Aperture**: 3 (hexagonal)
- **Topology**: HEXAGON
- **Pole**: Latitude 31°, Longitude -56°
- **Azimuth**: 0°

## Code Improvements

### Professional Standards
- **Documentation**: Comprehensive roxygen2-style documentation
- **Error Handling**: Robust exception handling with meaningful messages
- **Logging**: Structured logging with timestamps and levels
- **Input Validation**: Parameter validation with helpful error messages
- **Code Organization**: Clear separation of concerns and modular design

### Performance Optimizations
- **Efficient Algorithms**: Optimized coordinate transformations
- **Memory Management**: Proper cleanup of spatial objects
- **I/O Optimization**: Minimizes file operations
- **Resource Management**: Proper handling of spatial data

### Maintainability
- **Consistent Naming**: R style guide compliant variable and function names
- **Modular Design**: Reusable functions and clear interfaces
- **Configuration**: Centralized constants and lookup tables
- **Documentation**: Inline comments and comprehensive README

## Error Handling

The scripts include comprehensive error handling:
- File existence checks
- Parameter validation
- Exception catching with fallback values
- Graceful degradation for edge cases
- Informative error messages

## Logging

All scripts use structured logging:
- INFO level for progress updates
- WARNING level for non-critical issues
- ERROR level for failures
- Timing information for performance monitoring

## Coordinate Systems

### Geographic Coordinates
- **Input**: WGS84 geographic coordinates (longitude, latitude)
- **Output**: WGS84 geographic coordinates for centroids

### DGGS Coordinates
- **System**: Q2DI (Quad-2D-Index)
- **Components**: `quad`, `i`, `j` indices
- **Conversion**: Using dggridR package functions

## Data Flow

1. **Input**: Study area shapefile (WGS84)
2. **Processing**: 
   - Extract bounding box
   - Generate regular grid
   - Convert to DGGS coordinates
   - Remove duplicates
3. **Output**: CSV/Shapefile with DGGS indices

## Performance Considerations

- **Grid Generation**: Efficient hexagonal grid creation
- **Coordinate Conversion**: Optimized DGGS transformations
- **Memory Usage**: Minimal memory footprint for large areas
- **I/O Operations**: Efficient file reading and writing

## Authors

- Mingke Li
- Heather McGrath  
- Emmanuel Stefanakis

## Citation

If you use this code in your research, please cite:
```
Li, M., McGrath, H., & Stefanakis, E. (2022). Multi-resolution topographic analysis in hexagonal Discrete Global Grid Systems. International Journal of Applied Earth Observation and Geoinformation, 107, 102985.
```
