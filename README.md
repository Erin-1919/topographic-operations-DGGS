# Multi-resolution topographic analysis in hexagonal Discrete Global Grid Systems -- Source Code

This repository contains the source code for research on multi-resolution topographic analysis in hexagonal Discrete Global Grid Systems (DGGS). The code has been professionally improved with modern programming practices, comprehensive documentation, error handling, and logging.

## Overview

This work develops analytical functions based on terrain data in a pure environment of ISEA3H DGGS. The developed operations include descriptive statistics, topographic analysis, and topographic indices, classified into three categories: local, zonal, and focal operations. The open-sourced library [*dggridR*](https://github.com/r-barnes/dggridR) is used to complete conversion between geographic locations and ISEA3H DGGS cell indices.

The experiment was carried out using a hybrid of Python 3.7.7 and R 3.6.2 environments. The code has been significantly improved with professional coding standards and is available in the folders [*R_script*](https://github.com/Erin-1919/Topographic-operations-DGGS/tree/main/R_script) and [*Python_script*](https://github.com/Erin-1919/Topographic-operations-DGGS/tree/main/Python_script).

## Manuscript Information

### Title of Manuscript
Multi-resolution topographic analysis in hexagonal Discrete Global Grid Systems

### Keywords
Discrete Global Grid Systems; topographical analysis; multi-resolution; map algebra

### DOI
10.1016/j.jag.2022.102985

### Authors
Mingke Li, Heather McGrath, and Emmanuel Stefanakis

### Corresponding Author
[Mingke Li](https://erin-1919.github.io/) (mingke.li@ucalgary.ca)

[ORCID](https://orcid.org/0000-0001-6310-4964)

### Abstract
Discrete Global Grid Systems (DGGS) have been increasingly adopted as a standard framework for multi-source geospatial data. Previous research largely studied the mathematical foundation of discrete global grids, developed open-source libraries, and explored their application as data integration platforms. This study investigated the multi-resolution terrain analysis in a pure hexagonal DGGS environment, including descriptive statistics, topographic parameters, and topographic indices. Experiments across multiple grid resolutions were carried out in three study areas with different terrain roughness in Alberta, Canada. Five algorithms were proposed to calculate both the slope gradient and terrain aspect. A cell-based pair-wise comparison showed a strong positive correlation between the gradient values as calculated from five algorithms. The grid resolutions as well as the terrain roughness had a clear effect on the computed slope gradient and topographic indices. This research aims to enhance the analytical functionality of hexagonal DGGS to better support decision-making in real world problems.

### Code Repository
https://github.com/Erin-1919/Topographic-operations-DGGS

## Code Structure

### R Scripts (`R_script/`)
- **Centroids_q2di.R**: Generate centroids within areas of interest for DGGS analysis
- **Grids_q2di.R**: Generate hexagonal grids with DGGS coordinates

### Python Scripts (`Python_script/`)
- **DgBaseFunc.py**: Core utility functions for DGGS operations
- **DgTopogAnaly.py**: Slope analysis using multiple algorithms (MAG, MDG, MDN, FDA, BFP)
- **DgTopogIndex.py**: Calculate topographic indices (TRI, TPI)
- **DgZonalStats.py**: Zonal statistics and raster operations
- **DgFocalStats.py**: Focal statistics for neighborhood analysis
- **DgQuantization.py**: Elevation quantization and DEM resampling

## Professional Improvements

### Code Quality Standards
- **Type Hints**: All Python functions include proper type annotations
- **Documentation**: Comprehensive docstrings with Args/Returns sections
- **Error Handling**: Robust exception handling with meaningful messages
- **Logging**: Structured logging with different levels (INFO, WARNING, ERROR)
- **Input Validation**: Parameter validation with helpful error messages
- **Code Organization**: Clear separation of concerns and modular design

### Performance Optimizations
- **Parallel Processing**: Multiprocessing for large datasets
- **Memory Management**: Context managers for file handling
- **Efficient Algorithms**: Optimized mathematical operations
- **Resource Cleanup**: Proper pool management and file closing

### Maintainability
- **Consistent Naming**: PEP 8 (Python) and R style guide compliant names
- **Modular Design**: Reusable functions and clear interfaces
- **Configuration**: Environment-based settings and centralized constants
- **Documentation**: Inline comments and comprehensive README files

## Quick Start

### Prerequisites
```bash
# Python dependencies
pip install pandas numpy scipy rasterio gdal multiprocess

# R dependencies
R -e "install.packages(c('dggridR', 'rgdal', 'rgeos', 'dplyr', 'geosphere'))"
```

### Basic Usage

1. **Generate DGGS Components (R)**:
```bash
# Generate centroids
Rscript R_script/Centroids_q2di.R 20 Calgary

# Generate hexagonal grids
Rscript R_script/Grids_q2di.R 20 Calgary
```

2. **Perform Topographic Analysis (Python)**:
```bash
# Elevation quantization
python Python_script/DgQuantization.py 20 Calgary

# Slope analysis
python Python_script/DgTopogAnaly.py 20 Calgary FDA

# Topographic indices
python Python_script/DgTopogIndex.py 20 Calgary FDA

# Zonal statistics
python Python_script/DgZonalStats.py 20 Calgary

# Focal statistics
python Python_script/DgFocalStats.py 20 1 Calgary
```

## Input Parameters

### Resolution Levels
- **Valid range**: 16-29 (Python), 19-24 (R)
- **Recommended**: 20-24 for analysis
- **Higher resolution**: Smaller cells, more detail

### Study Areas
- Calgary
- Canmore  
- BuffaloLake

### Slope Analysis Methods
- **MAG**: Maximum Absolute Gradient
- **MDG**: Maximum Downward Gradient
- **MDN**: Multiple Downhill Neighbors
- **FDA**: Finite-Difference Algorithm
- **BFP**: Best Fit Plane

## Output Files

### R Scripts Output
- `Result/Area_{study_area}_centroids_{resolution}.csv`: Centroid coordinates
- `Result/Area_{study_area}_hex_{resolution}.shp`: Hexagonal grid shapefile

### Python Scripts Output
- `Result/{area}_elev_{resolution}.csv`: Quantized elevation data
- `Result/{area}_elev_{method}_{resolution}.csv`: Slope and aspect results
- `Result/{area}_curvature_{resolution}.csv`: Curvature values
- `Result/{area}_hillshade_{resolution}.csv`: Hillshade values
- `Result/{area}_TopoIndex_{method}_{resolution}.csv`: Topographic indices
- `Result/{area}_elev_zonal_{resolution}.csv`: Zonal statistics
- `Result/{area}_elev_focal_{resolution}_ring{rings}.csv`: Focal statistics

## Libraries Used

### Python
- numpy 1.19.4
- scipy 1.5.3
- rasterio 1.2.1
- gdal 3.1.4
- pandas 1.1.4
- multiprocess 0.70.12.2

### R
- dggridR 2.0.4
- rgdal 1.5.16
- rgeos 0.5.5
- dplyr 1.0.2
- geosphere 1.5.10

## Data Availability

The original Canadian Digital Elevation Model (CDEM) data can be downloaded via the Geospatial-Data Extraction tool in [Canada's Open Government Portal](https://maps.canada.ca/czs/index-en.html), or they can be obtained through the STAC API. The 2020 Annual Crop Inventory produced by the Agriculture and Agri-Food Canada, used as an example in zonal statistics, is accessible [here](https://open.canada.ca/data/en/dataset/ba2645d5-4458-414d-b196-6303ac06c1c9).

## Error Handling & Logging

All scripts include comprehensive error handling:
- File existence checks
- Parameter validation
- Exception catching with fallback values
- Graceful degradation for edge cases
- Informative error messages

Structured logging provides:
- Progress updates during processing
- Performance monitoring
- Clear success/failure reporting
- Timing information for optimization


## Citation

If you use this code in your research, please cite:
```
Li, M., McGrath, H., & Stefanakis, E. (2022). Multi-resolution topographic analysis in hexagonal Discrete Global Grid Systems. International Journal of Applied Earth Observation and Geoinformation, 107, 102985.
```
