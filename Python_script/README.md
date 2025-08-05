# DGGS Topographic Analysis - Python Scripts

This directory contains professional Python scripts for performing topographic analysis in hexagonal Discrete Global Grid Systems (DGGS). The scripts have been improved with modern Python practices, comprehensive documentation, error handling, and logging.

## Overview

The scripts implement various topographic analysis algorithms for hexagonal DGGS, including:
- Slope calculation using multiple algorithms (MAG, MDG, MDN, FDA, BFP)
- Curvature analysis
- Hillshade generation
- Topographic indices (TRI, TPI)
- Zonal statistics
- Focal statistics
- Elevation quantization

## Scripts

### 1. DgBaseFunc.py
**Purpose**: Core utility functions for DGGS operations

**Key Features**:
- DGGS resolution lookup table
- Neighbor coordinate calculations
- Edge cell detection
- Mathematical operations (cross product, derivatives)
- Aspect calculations

**Main Functions**:
- `look_up_table()`: DGGS resolution parameters
- `neighbor_coords()`: Calculate neighbor coordinates
- `neighbor_navig()`: Extract elevation values from neighbors
- `first_derivative()` / `second_derivative()`: Gradient calculations
- `aspect_restricted()` / `aspect_unrestricted()`: Aspect calculations

### 2. DgTopogAnaly.py
**Purpose**: Slope and topographic analysis using multiple algorithms

**Algorithms Implemented**:
- **MAG** (Maximum Absolute Gradient): Absolute maximum differences between center and neighbors
- **MDG** (Maximum Downward Gradient): Maximum downward gradient with pit-filling
- **MDN** (Multiple Downhill Neighbors): Average of all downslope neighbors
- **FDA** (Finite-Difference Algorithm): Project gradient to orthogonal axes
- **BFP** (Best Fit Plane): Surface fitting using least squares

**Usage**:
```bash
python DgTopogAnaly.py <resolution> <area> <method>
```

**Example**:
```bash
python DgTopogAnaly.py 20 Calgary FDA
```

### 3. DgTopogIndex.py
**Purpose**: Calculate topographic indices

**Indices**:
- **TRI** (Terrain Roughness Index): Square root of sum of squared differences
- **TPI** (Topographic Position Index): Difference from mean of surrounding cells

**Usage**:
```bash
python DgTopogIndex.py <resolution> <area> <method>
```

### 4. DgZonalStats.py
**Purpose**: Zonal statistics and raster operations

**Features**:
- Raster projection to NAD83 CSRS
- Vegetation classification
- Statistical aggregation by zones
- Volume calculations

**Usage**:
```bash
python DgZonalStats.py <resolution> <area>
```

### 5. DgFocalStats.py
**Purpose**: Focal statistics for neighborhood analysis

**Statistics**:
- Mean, max, min, median, standard deviation, range
- Configurable ring expansion
- Void handling

**Usage**:
```bash
python DgFocalStats.py <resolution> <ring_number> <area>
```

### 6. DgQuantization.py
**Purpose**: Elevation quantization and DEM resampling

**Features**:
- Bilinear interpolation
- No-data handling
- Precision control
- Parallel processing

**Usage**:
```bash
python DgQuantization.py <resolution> <area>
```

## Requirements

### Python Dependencies
```bash
pip install pandas numpy scipy rasterio gdal multiprocess
```

### Required Data Structure
```
Data/
├── Calgary.tif
├── Canmore.tif
├── BuffaloLake.tif
└── aci_2020_ab.tif

Result/
├── Area_Calgary_centroids_20.csv
├── Calgary_elev_20.csv
└── ...
```

## Input Parameters

### Resolution Levels
- Valid range: 16-29
- Recommended: 20-24 for analysis
- Higher resolution = smaller cells, more detail

### Study Areas
- Calgary
- Canmore  
- BuffaloLake

### Methods (for slope analysis)
- MAG: Maximum Absolute Gradient
- MDG: Maximum Downward Gradient
- MDN: Multiple Downhill Neighbors
- FDA: Finite-Difference Algorithm
- BFP: Best Fit Plane

## Output Files

### Slope Analysis
- `{area}_elev_{method}_{resolution}.csv`: Slope and aspect results
- `{area}_curvature_{resolution}.csv`: Curvature values
- `{area}_hillshade_{resolution}.csv`: Hillshade values

### Topographic Indices
- `{area}_TopoIndex_{method}_{resolution}.csv`: TRI and TPI values

### Zonal Statistics
- `{area}_elev_zonal_{resolution}.csv`: Statistics by vegetation zones

### Focal Statistics
- `{area}_elev_focal_{resolution}_ring{rings}.csv`: Neighborhood statistics

### Elevation Data
- `{area}_elev_{resolution}.csv`: Quantized elevation values

## Code Improvements

### Professional Standards
- **Type Hints**: All functions include proper type annotations
- **Docstrings**: Comprehensive documentation with Args/Returns sections
- **Error Handling**: Robust exception handling with meaningful messages
- **Logging**: Structured logging with different levels
- **Input Validation**: Parameter validation with helpful error messages
- **Code Organization**: Clear separation of concerns and modular design

### Performance Optimizations
- **Parallel Processing**: Multiprocessing for large datasets
- **Memory Management**: Context managers for file handling
- **Efficient Algorithms**: Optimized mathematical operations
- **Resource Cleanup**: Proper pool management and file closing

### Maintainability
- **Consistent Naming**: PEP 8 compliant variable and function names
- **Modular Design**: Reusable functions and clear interfaces
- **Configuration**: Environment-based settings (SLURM_CPUS_PER_TASK)
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

## Performance Considerations

- **Parallel Processing**: Automatically detects available cores
- **Memory Efficient**: Processes data in chunks
- **I/O Optimization**: Minimizes file operations
- **Algorithm Efficiency**: Optimized mathematical calculations

## Authors

- Mingke Li
- Heather McGrath  
- Emmanuel Stefanakis

## Citation

If you use this code in your research, please cite:
```
Li, M., McGrath, H., & Stefanakis, E. (2022). Multi-resolution topographic analysis in hexagonal Discrete Global Grid Systems. International Journal of Applied Earth Observation and Geoinformation, 107, 102985.
``` 