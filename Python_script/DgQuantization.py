"""
Discrete Global Grid System (DGGS) Quantization

This module provides elevation quantization functions for hexagonal DGGS,
including DEM resampling and interpolation.

Author: Mingke Li
"""

import pandas as pd
import numpy as np
import rasterio
import time
import sys
import os
import warnings
import multiprocess as mp
import logging
from typing import List, Tuple, Optional, Union, Any
from pathlib import Path
from scipy import interpolate

import DgBaseFunc

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

warnings.simplefilter('error', RuntimeWarning)


def find_neighbor(x: float, y: float, dem_tif: rasterio.DatasetReader) -> Tuple[List[float], List[float], List[float]]:
    """
    Find neighbors for interpolation.
    
    Determine the DEM source, find out the 4 neighbors geographic coords,
    extract the elevations at these 4 coords, convert 4 coords to array index 
    then back to 4 coords of grid mesh centers.
    
    Args:
        x: X coordinate
        y: Y coordinate
        dem_tif: Rasterio dataset reader for DEM
        
    Returns:
        Tuple[List[float], List[float], List[float]]: X, Y, and Z arrays for interpolation
    """
    x_index, y_index = rasterio.transform.rowcol(dem_tif.transform, x, y)
    xc, yc = rasterio.transform.xy(dem_tif.transform, x_index, y_index)
    
    if x > xc and y > yc:
        x_index_array = [x_index - 1, x_index - 1, x_index, x_index]
        y_index_array = [y_index, y_index + 1, y_index, y_index + 1]
    elif x > xc and y < yc:
        x_index_array = [x_index, x_index, x_index + 1, x_index + 1]
        y_index_array = [y_index, y_index + 1, y_index, y_index + 1]
    elif x < xc and y > yc:
        x_index_array = [x_index - 1, x_index - 1, x_index, x_index]
        y_index_array = [y_index - 1, y_index, y_index - 1, y_index]
    elif x < xc and y < yc:
        x_index_array = [x_index, x_index, x_index + 1, x_index + 1]
        y_index_array = [y_index - 1, y_index, y_index - 1, y_index]
    
    x_array, y_array = rasterio.transform.xy(dem_tif.transform, x_index_array, y_index_array)
    coords = [(lon, lat) for lon, lat in zip(x_array, y_array)]
    z_array = [elev[0] for elev in dem_tif.sample(coords)]
    
    return x_array, y_array, z_array


def dggs_elevation_cdem(x: float, y: float, dem_tif: rasterio.DatasetReader, 
                       vertical_res: int, interp: str = 'linear') -> float:
    """
    Resample CDEM elevation for DGGS cell.
    
    If an error is raised then return -32767 as its final elevation.
    If the point or any of its neighbors has the value -32767 then return -32767 as its final elevation.
    If none of its neighbor has value -32767 then interpolate elevation.
    Restrict the decimal places according to the look-up table defined earlier.
    
    Args:
        x: X coordinate
        y: Y coordinate
        dem_tif: Rasterio dataset reader for DEM
        vertical_res: Vertical resolution for rounding
        interp: Interpolation method (default: 'linear')
        
    Returns:
        float: Interpolated elevation value
    """
    try:
        x_array, y_array, z_array = find_neighbor(x, y, dem_tif)
        
        if -32767 in z_array:
            return -32767
        
        cdem_interp = interpolate.interp2d(x_array, y_array, z_array, kind=interp)
        elevation = cdem_interp(x, y)[0]
        elevation = round(elevation, vertical_res)
        
        return elevation
        
    except Exception:
        return -32767


def dggs_elevation_df(dataframe: pd.DataFrame, dem_tif: rasterio.DatasetReader, 
                     vertical_res: int) -> pd.DataFrame:
    """
    Calculate DGGS elevation for all cells in dataframe.
    
    Args:
        dataframe: DataFrame containing coordinates
        dem_tif: Rasterio dataset reader for DEM
        vertical_res: Vertical resolution for rounding
        
    Returns:
        pd.DataFrame: DataFrame with calculated elevation column
    """
    dataframe['model_elev'] = [
        dggs_elevation_cdem(lon, lat, dem_tif, vertical_res) 
        for lon, lat in zip(dataframe.lon_c, dataframe.lat_c)
    ]
    dataframe = dataframe.drop(columns=['lon_c', 'lat_c'])
    
    return dataframe


def main():
    """Main execution function."""
    if len(sys.argv) != 3:
        logger.error("Usage: python DgQuantization.py <resolution> <area>")
        sys.exit(1)
    
    # Set resolution level and study area
    dggs_res = int(sys.argv[1])  # 20-24
    area = sys.argv[2]  # Calgary / Canmore / BuffaloLake
    
    # Validate inputs
    valid_areas = ['Calgary', 'Canmore', 'BuffaloLake']
    if area not in valid_areas:
        logger.error(f"Invalid area: {area}. Must be one of {valid_areas}")
        sys.exit(1)
    
    # Look up cell size and vertical resolution
    look_up = DgBaseFunc.look_up_table()
    dggs_cellsize = look_up.loc[dggs_res, 'cell_size']
    vertical_res = look_up.loc[dggs_res, 'verti_res']
    
    logger.info(f"Processing: Resolution={dggs_res}, Area={area}")
    
    # Read input data
    input_csv_path = f'Result/Area_{area}_centroids_{dggs_res}.csv'
    if not Path(input_csv_path).exists():
        logger.error(f"Centroid file not found: {input_csv_path}")
        sys.exit(1)
    
    centroid_df = pd.read_csv(input_csv_path, sep=',')
    
    # Read DEM
    dem_path = f'Data/{area}.tif'
    if not Path(dem_path).exists():
        logger.error(f"DEM file not found: {dem_path}")
        sys.exit(1)
    
    logger.info("Processing elevation quantization...")
    start_time = time.time()
    
    # Process elevation quantization using parallel processing
    with rasterio.open(dem_path) as cdem_tif:
        n_cores = int(os.environ.get('SLURM_CPUS_PER_TASK', default=1))
        centroid_df_split = np.array_split(centroid_df, n_cores)
        
        with mp.Pool(processes=n_cores) as pool:
            centroid_df_output = pd.concat(pool.map(
                lambda df: dggs_elevation_df(df, cdem_tif, vertical_res), 
                centroid_df_split
            ))
    
    processing_time = time.time() - start_time
    logger.info(f"Processing time: {processing_time:.2f} seconds")
    
    # Save the results as csv
    output_csv_path = f'Result/{area}_elev_{dggs_res}.csv'
    centroid_df_output.to_csv(output_csv_path, index=False)
    logger.info(f"Saved elevation data to: {output_csv_path}")
    
    logger.info("Processing completed successfully!")


if __name__ == '__main__':
    main()
