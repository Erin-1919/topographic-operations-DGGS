"""
Discrete Global Grid System (DGGS) Zonal Statistics

This module provides zonal statistics calculations for hexagonal DGGS,
including raster projection and statistical aggregation by zones.

Author: Mingke Li
"""

import rasterio
import gdal
import pandas as pd
import numpy as np
import time
import sys
import os
import multiprocess as mp
import logging
from typing import List, Tuple, Optional, Union, Any
from pathlib import Path

import DgBaseFunc

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def reverse_projection(input_dem: str, output_dem: str) -> None:
    """
    Reverse project raster to NAD83 CSRS coordinate system.
    
    Args:
        input_dem: Path to input DEM file
        output_dem: Path to output projected DEM file
    """
    try:
        with rasterio.open(input_dem) as raster:
            wrap_option = gdal.WarpOptions(
                format='GTiff',
                outputType=gdal.GDT_Float32,
                srcSRS=raster.meta.get('crs'),
                dstSRS='EPSG:4617',  # NAD83(CSRS)
                dstNodata=raster.meta.get('nodata'),
                creationOptions=['COMPRESS=LZW']
            )
            gdal.Warp(output_dem, input_dem, options=wrap_option)
            logger.info(f"Successfully projected {input_dem} to {output_dem}")
    except Exception as e:
        logger.error(f"Error projecting raster: {e}")
        raise


def reclassify_df(dataframe: pd.DataFrame, vege_tif: rasterio.DatasetReader) -> pd.DataFrame:
    """
    Resample raster by nearest neighbor method for dataframe coordinates.
    
    Args:
        dataframe: DataFrame containing coordinates
        vege_tif: Rasterio dataset reader for vegetation raster
        
    Returns:
        pd.DataFrame: DataFrame with vegetation class column added
    """
    coords = [(lon, lat) for lon, lat in zip(dataframe.lon_c, dataframe.lat_c)]
    dataframe['vege_class'] = [DgBaseFunc.catch(lambda: vege[0]) for vege in vege_tif.sample(coords)]
    return dataframe


def calculate_zonal_stats(elev_df: pd.DataFrame, centroid_df: pd.DataFrame, 
                         cell_area: float) -> pd.DataFrame:
    """
    Calculate zonal statistics according to vegetation zones.
    
    Args:
        elev_df: DataFrame containing elevation data
        centroid_df: DataFrame containing centroid data with vegetation classes
        cell_area: Cell area in square meters
        
    Returns:
        pd.DataFrame: DataFrame with calculated statistics by zone
    """
    # Filter out no-data values
    elev_df = elev_df[elev_df['model_elev'] != -32767]
    
    # Join elevation and centroid data
    centroid_df_join = pd.merge(
        left=elev_df, 
        right=centroid_df, 
        how="inner", 
        on=['i', 'j']
    )
    
    # Calculate statistics by vegetation class
    stats_df = centroid_df_join.groupby(["vege_class"]).agg({
        'model_elev': ['mean', 'max', 'min', 'median', 'std', 'sum']
    })
    
    # Flatten column names
    stats_df.columns = ['elev_mean', 'elev_max', 'elev_min', 'elev_median', 'elev_std', 'elev_sum']
    
    # Handle NaN values in standard deviation
    stats_df['elev_std'] = stats_df['elev_std'].fillna(0)
    
    # Calculate additional statistics
    stats_df['elev_range'] = stats_df['elev_max'] - stats_df['elev_min']
    stats_df['volume'] = stats_df['elev_sum'] * cell_area / 1000000000  # Convert to km³
    
    return stats_df


def main():
    """Main execution function."""
    if len(sys.argv) != 3:
        logger.error("Usage: python DgZonalStats.py <resolution> <area>")
        sys.exit(1)
    
    # Set resolution level
    dggs_res = int(sys.argv[1])  # 20-24
    area = sys.argv[2]  # Calgary / Canmore / BuffaloLake
    
    # Validate inputs
    valid_areas = ['Calgary', 'Canmore', 'BuffaloLake']
    if area not in valid_areas:
        logger.error(f"Invalid area: {area}. Must be one of {valid_areas}")
        sys.exit(1)
    
    # Look up cell area
    look_up = DgBaseFunc.look_up_table()
    cell_area = look_up.loc[dggs_res, 'cell_area'] * 1000000
    
    logger.info(f"Processing: Resolution={dggs_res}, Area={area}")
    
    # Read input data
    centroid_csv_path = f'Result/Area_{area}_centroids_{dggs_res}.csv'
    elev_csv_path = f'Result/{area}_elev_{dggs_res}.csv'
    
    if not Path(centroid_csv_path).exists():
        logger.error(f"Centroid file not found: {centroid_csv_path}")
        sys.exit(1)
    
    if not Path(elev_csv_path).exists():
        logger.error(f"Elevation file not found: {elev_csv_path}")
        sys.exit(1)
    
    centroid_df = pd.read_csv(centroid_csv_path, sep=',')
    elev_df = pd.read_csv(elev_csv_path, sep=',')
    
    # Process raster projection
    input_tif = 'Data/aci_2020_ab.tif'
    output_tif = 'Data/aci_2020_ab_NAD83.tif'
    
    if not Path(input_tif).exists():
        logger.error(f"Input vegetation raster not found: {input_tif}")
        sys.exit(1)
    
    logger.info("Projecting vegetation raster...")
    reverse_projection(input_tif, output_tif)
    
    # Process vegetation classification
    logger.info("Processing vegetation classification...")
    start_time = time.time()
    
    # Extract raster values by nearest method using parallel processing
    with rasterio.open(output_tif) as vege_tif:
        n_cores = int(os.environ.get('SLURM_CPUS_PER_TASK', default=1))
        centroid_df_split = np.array_split(centroid_df, n_cores)
        
        with mp.Pool(processes=n_cores) as pool:
            centroid_df_output = pd.concat(pool.map(
                lambda df: reclassify_df(df, vege_tif), centroid_df_split
            ))
    
    processing_time = time.time() - start_time
    logger.info(f"Vegetation classification processing time: {processing_time:.2f} seconds")
    
    # Calculate zonal statistics
    logger.info("Calculating zonal statistics...")
    stats_df = calculate_zonal_stats(elev_df, centroid_df_output, cell_area)
    
    # Save the results as csv
    output_csv_path = f'Result/{area}_elev_zonal_{dggs_res}.csv'
    stats_df.to_csv(output_csv_path, index=True)
    logger.info(f"Saved zonal statistics to: {output_csv_path}")
    
    logger.info("Processing completed successfully!")


if __name__ == '__main__':
    main()
    
