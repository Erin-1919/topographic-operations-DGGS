"""
Discrete Global Grid System (DGGS) Topographic Indices

This module provides topographic index calculations for hexagonal DGGS,
including Terrain Roughness Index (TRI) and Topographic Position Index (TPI).

Author: Mingke Li
"""

import pandas as pd
import numpy as np
import math
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


def tri_calculate(coords: Tuple[int, int], res: int, elev_df: pd.DataFrame) -> float:
    """
    Calculate Terrain Roughness Index (TRI).
    
    TRI presents the square root of sum of squared difference between a central cell 
    and the adjacent cells.
    
    Args:
        coords: Cell coordinates (i, j)
        res: DGGS resolution
        elev_df: DataFrame containing elevation data
        
    Returns:
        float: TRI value
    """
    elev_neighbor = DgBaseFunc.neighbor_navig(res, coords, elev_df)
    
    if DgBaseFunc.edge_cell_exist(elev_neighbor):
        return np.nan
    
    tri_values = [(n - elev_neighbor[0]) ** 2 for n in elev_neighbor[1:]]
    tri_value = math.sqrt(sum(tri_values))
    
    return tri_value


def tpi_calculate(coords: Tuple[int, int], res: int, elev_df: pd.DataFrame) -> float:
    """
    Calculate Topographic Position Index (TPI).
    
    TPI presents the difference between a central cell and the mean of its surrounding cells.
    
    Args:
        coords: Cell coordinates (i, j)
        res: DGGS resolution
        elev_df: DataFrame containing elevation data
        
    Returns:
        float: TPI value
    """
    elev_neighbor = DgBaseFunc.neighbor_navig(res, coords, elev_df)
    
    if DgBaseFunc.edge_cell_exist(elev_neighbor):
        return np.nan
    
    neighbor_elevations = elev_neighbor[1:]
    tpi_value = elev_neighbor[0] - (sum(neighbor_elevations) / len(neighbor_elevations))
    
    return tpi_value


def topo_index_df(dataframe: pd.DataFrame, res: int, elev_df: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate topographic indices for all cells in dataframe.
    
    Args:
        dataframe: DataFrame to process
        res: DGGS resolution
        elev_df: DataFrame containing elevation data
        
    Returns:
        pd.DataFrame: DataFrame with calculated TRI and TPI columns
    """
    dataframe['TRI'] = [tri_calculate(ij, res, elev_df) for ij in dataframe.index.values]
    dataframe['TPI'] = [tpi_calculate(ij, res, elev_df) for ij in dataframe.index.values]
    
    return dataframe


def main():
    """Main execution function."""
    if len(sys.argv) != 4:
        logger.error("Usage: python DgTopogIndex.py <resolution> <area> <method>")
        sys.exit(1)
    
    # Set resolution level and study area
    dggs_res = int(sys.argv[1])  # 20-24
    area = sys.argv[2]  # Calgary / Canmore / BuffaloLake
    method = sys.argv[3]  # MAG, MDG, MDN, FDA, BFP
    
    # Validate inputs
    valid_methods = ['MAG', 'MDG', 'MDN', 'FDA', 'BFP']
    if method not in valid_methods:
        logger.error(f"Invalid method: {method}. Must be one of {valid_methods}")
        sys.exit(1)
    
    # Look up cell spacing and vertical resolution
    look_up = DgBaseFunc.look_up_table()
    cell_spacing = look_up.loc[dggs_res, 'cell_spacing'] * 1000
    vertical_res = look_up.loc[dggs_res, 'verti_res']
    
    logger.info(f"Processing: Resolution={dggs_res}, Area={area}, Method={method}")
    
    # Read input data
    input_csv_path = f'Result/{area}_flow_{method}_{dggs_res}.csv'
    if not Path(input_csv_path).exists():
        logger.error(f"Input file not found: {input_csv_path}")
        sys.exit(1)
    
    elev_df = pd.read_csv(input_csv_path, sep=',')
    elev_df = elev_df.set_index(['i', 'j'])
    
    # Process topographic indices
    logger.info("Processing topographic indices...")
    start_time = time.time()
    
    # Call the function by parallel processing
    n_cores = int(os.environ.get('SLURM_CPUS_PER_TASK', default=1))
    elev_df_split = np.array_split(elev_df, n_cores)
    
    with mp.Pool(processes=n_cores) as pool:
        elev_df_output = pd.concat(pool.map(
            lambda df: topo_index_df(df, dggs_res, elev_df), elev_df_split
        ))
    
    processing_time = time.time() - start_time
    logger.info(f"Processing time: {processing_time:.2f} seconds")
    
    # Save the results as csv
    output_csv_path = f'Result/{area}_TopoIndex_{method}_{dggs_res}.csv'
    elev_df_output.to_csv(output_csv_path, index=True)
    logger.info(f"Saved topographic indices to: {output_csv_path}")
    
    logger.info("Processing completed successfully!")


if __name__ == '__main__':
    main()

