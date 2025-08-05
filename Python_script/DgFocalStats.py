"""
Discrete Global Grid System (DGGS) Focal Statistics

This module provides focal statistics calculations for hexagonal DGGS,
including descriptive statistics for neighborhood analysis.

Author: Mingke Li
"""

import pandas as pd
import numpy as np
import statistics
import time
import sys
import os
import multiprocess as mp
import logging
from typing import List, Tuple, Optional, Union, Any
from pathlib import Path
from itertools import product

import DgBaseFunc

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def focal_elev_stats(coords: Tuple[int, int], df: pd.DataFrame, rings: int = 1, 
                    res: int = 20) -> List[float]:
    """
    Calculate focal elevation statistics for a cell.
    
    Remove voids in neighbors and calculate descriptive statistics:
    mean, max, min, median, std, range.
    
    Args:
        coords: Cell coordinates (i, j)
        df: DataFrame containing elevation data
        rings: Number of rings to expand (default: 1)
        res: DGGS resolution (default: 20)
        
    Returns:
        List[float]: List of statistics [mean, max, min, median, std, range]
    """
    elev_neighbor = DgBaseFunc.neighbor_navig_by_rings(res, coords, df, rings)
    elev_neighbor = [x for x in elev_neighbor if x != -32767 and not np.isnan(x)]
    
    if len(elev_neighbor) == 0:
        return [np.nan] * 6
    
    if len(elev_neighbor) == 1:
        return elev_neighbor * 4 + [0] * 2
    
    # Calculate statistics
    elev_mean = round(statistics.mean(elev_neighbor), 3)
    elev_max = round(max(elev_neighbor), 3)
    elev_min = round(min(elev_neighbor), 3)
    elev_median = round(statistics.median(elev_neighbor), 3)
    elev_std = round(statistics.stdev(elev_neighbor), 3)
    elev_range = round(elev_max - elev_min, 3)
    
    return [elev_mean, elev_max, elev_min, elev_median, elev_std, elev_range]


def focal_elev_stats_df(dataframe: pd.DataFrame, elev_df: pd.DataFrame, 
                       rings: int = 1, res: int = 20) -> pd.DataFrame:
    """
    Calculate focal elevation statistics for all cells in dataframe.
    
    Args:
        dataframe: DataFrame to process
        elev_df: DataFrame containing elevation data
        rings: Number of rings to expand (default: 1)
        res: DGGS resolution (default: 20)
        
    Returns:
        pd.DataFrame: DataFrame with calculated focal statistics columns
    """
    dataframe[['mean', 'max', 'min', 'median', 'std', 'range']] = [
        focal_elev_stats(ij, elev_df, rings, res) for ij in dataframe.index.values
    ]
    return dataframe


def main():
    """Main execution function."""
    if len(sys.argv) != 4:
        logger.error("Usage: python DgFocalStats.py <resolution> <ring_number> <area>")
        sys.exit(1)
    
    # Set resolution level, study area, and ring number
    dggs_res = int(sys.argv[1])  # 20-24
    ring_n = int(sys.argv[2])  # suggest 1-3
    area = sys.argv[3]  # Calgary / Canmore / BuffaloLake
    
    # Validate inputs
    valid_areas = ['Calgary', 'Canmore', 'BuffaloLake']
    if area not in valid_areas:
        logger.error(f"Invalid area: {area}. Must be one of {valid_areas}")
        sys.exit(1)
    
    if ring_n < 1 or ring_n > 5:
        logger.warning(f"Ring number {ring_n} is outside recommended range (1-3)")
    
    # Look up cell size and vertical resolution
    look_up = DgBaseFunc.look_up_table()
    dggs_cellsize = look_up.loc[dggs_res, 'cell_size']
    vertical_res = look_up.loc[dggs_res, 'verti_res']
    
    logger.info(f"Processing: Resolution={dggs_res}, Rings={ring_n}, Area={area}")
    
    # Read input data
    input_csv_path = f'Result/{area}_elev_{dggs_res}.csv'
    if not Path(input_csv_path).exists():
        logger.error(f"Input file not found: {input_csv_path}")
        sys.exit(1)
    
    elev_df = pd.read_csv(input_csv_path, sep=',')
    elev_df = elev_df.set_index(['i', 'j'])
    
    # Process focal statistics
    logger.info("Processing focal statistics...")
    start_time = time.time()
    
    # Call the function by parallel processing
    n_cores = int(os.environ.get('SLURM_CPUS_PER_TASK', default=1))
    elev_df_split = np.array_split(elev_df, n_cores)
    
    with mp.Pool(processes=n_cores) as pool:
        elev_df_output = pd.concat(pool.starmap(
            focal_elev_stats_df, 
            product(elev_df_split, [ring_n] * len(elev_df_split))
        ))
    
    processing_time = time.time() - start_time
    logger.info(f"Processing time: {processing_time:.2f} seconds")
    
    # Save the results as csv
    output_csv_path = f'Result/{area}_elev_focal_{dggs_res}_ring{ring_n}.csv'
    elev_df_output.to_csv(output_csv_path, index=True)
    logger.info(f"Saved focal statistics to: {output_csv_path}")
    
    logger.info("Processing completed successfully!")


if __name__ == '__main__':
    main()
