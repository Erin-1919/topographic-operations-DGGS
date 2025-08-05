"""
Discrete Global Grid System (DGGS) Topographic Analysis

This module provides topographic analysis functions for hexagonal DGGS,
including slope calculation using various algorithms and curvature analysis.

Author: Mingke Li
"""

import pandas as pd
import numpy as np
import math
import time
import sys
import os
import warnings
import functools
import multiprocess as mp
import logging
from typing import List, Tuple, Optional, Union, Any
from pathlib import Path

import DgBaseFunc as dbfc

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

warnings.simplefilter('error', RuntimeWarning)


def slope_MAG(coords: Tuple[int, int], df: pd.DataFrame, res: int, cell_spacing: float) -> Tuple[float, float]:
    """
    Calculate slope using Maximum Absolute Gradient (MAG) method.
    
    Absolute maximum differences between the center cell and its six neighbours.
    If edge cell, assign nan values.
    If flat, slope = 0; aspect = -1.
    If not flat, slope angle is the normalized gradient, slope aspect is one of the restricted aspect values.
    If multiple equal gradient found, then choose the first neighbour encountered clockwise from north.
    If the slope is uphill, the downhill aspect is considered to be in the directly opposite direction.
    
    Args:
        coords: Cell coordinates (i, j)
        df: DataFrame containing elevation data
        res: DGGS resolution
        cell_spacing: Cell spacing in meters
        
    Returns:
        Tuple[float, float]: (gradient_mag, aspect_mag) in degrees
    """
    elev_neighbor = dbfc.neighbor_navig(res, coords, df)
    
    if dbfc.edge_cell_exist(elev_neighbor):
        return np.nan, np.nan
    
    gradient = [e - elev_neighbor[0] for e in elev_neighbor]
    gradient_abs = [abs(e) for e in gradient]
    gradient_max = max(gradient_abs)
    
    if gradient_max == 0:
        gradient_mag = 0
        aspect_mag = -1
    else:
        gradient_rad_mag = math.atan(gradient_max / cell_spacing)
        gradient_mag = math.degrees(gradient_rad_mag)
        aspect_index = gradient_abs.index(gradient_max)
        
        if gradient[aspect_index] > 0:
            aspect_mag = dbfc.aspect_restricted_oppo(res)[aspect_index - 1]
        elif gradient[aspect_index] < 0:
            aspect_mag = dbfc.aspect_restricted(res)[aspect_index - 1]
    
    return gradient_mag, aspect_mag


def slope_MDG(coords: Tuple[int, int], df: pd.DataFrame, res: int, cell_spacing: float) -> Tuple[float, float]:
    """
    Calculate slope using Maximum Downward Gradient (MDG) method.
    
    Maximum downward gradient between the center cell and its six neighbours.
    Need pit-filling beforehand without saving the altered elevations.
    If edge cell, assign nan values.
    If flat, slope = 0; aspect = -1.
    If not flat, slope angle is the normalized gradient, slope aspect is one of the restricted aspect values.
    If multiple equal gradient found, then choose the first neighbour encountered clockwise from north.
    
    Args:
        coords: Cell coordinates (i, j)
        df: DataFrame containing elevation data
        res: DGGS resolution
        cell_spacing: Cell spacing in meters
        
    Returns:
        Tuple[float, float]: (gradient_mdg, aspect_mdg) in degrees
    """
    elev_neighbor = dbfc.neighbor_navig(res, coords, df)
    
    if dbfc.edge_cell_exist(elev_neighbor):
        return np.nan, np.nan
    
    elev_neighbor_ = elev_neighbor[1:]
    if all(i > elev_neighbor[0] for i in elev_neighbor_):
        elev_neighbor[0] = min(elev_neighbor_)
    
    gradient = [e - elev_neighbor[0] for e in elev_neighbor]
    gradient_min = min(gradient)
    
    if gradient_min >= 0:
        gradient_mdg = 0
        aspect_mdg = -1
    else:
        gradient_rad_mdg = math.atan(abs(gradient_min) / cell_spacing)
        gradient_mdg = math.degrees(gradient_rad_mdg)
        aspect_index = gradient.index(gradient_min)
        aspect_mdg = dbfc.aspect_restricted(res)[aspect_index - 1]
    
    return gradient_mdg, aspect_mdg


def slope_MDN(coords: Tuple[int, int], df: pd.DataFrame, res: int, cell_spacing: float) -> Tuple[float, float]:
    """
    Calculate slope using Multiple Downhill Neighbours (MDN) method.
    
    Multiple downhill neighbours - distributing flow from a pixel amongst all of its lower elevation neighbor pixels.
    Need pit-filling beforehand without saving the altered elevations.
    If edge cell, assign nan values.
    If flat, slope = 0; aspect = -1.
    If not flat, slope and aspect are calculated by taking the average of all the downslope neighbors.
    Aspect is represented by the mean norm vector which is a "mean" surface orientation.
    
    Args:
        coords: Cell coordinates (i, j)
        df: DataFrame containing elevation data
        res: DGGS resolution
        cell_spacing: Cell spacing in meters
        
    Returns:
        Tuple[float, float]: (gradient_mdn, aspect_mdn) in degrees
    """
    elev_neighbor = dbfc.neighbor_navig(res, coords, df)
    
    if dbfc.edge_cell_exist(elev_neighbor):
        return np.nan, np.nan
    
    elev_neighbor_ = elev_neighbor[1:]
    if all(i > elev_neighbor[0] for i in elev_neighbor_):
        elev_neighbor[0] = min(elev_neighbor_)
    
    gradient = [e - elev_neighbor[0] for e in elev_neighbor]
    
    if min(gradient) >= 0:
        gradient_mdn = 0
        aspect_mdn = -1
    else:
        gradient_rad_mdn_ls = [math.atan(abs(g) / cell_spacing) for g in gradient if g < 0]
        gradient_rad_mdn = sum(gradient_rad_mdn_ls) / len(gradient_rad_mdn_ls)
        gradient_mdn = math.degrees(gradient_rad_mdn)
        
        cell_ls = []
        for i, g in enumerate(gradient):
            if g < 0:
                cell_ls.append(i)
        
        if len(cell_ls) == 1:
            aspect_index = cell_ls[0]
            aspect_mdn = dbfc.aspect_restricted(res)[aspect_index - 1]
        elif len(cell_ls) == 2 and abs(cell_ls[0] - cell_ls[1]) == 3:
            aspect_index_1, aspect_index_2 = cell_ls[0], cell_ls[1]
            gradient_1, gradient_2 = gradient[aspect_index_1], gradient[aspect_index_2]
            aspect_index = aspect_index_1 if gradient_1 <= gradient_2 else aspect_index_2
            aspect_mdn = dbfc.aspect_restricted(res)[aspect_index - 1]
        else:
            avg_norm = dbfc.mean_norm_vector(res, elev_neighbor, *cell_ls)
            aspect_mdn = dbfc.aspect_unrestricted(avg_norm[0], avg_norm[1])
            aspect_mdn = math.degrees(aspect_mdn)
        
        if aspect_mdn < 0:
            aspect_mdn += 360
    
    return gradient_mdn, aspect_mdn


def slope_FDA(coords: Tuple[int, int], df: pd.DataFrame, res: int, cell_spacing: float) -> Tuple[float, float]:
    """
    Calculate slope using Finite-Difference Algorithm (FDA).
    
    Project the non-normalized gradient in three directions to orthogonal x y axes.
    If edge cell, assign nan values.
    If flat, slope = 0; aspect = -1.
    If not flat, calculate slope and aspect by using finite difference algorithms 
    and combine the two component partial derivatives.
    
    Args:
        coords: Cell coordinates (i, j)
        df: DataFrame containing elevation data
        res: DGGS resolution
        cell_spacing: Cell spacing in meters
        
    Returns:
        Tuple[float, float]: (gradient_fda, aspect_fda) in degrees
    """
    elev_neighbor = dbfc.neighbor_navig(res, coords, df)
    
    if dbfc.edge_cell_exist(elev_neighbor):
        return np.nan, np.nan
    
    dzx, dzy = dbfc.first_derivative(res, elev_neighbor)
    
    if dzx == 0 and dzy == 0:
        gradient_fda = 0
        aspect_fda = -1
    else:
        gradient_rad_fda = math.atan(math.sqrt(dzx**2 + dzy**2) / cell_spacing)
        gradient_fda = math.degrees(gradient_rad_fda)
        aspect_fda = dbfc.aspect_unrestricted(dzx, dzy)
        aspect_fda = math.degrees(aspect_fda)
    
    if aspect_fda < 0:
        aspect_fda += 360
    
    return gradient_fda, aspect_fda


def slope_BFP(coords: Tuple[int, int], df: pd.DataFrame, res: int, cell_spacing: float) -> Tuple[float, float]:
    """
    Calculate slope using Best Fit Plane (BFP) method.
    
    If edge cell, assign nan values.
    If flat, slope = 0; aspect = -1.
    If not flat, a surface is fitted to seven centroids by multiple linear regression models, 
    using least squares to minimize the sum of distances from the surface to the cells.
    
    Args:
        coords: Cell coordinates (i, j)
        df: DataFrame containing elevation data
        res: DGGS resolution
        cell_spacing: Cell spacing in meters
        
    Returns:
        Tuple[float, float]: (gradient_bfp, aspect_bfp) in degrees
    """
    elev_neighbor = dbfc.neighbor_navig(res, coords, df)
    
    if dbfc.edge_cell_exist(elev_neighbor):
        return np.nan, np.nan
    
    try:
        norm_vec = np.array(dbfc.fit_plane_norm_vector(res, elev_neighbor))
        norm_vec_mag = np.linalg.norm(norm_vec)
        unit_norm_vec = np.array([i / norm_vec_mag for i in norm_vec])
        ref_vec = np.array([0, 0, -1])
        ref_vec_proj = ref_vec - (np.dot(ref_vec, unit_norm_vec) * unit_norm_vec)
        
        gradient_rad_bfp = math.atan((ref_vec_proj[2]**2 / math.sqrt(ref_vec_proj[0]**2 + ref_vec_proj[1]**2)) / cell_spacing)
        gradient_bfp = math.degrees(gradient_rad_bfp)
        aspect_bfp = dbfc.aspect_unrestricted(ref_vec_proj[0], ref_vec_proj[1])
        aspect_bfp = math.degrees(aspect_bfp)
        
        if aspect_bfp < 0:
            aspect_bfp += 360
            
    except Exception:
        gradient_bfp = 0
        aspect_bfp = -1
    
    return gradient_bfp, aspect_bfp


def curvature(coords: Tuple[int, int], df: pd.DataFrame, res: int, cell_spacing: float) -> float:
    """
    Calculate curvature (slope rate of change of landform).
    
    Second derivative of DEM, first derivative of slope.
    
    Args:
        coords: Cell coordinates (i, j)
        df: DataFrame containing elevation data
        res: DGGS resolution
        cell_spacing: Cell spacing in meters
        
    Returns:
        float: Curvature value
    """
    elev_neighbor = dbfc.neighbor_navig(res, coords, df)
    
    if dbfc.edge_cell_exist(elev_neighbor):
        return np.nan
    
    dzx2, dzy2 = dbfc.second_derivative(res, elev_neighbor)
    curv = math.sqrt(dzx2**2 + dzy2**2) / cell_spacing**2
    return curv


def hillshade(coords: Tuple[int, int], df: pd.DataFrame, res: int, 
               altitude: float = 45, azimuth: float = 315) -> float:
    """
    Calculate hillshade value.
    
    A hillshade is a grayscale 3D representation of the surface, with the sun's relative position 
    taken into account for shading. This function uses the altitude and azimuth properties to 
    specify the sun's position.
    
    Azimuth is the angular direction of the sun, measured from north in clockwise degrees from 0 to 360. 
    An azimuth of 90 degrees is east. The default azimuth is 315 degrees (NW).
    Altitude is the slope or angle of the illumination source above the horizon. 
    The units are in degrees, from 0 (on the horizon) to 90 (overhead). The default is 45 degrees.
    
    Args:
        coords: Cell coordinates (i, j)
        df: DataFrame containing elevation data
        res: DGGS resolution
        altitude: Sun altitude in degrees (default: 45)
        azimuth: Sun azimuth in degrees (default: 315)
        
    Returns:
        float: Hillshade value (0-255)
    """
    slope_deg, aspect_deg = slope_FDA(coords, df, res, 1000)  # Using default cell_spacing
    slope_rad = math.radians(slope_deg)
    aspect_rad = math.radians(aspect_deg)
    
    if slope_rad == 0:
        hs = 255.0
    else:
        zenith_deg = 90.0 - altitude
        zenith_rad = math.radians(zenith_deg)
        azimuth_math = 360.0 - azimuth + 90.0
        
        if azimuth_math >= 360.0:
            azimuth_math = azimuth_math - 360.0
            
        azimuth_rad = math.radians(azimuth_math)
        hs = 255.0 * ((math.cos(zenith_rad) * math.cos(slope_rad)) + 
                       (math.sin(zenith_rad) * math.sin(slope_rad) * math.cos(azimuth_rad - aspect_rad)))
    
    return hs


def slope_aspect_df(dataframe: pd.DataFrame, elev_df: pd.DataFrame, method: str, 
                   res: int, cell_spacing: float) -> pd.DataFrame:
    """
    Calculate slope and aspect by specified method.
    
    Args:
        dataframe: DataFrame to process
        elev_df: DataFrame containing elevation data
        method: Slope calculation method ('MAG', 'MDG', 'MDN', 'FDA', 'BFP')
        res: DGGS resolution
        cell_spacing: Cell spacing in meters
        
    Returns:
        pd.DataFrame: DataFrame with calculated slope and aspect columns
    """
    method_functions = {
        'MAG': slope_MAG,
        'MDG': slope_MDG,
        'MDN': slope_MDN,
        'FDA': slope_FDA,
        'BFP': slope_BFP
    }
    
    if method not in method_functions:
        raise ValueError(f"Invalid method: {method}. Must be one of {list(method_functions.keys())}")
    
    func = method_functions[method]
    dataframe[['gradient_deg', 'aspect_deg']] = [
        func(ij, elev_df, res, cell_spacing) for ij in dataframe.index.values
    ]
    
    return dataframe


def curvature_df(dataframe: pd.DataFrame, elev_df: pd.DataFrame, res: int, cell_spacing: float) -> pd.DataFrame:
    """
    Calculate curvature for all cells in dataframe.
    
    Args:
        dataframe: DataFrame to process
        elev_df: DataFrame containing elevation data
        res: DGGS resolution
        cell_spacing: Cell spacing in meters
        
    Returns:
        pd.DataFrame: DataFrame with calculated curvature column
    """
    dataframe['curv'] = [curvature(ij, elev_df, res, cell_spacing) for ij in dataframe.index.values]
    return dataframe


def hillshade_df(dataframe: pd.DataFrame, elev_df: pd.DataFrame, res: int) -> pd.DataFrame:
    """
    Calculate hillshade for all cells in dataframe.
    
    Args:
        dataframe: DataFrame to process
        elev_df: DataFrame containing elevation data
        res: DGGS resolution
        
    Returns:
        pd.DataFrame: DataFrame with calculated hillshade column
    """
    dataframe['hs'] = [hillshade(ij, elev_df, res) for ij in dataframe.index.values]
    return dataframe


def main():
    """Main execution function."""
    if len(sys.argv) != 4:
        logger.error("Usage: python DgTopogAnaly.py <resolution> <area> <method>")
        sys.exit(1)
    
    # Set resolution level
    dggs_res = int(sys.argv[1])
    area = sys.argv[2]  # Calgary / Canmore / BuffaloLake
    method = sys.argv[3]  # MAG, MDG, MDN, FDA, BFP
    
    # Validate inputs
    valid_methods = ['MAG', 'MDG', 'MDN', 'FDA', 'BFP']
    if method not in valid_methods:
        logger.error(f"Invalid method: {method}. Must be one of {valid_methods}")
        sys.exit(1)
    
    # Look up cell spacing and vertical resolution
    look_up = dbfc.look_up_table()
    cell_spacing = look_up.loc[dggs_res, 'cell_spacing'] * 1000
    vertical_res = look_up.loc[dggs_res, 'verti_res']
    
    logger.info(f"Processing: Resolution={dggs_res}, Area={area}, Method={method}")
    
    # Read input data
    input_csv_path = f'Result/{area}_elev_{dggs_res}.csv'
    if not Path(input_csv_path).exists():
        logger.error(f"Input file not found: {input_csv_path}")
        sys.exit(1)
    
    elev_df = pd.read_csv(input_csv_path, sep=',')
    elev_df = elev_df.set_index(['i', 'j'])
    
    # Process slope and aspect
    logger.info("Processing slope and aspect...")
    start_time = time.time()
    
    elev_df_copy = elev_df.copy()
    slope_aspect_df_p = functools.partial(slope_aspect_df, elev_df=elev_df, 
                                         method=method, res=dggs_res, cell_spacing=cell_spacing)
    
    n_cores = int(os.environ.get('SLURM_CPUS_PER_TASK', default=1))
    elev_df_split = np.array_split(elev_df_copy, n_cores)
    
    with mp.Pool(processes=n_cores) as pool:
        elev_df_output = pd.concat(pool.map(slope_aspect_df_p, elev_df_split))
    
    processing_time = time.time() - start_time
    logger.info(f"Slope and aspect processing time: {processing_time:.2f} seconds")
    
    # Save results
    output_csv_path = f'Result/{area}_elev_{method}_{dggs_res}.csv'
    elev_df_output.to_csv(output_csv_path, index=True)
    logger.info(f"Saved slope and aspect results to: {output_csv_path}")
    
    # Process curvature
    logger.info("Processing curvature...")
    start_time = time.time()
    
    elev_df_copy = elev_df.copy()
    curvature_df_p = functools.partial(curvature_df, elev_df=elev_df, 
                                      res=dggs_res, cell_spacing=cell_spacing)
    
    with mp.Pool(processes=n_cores) as pool:
        elev_df_output = pd.concat(pool.map(curvature_df_p, elev_df_split))
    
    processing_time = time.time() - start_time
    logger.info(f"Curvature processing time: {processing_time:.2f} seconds")
    
    # Save curvature results
    output_csv_path = f'Result/{area}_curvature_{dggs_res}.csv'
    elev_df_output.to_csv(output_csv_path, index=True)
    logger.info(f"Saved curvature results to: {output_csv_path}")
    
    # Process hillshade
    logger.info("Processing hillshade...")
    start_time = time.time()
    
    elev_df_copy = elev_df.copy()
    hillshade_df_p = functools.partial(hillshade_df, elev_df=elev_df, res=dggs_res)
    
    with mp.Pool(processes=n_cores) as pool:
        elev_df_output = pd.concat(pool.map(hillshade_df_p, elev_df_split))
    
    processing_time = time.time() - start_time
    logger.info(f"Hillshade processing time: {processing_time:.2f} seconds")
    
    # Save hillshade results
    output_csv_path = f'Result/{area}_hillshade_{dggs_res}.csv'
    elev_df_output.to_csv(output_csv_path, index=True)
    logger.info(f"Saved hillshade results to: {output_csv_path}")
    
    logger.info("Processing completed successfully!")


if __name__ == '__main__':
    main()
