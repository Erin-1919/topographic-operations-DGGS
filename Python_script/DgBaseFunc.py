"""
Discrete Global Grid System (DGGS) Base Functions

This module provides fundamental functions for working with hexagonal DGGS,
including neighbor calculations, coordinate transformations, and mathematical operations.

Author: Mingke Li
"""

import pandas as pd
import numpy as np
import scipy.linalg
import math
from typing import List, Tuple, Optional, Union, Any
import warnings


def look_up_table() -> pd.DataFrame:
    """
    Construct a look-up table for DGGS resolution parameters.
    
    Returns:
        pd.DataFrame: DataFrame containing resolution, cell size, cell area, 
                     vertical resolution, and cell spacing for each DGGS resolution level.
    """
    res_list = [16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29]
    cell_size_list = [0.005, 0.003, 0.001, 0.0009, 0.0008, 0.0003, 0.0003, 
                      0.0001, 0.0001, 0.00006, 0.00003, 0.00002, 0.00001, 0.000005]
    cell_spacing_list = [1.0750880097, 0.6207023518, 0.3583626699, 0.2069007839, 
                         0.1194542233, 0.0689669280, 0.0398180744, 0.0229889760, 
                         0.0132726915, 0.0076629920, 0.0044242305, 0.0025543307, 
                         0.0014747435, 0.0008514436]
    cell_area_list = [1.1849116724224, 0.3949705574741, 0.1316568524914, 0.0438856174971, 
                      0.0146285391657, 0.0048761797219, 0.0016253932406, 0.0005417977469, 
                      0.0001805992490, 0.0000601997497, 0.0000200665832, 0.0000066888611, 
                      0.0000022296204, 0.0000007432068]
    vertical_res_list = [0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6]
    
    look_up = pd.DataFrame({
        'res': res_list, 
        'cell_size': cell_size_list, 
        'verti_res': vertical_res_list, 
        'cell_area': cell_area_list, 
        'cell_spacing': cell_spacing_list
    }, index=res_list)
    
    return look_up


def check_if_duplicates(list_of_elements: List[Any]) -> bool:
    """
    Check if given list contains any duplicates.
    
    Args:
        list_of_elements: List to check for duplicates
        
    Returns:
        bool: True if duplicates exist, False otherwise
    """
    return len(list_of_elements) != len(set(list_of_elements))


def common_member(a: List[Any], b: List[Any]) -> set:
    """
    Check if two lists have at least one element in common.
    
    Args:
        a: First list
        b: Second list
        
    Returns:
        set: Set of common elements between the two lists
    """
    return set(a) & set(b)


def catch(func: callable, handle: callable = lambda e: e, *args, **kwargs) -> Any:
    """
    Handle try-except in a general function.
    
    Args:
        func: Function to execute
        handle: Error handler function
        *args: Positional arguments for func
        **kwargs: Keyword arguments for func
        
    Returns:
        Any: Result of func or handle function
    """
    try:
        return func(*args, **kwargs)
    except Exception:
        return np.nan


def edge_cell_exist(elevation_list: List[float]) -> bool:
    """
    Determine if edge cell exists in the neighborhood.
    
    Args:
        elevation_list: List of elevation values (center + 6 neighbors)
        
    Returns:
        bool: True if edge cell exists (less than 7 valid cells)
    """
    cell_list = [x for x in elevation_list if x != -32767 and not np.isnan(x)]
    return len(cell_list) != 7


def edge_cell(res: int, coords: Tuple[int, int], all_cells: set) -> bool:
    """
    Check if a cell is an edge cell.
    
    Args:
        res: DGGS resolution
        coords: Cell coordinates (i, j)
        all_cells: Set of all valid cell coordinates
        
    Returns:
        bool: True if cell is an edge cell
    """
    neighbors = neighbor_coords(res, coords)[1:]
    neighbor_in = [i in all_cells for i in neighbors]
    neighbor_filt = [x for x in neighbor_in if x is True]
    return len(neighbor_filt) != 6


def edge_cell_chain(res: int, coords: Tuple[int, int], all_cells: set) -> List[Tuple[int, int]]:
    """
    Find chain of edge cells.
    
    Args:
        res: DGGS resolution
        coords: Starting cell coordinates
        all_cells: Set of all valid cell coordinates
        
    Returns:
        List[Tuple[int, int]]: Chain of edge cell coordinates
    """
    chain = [coords]
    next_cell = np.nan
    
    while next_cell != coords:
        current = chain[-1]
        neighbors = neighbor_coords(res, current)[1:]
        neighbor_in = [i for i in neighbors if i in all_cells]
        
        for neighbor in neighbor_in:
            if edge_cell(res, neighbor, all_cells) and (neighbor not in chain[1:]):
                next_cell = neighbor
                if next_cell != coords:
                    chain.append(next_cell)
                break
                
    return chain


def neighbor_coords(res: int, coords: Tuple[int, int]) -> List[Tuple[int, int]]:
    """
    Find the ij coordinates within the neighborhood (6 neighbors + center cell).
    
    Args:
        res: DGGS resolution
        coords: Center cell coordinates (i, j)
        
    Returns:
        List[Tuple[int, int]]: List of neighbor coordinates including center
    """
    i, j = coords[0], coords[1]
    
    if res % 2 == 0:
        index_i = [i, i-1, i-1, i, i+1, i+1, i]
        index_j = [j, j-1, j, j+1, j+1, j, j-1]
    else:  # res % 2 == 1
        index_i = [i, i-2, i-1, i+1, i+2, i+1, i-1]
        index_j = [j, j-1, j+1, j+2, j+1, j-1, j-2]
        
    coords_neighbor = [(ih, jh) for ih, jh in zip(index_i, index_j)]
    return coords_neighbor


def neighbor_navig_by_rings(res: int, coords: Tuple[int, int], df: pd.DataFrame, rings: int = 1) -> List[float]:
    """
    Navigate among neighborhood defined by ring numbers.
    
    Args:
        res: DGGS resolution
        coords: Center cell coordinates
        df: DataFrame containing elevation data
        rings: Number of rings to expand (default: 1)
        
    Returns:
        List[float]: Elevation values of neighbors + center cell
    """
    i = 1
    coords_edge = [coords]
    coords_neighbor = [coords]
    
    while i <= rings:
        i += 1
        coords_neighbor_prev = coords_neighbor.copy()
        
        for c in coords_edge:
            coords_neighbor.extend(neighbor_coords(res, c)[1:])
            
        coords_neighbor = list(set(coords_neighbor))
        coords_edge = [x for x in coords_neighbor if x not in coords_neighbor_prev]
        
    elev_neighbor = [catch(lambda: df['model_elev'].loc[ij]) for ij in coords_neighbor]
    return elev_neighbor


def edge_navig_by_rings(res: int, coords: Tuple[int, int], rings: int = 1) -> List[Tuple[int, int]]:
    """
    Navigate among neighborhood defined by ring numbers, return coordinates only.
    
    Args:
        res: DGGS resolution
        coords: Center cell coordinates
        rings: Number of rings to expand (default: 1)
        
    Returns:
        List[Tuple[int, int]]: Coordinates of edge cells
    """
    i = 1
    coords_edge = [coords]
    coords_neighbor = [coords]
    
    while i <= rings:
        i += 1
        coords_neighbor_prev = coords_neighbor.copy()
        
        for c in coords_edge:
            coords_neighbor.extend(neighbor_coords(res, c)[1:])
            
        coords_neighbor = list(set(coords_neighbor))
        coords_edge = [x for x in coords_neighbor if x not in coords_neighbor_prev]
        
    return coords_edge


def neighbor_navig(res: int, coords: Tuple[int, int], df: pd.DataFrame) -> List[float]:
    """
    Extract elevation values among the neighborhood (1 center + 6 neighbors).
    
    Args:
        res: DGGS resolution
        coords: Center cell coordinates
        df: DataFrame containing elevation data
        
    Returns:
        List[float]: Elevation values of center and neighbors
    """
    coords_neighbor = neighbor_coords(res, coords)
    elev_neighbor = [catch(lambda: df['model_elev'].loc[ij]) for ij in coords_neighbor]
    return elev_neighbor


def neighbor_slope_navig(res: int, coords: Tuple[int, int], df: pd.DataFrame) -> List[float]:
    """
    Extract slope values among the neighborhood (1 center + 6 neighbors).
    
    Args:
        res: DGGS resolution
        coords: Center cell coordinates
        df: DataFrame containing slope data
        
    Returns:
        List[float]: Slope values of center and neighbors
    """
    coords_neighbor = neighbor_coords(res, coords)
    slp_neighbor = [catch(lambda: df['gradient_deg'].loc[ij]) for ij in coords_neighbor]
    return slp_neighbor


def neighbor_direc_navig(res: int, coords: Tuple[int, int], df: pd.DataFrame) -> List[float]:
    """
    Extract direction codes among the neighborhood (1 center + 6 neighbors).
    
    Args:
        res: DGGS resolution
        coords: Center cell coordinates
        df: DataFrame containing direction data
        
    Returns:
        List[float]: Direction codes of center and neighbors
    """
    coords_neighbor = neighbor_coords(res, coords)
    direc_neighbor = [catch(lambda: df['direction_code'].loc[ij]) for ij in coords_neighbor]
    return direc_neighbor


def neighbor_multi_direc_navig(res: int, coords: Tuple[int, int], df: pd.DataFrame) -> List[float]:
    """
    Extract multi-direction codes among the neighborhood (1 center + 6 neighbors).
    
    Args:
        res: DGGS resolution
        coords: Center cell coordinates
        df: DataFrame containing multi-direction data
        
    Returns:
        List[float]: Multi-direction codes of center and neighbors
    """
    coords_neighbor = neighbor_coords(res, coords)
    direc_multi_neighbor = [catch(lambda: df['direction_code'].loc[ij]) for ij in coords_neighbor]
    return direc_multi_neighbor


def cross_product(p0: Tuple[float, float, float], 
                  p1: Tuple[float, float, float], 
                  p2: Tuple[float, float, float]) -> np.ndarray:
    """
    Calculate the cross product of two vectors composed by three points.
    
    Args:
        p0: First point (x, y, z)
        p1: Second point (x, y, z)
        p2: Third point (x, y, z)
        
    Returns:
        np.ndarray: Cross product vector
    """
    x0, y0, z0 = p0
    x1, y1, z1 = p1
    x2, y2, z2 = p2
    
    # Vector u
    ux, uy, uz = x1 - x0, y1 - y0, z1 - z0
    # Vector v
    vx, vy, vz = x2 - x0, y2 - y0, z2 - z0
    
    u_cross_v = np.array([uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx])
    
    if u_cross_v[2] < 0:
        u_cross_v = u_cross_v * (-1)
        
    return u_cross_v


def aspect_restricted(res: int) -> List[int]:
    """
    Get restricted aspect lists based on resolution.
    
    Args:
        res: DGGS resolution
        
    Returns:
        List[int]: List of restricted aspect angles
    """
    if res % 2 == 1:
        asp_list = [30, 90, 150, 210, 270, 330]
    else:  # res % 2 == 0
        asp_list = [0, 60, 120, 180, 240, 300]
    return asp_list


def aspect_restricted_oppo(res: int) -> List[int]:
    """
    Get restricted aspect in opposite directions.
    
    Args:
        res: DGGS resolution
        
    Returns:
        List[int]: List of restricted aspect angles in opposite directions
    """
    if res % 2 == 1:
        asp_list = [210, 270, 330, 30, 90, 150]
    else:  # res % 2 == 0
        asp_list = [180, 240, 300, 0, 60, 120]
    return asp_list


def aspect_unrestricted(dzx: float, dzy: float) -> float:
    """
    Calculate the unrestricted direction to which that cell is oriented.
    Algorithm from Hodgson 1998.
    
    Args:
        dzx: X-component of gradient
        dzy: Y-component of gradient
        
    Returns:
        float: Aspect angle in radians
    """
    if dzx == 0:
        if dzy < 0:
            asp = math.pi
        elif dzy > 0:
            asp = 0
        else:
            asp = -1
    elif dzx > 0:
        asp = math.pi / 2 - math.atan(dzy / dzx)
    else:
        asp = math.pi / 2 * 3 - math.atan(dzy / dzx)
    return asp


def mean_norm_vector(res: int, z_list: List[float], *cells: int) -> np.ndarray:
    """
    Calculate the mean of norm vectors, vector count > 1.
    
    Args:
        res: DGGS resolution
        z_list: List of elevation values
        *cells: Cell indices to include in calculation
        
    Returns:
        np.ndarray: Mean normal vector
    """
    if res % 2 == 1:
        x_list = [0, math.sqrt(3)/3, math.sqrt(3)/3*2, math.sqrt(3)/3, 
                  math.sqrt(3)/-3, math.sqrt(3)/-3*2, math.sqrt(3)/-3]
        y_list = [0, 1, 0, -1, -1, 0, 1]
    else:  # res % 2 == 0
        x_list = [0, 0, 1, 1, 0, -1, -1]
        y_list = [0, math.sqrt(3)/3*2, math.sqrt(3)/3, math.sqrt(3)/-3, 
                  math.sqrt(3)/-3*2, math.sqrt(3)/-3, math.sqrt(3)/3]
        
    pt_c = [x_list[0], y_list[0], z_list[0]]
    points = []
    x_sum = y_sum = z_sum = 0
    
    for i in cells:
        points.append([x_list[i], y_list[i], z_list[i]])
        
    for i, pt in enumerate(points[:-1]):
        norm = cross_product(pt_c, pt, points[i + 1])
        x_sum += norm[0]
        y_sum += norm[1]
        z_sum += norm[2]
        
    return np.array([x_sum, y_sum, z_sum])


def fit_plane_norm_vector(res: int, z_list: List[float]) -> np.ndarray:
    """
    Fit a plane by using least squares and find out its norm vector.
    
    Args:
        res: DGGS resolution
        z_list: List of elevation values
        
    Returns:
        np.ndarray: Normal vector of fitted plane
    """
    if res % 2 == 1:
        x_list = [0, math.sqrt(3)/3, math.sqrt(3)/3*2, math.sqrt(3)/3, 
                  math.sqrt(3)/-3, math.sqrt(3)/-3*2, math.sqrt(3)/-3]
        y_list = [0, 1, 0, -1, -1, 0, 1]
    else:  # res % 2 == 0
        x_list = [0, 0, 1, 1, 0, -1, -1]
        y_list = [0, math.sqrt(3)/3*2, math.sqrt(3)/3, math.sqrt(3)/-3, 
                  math.sqrt(3)/-3*2, math.sqrt(3)/-3, math.sqrt(3)/3]
        
    G = np.c_[x_list, y_list, z_list]
    A = np.c_[G[:, 0], G[:, 1], np.ones(G.shape[0])]
    C, _, _, _ = scipy.linalg.lstsq(A, G[:, 2], rcond=None)
    norm_vector = np.array([C[0] * -1, C[1] * -1, 1])
    return norm_vector


def first_derivative(res: int, elev_neighbor: List[float]) -> Tuple[float, float]:
    """
    Project the non-normalized gradient to orthogonal x y axes.
    
    Args:
        res: DGGS resolution
        elev_neighbor: List of elevation values (center + 6 neighbors)
        
    Returns:
        Tuple[float, float]: X and Y components of gradient
    """
    if res % 2 == 1:
        d_zi = (elev_neighbor[3] - elev_neighbor[6]) / 2
        d_zj = (elev_neighbor[4] - elev_neighbor[1]) / 2
        d_zk = (elev_neighbor[5] - elev_neighbor[2]) / 2
        d_zx = d_zk + d_zj * math.sin(math.pi/6) - d_zi * math.sin(math.pi/6)
        d_zy = d_zi * math.cos(math.pi/6) + d_zj * math.cos(math.pi/6)
    else:  # res % 2 == 0
        d_zi = (elev_neighbor[4] - elev_neighbor[1]) / 2
        d_zj = (elev_neighbor[5] - elev_neighbor[2]) / 2
        d_zk = (elev_neighbor[6] - elev_neighbor[3]) / 2
        d_zx = d_zj * math.cos(math.pi/6) + d_zk * math.cos(math.pi/6)
        d_zy = d_zi + d_zj * math.sin(math.pi/6) - d_zk * math.sin(math.pi/6)
        
    return d_zx, d_zy


def second_derivative(res: int, elev_neighbor: List[float]) -> Tuple[float, float]:
    """
    Project the non-normalized second derivatives to orthogonal x y axes.
    
    Args:
        res: DGGS resolution
        elev_neighbor: List of elevation values (center + 6 neighbors)
        
    Returns:
        Tuple[float, float]: X and Y components of second derivative
    """
    if res % 2 == 1:
        d_zi2 = 2 * elev_neighbor[0] - elev_neighbor[3] - elev_neighbor[6]
        d_zj2 = 2 * elev_neighbor[0] - elev_neighbor[1] - elev_neighbor[4]
        d_zk2 = 2 * elev_neighbor[0] - elev_neighbor[2] - elev_neighbor[5]
        d_zx2 = d_zk2 + d_zj2 * math.sin(math.pi/6) - d_zi2 * math.sin(math.pi/6)
        d_zy2 = d_zi2 * math.cos(math.pi/6) + d_zj2 * math.cos(math.pi/6)
    else:  # res % 2 == 0
        d_zi2 = 2 * elev_neighbor[0] - elev_neighbor[1] - elev_neighbor[4]
        d_zj2 = 2 * elev_neighbor[0] - elev_neighbor[2] - elev_neighbor[5]
        d_zk2 = 2 * elev_neighbor[0] - elev_neighbor[3] - elev_neighbor[6]
        d_zx2 = d_zj2 * math.cos(math.pi/6) + d_zk2 * math.cos(math.pi/6)
        d_zy2 = d_zi2 + d_zj2 * math.sin(math.pi/6) - d_zk2 * math.sin(math.pi/6)
        
    return d_zx2, d_zy2


if __name__ == '__main__':
    # Test the look-up table function
    lookup = look_up_table()
    print("DGGS Look-up Table:")
    print(lookup.head())
