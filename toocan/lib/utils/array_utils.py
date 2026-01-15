def compute_index(t, y, x, Y, X):
    """
    Compute the flattened index in a 1D raveled array 
    corresponding to coordinates (t, y, x) in a 3D array (T, Y, X).

    Parameters:
        t (int): Time index
        y (int): Y (latitude) index
        x (int): X (longitude) index
        Y (int): Total size in Y dimension
        X (int): Total size in X dimension

    Returns:
        int: Flattened index
    """
    return t * (Y * X) + y * X + x


def valid_coords(t, y, x, T, Y, X):
    """
    Check whether a 3D coordinate (t, y, x) is within bounds.

    Parameters:
        t (int): Time index
        y (int): Y (latitude) index
        x (int): X (longitude) index
        T (int): Time dimension size
        Y (int): Y dimension size
        X (int): X dimension size

    Returns:
        bool: True if coordinates are within bounds
    """
    return 0 <= t < T and 0 <= y < Y and 0 <= x < X


def compute_area_per_timestep(mask3d, surface_area_2d):
    import numpy as np
    return np.array([
        np.sum(surface_area_2d[mask3d[t]]) for t in range(mask3d.shape[0])
    ])

def normalize_array(arr):
    arr_min = np.nanmin(arr)
    arr_max = np.nanmax(arr)
    return (arr - arr_min) / (arr_max - arr_min + 1e-6)

def pad_array(arr, pad_width, pad_value=0):
    return np.pad(arr, pad_width=pad_width, mode='constant', constant_values=pad_value)