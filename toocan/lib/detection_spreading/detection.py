# =============================================================================
# File        : detect.py
# Author      : Thomas Fiolleau
# Date        : 2025-07-04
# Description : Implements hierarchical detection of cloud systems for TOOCAN.
#               Includes temperature thresholding, and C-accelerated labeling
#               to detect mesoscale convective systems.
#
# Functions   :
#   - load_c_library()   : portable loading of the C shared library (label.so)
#   - prepare_cloudy_indices() : build seed mask and linear indices for C code
#   - detect_objects()   : calls C Label_region to label convective objects
#
# Project     : TOOCAN - Tracking Organized Deep Convection
# License     : (à préciser : MIT, CeCILL, etc.)
# =============================================================================

import numpy as np
import ctypes
import platform
from importlib import resources

from toocan.lib.struct.data_param import DataParam  # C struct definition
import toocan.lib.detection_spreading as ds_pkg     # package containing label.so


# =============================================================================
# Portable loader for label.so / label.dylib / label.dll
# =============================================================================

def load_c_library(lib_basename: str = "label") -> ctypes.CDLL:
    """
    Loads the low-level C library (label.so / dll / dylib) in a portable way.

    - Works after pip install
    - Works on Linux, macOS, Windows
    - No absolute path required

    Parameters
    ----------
    lib_basename : str
        Basename of the library, without extension (e.g. "label").

    Returns
    -------
    ctypes.CDLL
        Loaded C library.
    """
    system = platform.system()

    if system == "Linux":
        filename = f"{lib_basename}.so"
    elif system == "Darwin":  # macOS
        filename = f"{lib_basename}.dylib"
    else:  # Windows
        filename = f"{lib_basename}.dll"

    # Path inside the `toocan.lib.detection_spreading` package
    lib_path = resources.files(ds_pkg) / filename

    if not lib_path.exists():
        raise FileNotFoundError(f"❌ C library not found: {lib_path}")

    return ctypes.CDLL(str(lib_path))


# Load default C library at import time
liblabel = load_c_library("label")


# =============================================================================
# C function signature for Label_region
# =============================================================================

liblabel.Label_region.argtypes = [
    DataParam,                            # struct (passed by value)
    ctypes.c_void_p,                      # clusters (pointer to Blob array)
    ctypes.POINTER(ctypes.c_float),       # imIR
    ctypes.POINTER(ctypes.c_int),         # imlabel
    ctypes.POINTER(ctypes.c_byte),        # imseg
    ctypes.POINTER(ctypes.c_float),       # imsurf
    ctypes.POINTER(ctypes.c_ulong),       # indice_CloudyPix
    ctypes.POINTER(ctypes.c_ulong),       # nbPix_ConvSeed (in/out)
    ctypes.c_int,                         # NSEEDS (nb_ConvSeeds)
    ctypes.c_int,                         # labelMin
    ctypes.c_double,                      # threshold1
    ctypes.c_double                       # threshold2
]

liblabel.Label_region.restype = ctypes.c_int


# =============================================================================
# Constants
# =============================================================================

MIN_VOXELS = 50       # not used here (C handles filtering), kept for future use
MIN_TIMESTEPS = 3     # idem, placeholder for high-level temporal filtering


# =============================================================================
# Helper: build seed mask and linear indices
# =============================================================================

def prepare_cloudy_indices(volume_bt: np.ndarray, threshold: float):
    """
    Builds a mask with value -127 for 'convective candidate' pixels and returns
    their linear indices in C-order.

    Parameters
    ----------
    volume_bt : np.ndarray
        3D brightness temperature array (T, Y, X), float32.
    threshold : float
        BT threshold (e.g., 235 K) to select candidate convective pixels.

    Returns
    -------
    mask : np.ndarray (int8)
        3D mask with -127 for candidate pixels, 0 elsewhere.
    indice_CloudyPix : np.ndarray (uint64)
        1D array of linear indices (C-order) of candidate pixels.
    nbPix_ConvSeed : int
        Number of candidate pixels.
    seuil1 : float
        Main BT threshold (for C code).
    seuil2 : float
        Secondary threshold (kept at 0.0 here).
    """
    # Dimensions (assumées Z,Y,X)
    ZSIZE, YSIZE, XSIZE = volume_bt.shape

    # Création du masque
    mask = np.zeros_like(volume_bt, dtype=np.int8)
    valid = (volume_bt > 0) & (volume_bt <= threshold)
    mask[valid] = -127

    # Indices 3D
    zz, yy, xx = np.where(mask == -127)

    # Conversion en indices linéaires (style C, C-order)
    indice_CloudyPix = np.sort(
        (zz * XSIZE * YSIZE + yy * XSIZE + xx).astype(np.uint64)
    )

    # Nombre de pixels valides
    nbPix_ConvSeed = indice_CloudyPix.size

    seuil1 = float(threshold)
    seuil2 = 0.0

    return mask, indice_CloudyPix, nbPix_ConvSeed, seuil1, seuil2


# =============================================================================
# Main detection entry point
# =============================================================================

def detect_objects(
    data_param: DataParam,
    clusters,
    volume_bt: np.ndarray,
    surface_area_2d: np.ndarray,
    global_label_volume: np.ndarray,
    threshold: float,
    area_threshold_km2: float,
    kernel,
    nb_ConvSeeds: int,
    labelMin: int,
):
    """
    Detect new cloud systems (objects) below a temperature threshold and update
    the global_label_volume using the C-accelerated Label_region() function.

    Parameters
    ----------
    data_param : DataParam
        Configuration / geometry / thresholds structure passed also to C.
    clusters : ctypes array of Blob
        Array of Blob structs used by the C detection/spreading core.
    volume_bt : np.ndarray
        3D array of brightness temperature (T, Y, X) in K.
    surface_area_2d : np.ndarray
        2D array of pixel area in km² (Y, X), already cropped.
    global_label_volume : np.ndarray
        3D label array (T, Y, X) holding current labeled systems.
    threshold : float
        Current brightness temperature threshold for detection (e.g. 235 K).
    area_threshold_km2 : float
        Minimum area per timestep in km² (not yet applied here, handled in C or later).
    kernel : any
        Structuring element for 3D connectivity (currently not used here directly).
    nb_ConvSeeds : int
        Current number of convective seeds (first available label index in C).
    labelMin : int
        Minimum label index used in overlap handling.

    Returns
    -------
    global_label_volume : np.ndarray
        Updated 3D label array (T, Y, X).
    nb_ConvSeeds : int
        Updated number of convective seeds after C labeling.
    """

    # Ensure contiguous float32 volume
    volume_bt = np.ascontiguousarray(volume_bt, dtype=np.float32)

    # === Step 1 : build seed mask and linear indices ===
    imseg, indice_CloudyPix, nbPix_ConvSeed, seuil1, seuil2 = prepare_cloudy_indices(
        volume_bt, threshold
    )

    T, Y, X = volume_bt.shape

    # Ensure float32 + C-contiguous surface area
    if surface_area_2d.dtype != np.float32:
        surface_area_2d = surface_area_2d.astype(np.float32, copy=False)
    if not surface_area_2d.flags['C_CONTIGUOUS']:
        surface_area_2d = np.ascontiguousarray(surface_area_2d)

    # Update geometry in data_param for C code
    data_param.ZSIZE = int(T)
    data_param.YSIZE = int(Y)
    data_param.XSIZE = int(X)

    # --- Build C pointers ---
    ptr_volume_bt = volume_bt.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    ptr_imlabel = global_label_volume.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    ptr_imseg = imseg.ctypes.data_as(ctypes.POINTER(ctypes.c_byte))
    ptr_imsurf = surface_area_2d.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    ptr_indice_CloudyPix = indice_CloudyPix.ctypes.data_as(ctypes.POINTER(ctypes.c_ulong))
    ptr_nbPix = ctypes.pointer(ctypes.c_ulong(nbPix_ConvSeed))

    # Cast clusters array to void* for C
    ptr_clusters = ctypes.cast(clusters, ctypes.c_void_p)

    # --- Call C routine ---
    
    nb_ConvSeeds = liblabel.Label_region(
        data_param,
        ptr_clusters,
        ptr_volume_bt,
        ptr_imlabel,
        ptr_imseg,
        ptr_imsurf,
        ptr_indice_CloudyPix,
        ptr_nbPix,
        nb_ConvSeeds,
        labelMin,
        seuil1,
        seuil2,
    )

    # Reshape just in case (usually no-op)
    global_label_volume = global_label_volume.reshape(
        (data_param.ZSIZE, data_param.YSIZE, data_param.XSIZE)
    )

    return global_label_volume, nb_ConvSeeds























    #n_filtered = nseeds_updated - labelMin




#    # 3. Label connected components (3D clusters)
#    labeled_new, n_new = label(mask, structure=kernel)
#
#
#    print(n_voxels,n_new)
#    # Step 3: Compute neighborhood of already labeled objects
#    #label_neighborhood = binary_dilation(global_label_volume > 0, structure=kernel)
#
#    # 4. Prepare mask for filtered (valid) objects
#    valid_mask = np.zeros_like(labeled_new, dtype=bool)
#
#    for lbl in range(1, n_new + 1):
#        # 5. Extract boolean mask of current label across time
#        time_slices = (labeled_new == lbl)
#        
#        # 6. Compute area (in km²) per time step
#        area_per_timestep = np.array([
#            np.sum(surface_area_2d[time_slices[t]]) for t in range(volume_bt.shape[0])
#        ])
#        #area_per_timestep = np.sum(time_slices * surface_area_2d, axis=(1, 2))
#
#        # 7. Identify which time steps exceed minimum area
#        time_steps_ok = (area_per_timestep >= area_threshold_km2)
#        # 8. Check if object is valid for at least 3 time steps (non-consecutive)
#        if np.sum(time_steps_ok) >= 3:
#            valid_mask |= time_slices
#
#        # 8. Check if object lasts at least 3 **consecutive** time steps
#        #convolved = np.convolve(time_steps_ok.astype(int), np.ones(3, dtype=int), mode='valid')
#        #if np.any(convolved == 3):
#        #    valid_mask |= time_slices  # Keep this object
#
#    # 9. Relabel the surviving objects
#    labeled_filtered, n_filtered = label(valid_mask, structure=kernel)


