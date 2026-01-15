# =============================================================================
# File        : detect.py
# Author      : Thomas Fiolleau
# Date        : 2025-07-04
# Description : Implements hierarchical detection of cloud systems for TOOCAN.
#               Includes temperature thresholding, spatial-temporal filtering,
#               and minimum area/persistence validation to detect mesoscale convective systems.
#
# Functions   :
#   - detect_new_objects(): Identifies new cloud systems using BT thresholds,
#                           removes small/noisy detections, filters by area/time.
#
# Project     : TOOCAN - Tracking Organized Deep Convection
# License     : MIT (or specify your license)
# =============================================================================
import numpy as np
import ctypes
from scipy.ndimage import convolve

import os
import ctypes

lib_path = os.path.join(os.path.dirname(__file__), "label_sparse_voxels.so")
labeler = ctypes.CDLL(lib_path)

labeler.label_sparse_voxels.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.c_int,
    ctypes.c_int, ctypes.c_int, ctypes.c_int,
    ctypes.POINTER(ctypes.c_int)
]
labeler.label_sparse_voxels.restype = ctypes.c_int


import numpy as np
from scipy.ndimage import convolve, binary_dilation
from .kernel import get_custom_kernel_3d  # Ensure you import your kernel function

MIN_VOXELS = 50  # minimum object size in voxels
MIN_TIMESTEPS = 3  # minimum number of time steps the object must persist

def run_labeling_with_c(coords, shape):
    import ctypes

    n_voxels = coords.shape[0]
    if n_voxels == 0:
        return np.zeros(0, dtype=np.int32)

    coords_flat = coords.astype(np.int32).flatten()
    labels_out = np.zeros(n_voxels, dtype=np.int32)

    coords_c = coords_flat.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    labels_c = labels_out.ctypes.data_as(ctypes.POINTER(ctypes.c_int))

    ret = labeler.label_sparse_voxels(
        coords_c,
        ctypes.c_int(n_voxels),
        ctypes.c_int(shape[0]), ctypes.c_int(shape[1]), ctypes.c_int(shape[2]),
        labels_c
    )
    if ret != 0:
        raise RuntimeError(f"C labeling returned error code {ret}")

    return labels_out


def detect_new_objects(volume_bt, surface_area_2d, global_label_volume, threshold, area_threshold_km2, kernel):
    print("🟢 detect_new_objects(): started")

    # Step 1: create binary mask
    mask = (volume_bt > 150) & (volume_bt <= threshold) & (global_label_volume == 0)
    print(f"  ➤ mask created: {np.sum(mask)} voxels")

    # Step 2: remove isolated voxels (<2 neighbors)
    neighbor_count = convolve(mask.astype(int), np.ones((3, 3, 3)), mode='constant')
    mask &= (neighbor_count >= 2)
    print(f"  ➤ after filtering: {np.sum(mask)} voxels")

    # Step 3: extract coordinates
    coords = np.argwhere(mask)
    coords = coords[np.lexsort((coords[:,2], coords[:,1], coords[:,0]))]
    if len(coords) == 0:
        print("🟡 No candidates found")
        return np.zeros_like(volume_bt, dtype=np.int32), 0
    print(volume_bt.shape)
    # Step 4: call fast C labeling
    labels = run_labeling_with_c(coords, volume_bt.shape)
    print(f"  ✅ C labeling done: {np.max(labels)} labels")

    n_clusters = len(np.unique(labels[labels > 0]))
    print(f"  ✅ C labeling done: {n_clusters} clusters")

    # Step 5: create sparse labeled volume
    labeled_sparse = np.zeros_like(volume_bt, dtype=np.int32)
    for idx, (z, y, x) in enumerate(coords):
        labeled_sparse[z, y, x] = labels[idx]

    # Step 6: object-wise filtering
    valid_mask = np.zeros_like(volume_bt, dtype=bool)
    for lbl in np.unique(labels):
        if lbl == 0:
            continue

        time_slices = (labeled_sparse == lbl)

        # Skip small blobs
        if np.sum(time_slices) < MIN_VOXELS:
            continue

        # Skip objects touching already labeled ones
        touching = binary_dilation(time_slices, structure=kernel) & (global_label_volume > 0)
        if np.any(touching):
            continue

        # Compute area per timestep
        try:
            area_per_timestep = np.array([
                np.sum(surface_area_2d[time_slices[t]]) for t in range(volume_bt.shape[0])
            ])
        except Exception as e:
            print(f"  ❌ Area error for label {lbl}: {e}")
            continue

        # Check persistence
        time_steps_ok = (area_per_timestep >= area_threshold_km2)
        if np.sum(time_steps_ok) >= MIN_TIMESTEPS:
            valid_mask |= time_slices

    # Step 7: final relabeling
    final_coords = np.argwhere(valid_mask)
    if len(final_coords) == 0:
        print("🔴 No final objects detected")
        return np.zeros_like(volume_bt, dtype=np.int32), 0

    final_labels = run_labeling_with_c(final_coords, volume_bt.shape)
    final_labeled = np.zeros_like(volume_bt, dtype=np.int32)
    for idx, (z, y, x) in enumerate(final_coords):
        final_labeled[z, y, x] = final_labels[idx]

    print(f"🟢 detect_new_objects(): completed with {np.max(final_labels)} objects")
    return final_labeled, np.max(final_labels)


#
#def run_labeling_with_c(coords, shape):
#    import ctypes
#    import numpy as np
#
#    try:
#        n_voxels = coords.shape[0]
#        if n_voxels == 0:
#            return np.zeros(0, dtype=np.int32)
#
#        # Allocate ctypes arrays
#        coords_flat = coords.astype(np.int32).flatten()
#        labels_out = np.zeros(n_voxels, dtype=np.int32)
#
#        coords_c = coords_flat.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
#        labels_c = labels_out.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
#
#        print(f"    ➤ C labeling with {n_voxels} voxels...")
#        ret = labeler.label_sparse_voxels(
#            coords_c,
#            ctypes.c_int(n_voxels),
#            ctypes.c_int(shape[0]), ctypes.c_int(shape[1]), ctypes.c_int(shape[2]),
#            labels_c
#        )
#        if ret != 0:
#            print(f"    ❌ C function returned error code: {ret}")
#        else:
#            print(f"    ✅ C labeling completed")
#
#        return labels_out
#    except Exception as e:
#        print(f"    ❌ Exception in run_labeling_with_c: {e}")
#        raise
#
#def detect_new_objects(volume_bt, surface_area_2d, global_label_volume, threshold, area_threshold_km2, kernel):
#    print("🟢 detect_new_objects(): started")
#
#    mask = (volume_bt > 150) & (volume_bt <= threshold) & (global_label_volume == 0)
#    print(f"  ➤ mask created: {np.sum(mask)} voxels")
#
#    neighbor_count = convolve(mask.astype(int), np.ones((3,3,3)), mode='constant')
#    mask &= (neighbor_count >= 2)
#    print(f"  ➤ after filtering: {np.sum(mask)} voxels")
#
#    coords = np.argwhere(mask)
#    print(f"  ➤ np.argwhere(mask): {len(coords)} voxels")
#
#    if len(coords) == 0:
#        print("🟡 No candidates found")
#        return np.zeros_like(volume_bt, dtype=np.int32), 0
#
#    # Run C labeling
#    print("  ➤ Calling run_labeling_with_c...")
#    labels = run_labeling_with_c(coords, volume_bt.shape)
#    print(f"  ✅ C labeling done: {np.max(labels)} labels")
#
#    # Fill volume
#    labeled_sparse = np.zeros_like(volume_bt, dtype=np.int32)
#    for idx, (z, y, x) in enumerate(coords):
#        labeled_sparse[z, y, x] = labels[idx]
#
#    print("  ➤ Labeled sparse volume filled")
#
#    # Filtering
#    valid_mask = np.zeros_like(volume_bt, dtype=bool)
#    unique_labels = np.unique(labels)
#    print(f"  ➤ Checking {len(unique_labels)} labels")
#
#    for lbl in unique_labels:
#        if lbl == 0:
#            continue
#        time_slices = (labeled_sparse == lbl)
#
#        try:
#            area_per_timestep = np.array([
#                np.sum(surface_area_2d[time_slices[t]]) for t in range(volume_bt.shape[0])
#            ])
#        except Exception as e:
#            print(f"  ❌ Error in area computation for label {lbl}: {e}")
#            continue
#
#        time_steps_ok = (area_per_timestep >= area_threshold_km2)
#        if np.sum(time_steps_ok) >= 3:
#            valid_mask |= time_slices
#
#    # Final relabeling
#    final_coords = np.argwhere(valid_mask)
#    print(f"  ➤ {len(final_coords)} voxels after filtering")
#
#    if len(final_coords) == 0:
#        print("🔴 No final objects detected")
#        return np.zeros_like(volume_bt, dtype=np.int32), 0
#
#    final_labels = run_labeling_with_c(final_coords, volume_bt.shape)
#    final_labeled = np.zeros_like(volume_bt, dtype=np.int32)
#    for idx, (z, y, x) in enumerate(final_coords):
#        final_labeled[z, y, x] = final_labels[idx]
#
#    print(f"🟢 detect_new_objects(): completed with {np.max(final_labels)} objects")
#    return final_labeled, np.max(final_labels)

import numpy as np
from skimage.morphology import remove_small_objects
from scipy.ndimage import label
from toocan.utils.array_utils import compute_area_per_timestep
from scipy.ndimage import binary_dilation
from scipy.ndimage import convolve

def detect_objects(volume_bt, surface_area_2d, global_label_volume, threshold, area_threshold_km2, kernel):
    """
    Detect new cloud systems (objects) below a temperature threshold
    and filter them using area and temporal persistence criteria.

    Parameters:
        volume_bt (np.ndarray): 3D array of brightness temperature (time, lat, lon)
        surface_area_2d (np.ndarray): 2D array of pixel area in km² (lat, lon)
        global_label_volume (np.ndarray): Current label map (same shape as volume_bt)
        threshold (float): Current brightness temperature threshold for detection
        area_threshold_km2 (float): Minimum required area per timestep (in km²)
        kernel: Structuring element for 3D connectivity (e.g., 10-connectivity)

    Returns:
        labeled_filtered (np.ndarray): Filtered label map of new valid detections
        n_filtered (int): Number of new detections passing filtering
    """
    
    # 1. Create a binary mask where BT < threshold and pixel is not already labeled
    mask = (volume_bt > 150) & (volume_bt <= threshold) & (global_label_volume == 0)

    # 2. Remove isolated small objects (less than ~50 voxels)
    #mask = remove_small_objects(mask.astype(bool), min_size=50, connectivity=3)


    # 2. Remove isolated pixels (neighbors < 3)
    neighbor_count = convolve(mask.astype(int), np.ones((3,3,3)), mode='constant')
    mask &= (neighbor_count >= 2)
    
    coords = np.argwhere(mask)
    n_voxels = coords.shape[0]

    # 3. Label connected components (3D clusters)
    labeled_new, n_new = label(mask, structure=kernel)
    print(n_voxels,n_new)
    # Step 3: Compute neighborhood of already labeled objects
    #label_neighborhood = binary_dilation(global_label_volume > 0, structure=kernel)

    # 4. Prepare mask for filtered (valid) objects
    valid_mask = np.zeros_like(labeled_new, dtype=bool)

    for lbl in range(1, n_new + 1):
        # 5. Extract boolean mask of current label across time
        time_slices = (labeled_new == lbl)
        

        # Reject objects touching previous ones
        #if np.any(label_neighborhood & time_slices):
        #    print("Touching")
        #    continue
        
        # Check if this object touches any existing label
        #touching_mask = binary_dilation(time_slices, structure=kernel) & (global_label_volume > 0)
        #if np.any(touching_mask):
        #    print("!!! Touching")
        #    continue  # Skip this object — it's touching an existing one

        # 6. Compute area (in km²) per time step
        area_per_timestep = np.array([
            np.sum(surface_area_2d[time_slices[t]]) for t in range(volume_bt.shape[0])
        ])
        #area_per_timestep = np.sum(time_slices * surface_area_2d, axis=(1, 2))

        # 7. Identify which time steps exceed minimum area
        time_steps_ok = (area_per_timestep >= area_threshold_km2)
        # 8. Check if object is valid for at least 3 time steps (non-consecutive)
        if np.sum(time_steps_ok) >= 3:
            valid_mask |= time_slices

        # 8. Check if object lasts at least 3 **consecutive** time steps
        #convolved = np.convolve(time_steps_ok.astype(int), np.ones(3, dtype=int), mode='valid')
        #if np.any(convolved == 3):
        #    valid_mask |= time_slices  # Keep this object

    # 9. Relabel the surviving objects
    labeled_filtered, n_filtered = label(valid_mask, structure=kernel)

    return labeled_filtered, n_filtered


import numpy as np
from skimage.morphology import remove_small_objects
from scipy.ndimage import label
from scipy.ndimage import binary_dilation

def detect_new_objects_2(volume_bt, surface_area_2d, global_label_volume, threshold, area_threshold_km2, kernel, min_time_steps=3, n_jobs=1):
    """
    Detect new cloud systems under a temperature threshold, filter them by size and persistence.

    Parameters:
        volume_bt (np.ndarray): (T, H, W) Brightness temperature array
        surface_area_2d (np.ndarray): (H, W) Area in km² per pixel
        global_label_volume (np.ndarray): (T, H, W) Existing labeled systems
        threshold (float): Brightness temperature threshold
        area_threshold_km2 (float): Minimum area required per timestep (in km²)
        kernel: Structuring element for 3D connectivity (e.g., 10-connectivity)
        min_time_steps (int): Minimum # of timesteps the object must meet area condition
        n_jobs (int): Parallel jobs (1 = no parallelism)

    Returns:
        labeled_filtered (np.ndarray): New valid systems (relabeled)
        n_filtered (int): Number of valid detections
    """

    from joblib import Parallel, delayed

    # 1. Threshold-based mask for new cold objects not touching existing labels
    mask = (volume_bt > 150) & (volume_bt <= threshold) & (global_label_volume == 0)

    # 2. Remove small noise (< 50 voxels)
    mask = remove_small_objects(mask.astype(bool), min_size=50, connectivity=3)

    # 3. Label new connected components
    labeled_new, n_new = label(mask, structure=kernel)

    # Preallocate valid mask
    valid_mask = np.zeros_like(labeled_new, dtype=bool)

    def check_label(lbl):
        time_slices = (labeled_new == lbl)  # 3D bool mask for this object
        # Vectorized area calc over time
        area_per_timestep = np.sum(time_slices * surface_area_2d, axis=(1, 2))
        time_steps_ok = (area_per_timestep >= area_threshold_km2)
        if np.sum(time_steps_ok) >= min_time_steps:
            return time_slices
        return None

    # 4. Run check per label (parallel or serial)
    if n_jobs == 1:
        results = [check_label(lbl) for lbl in range(1, n_new + 1)]
    else:
        results = Parallel(n_jobs=n_jobs, prefer="threads")(
            delayed(check_label)(lbl) for lbl in range(1, n_new + 1)
        )

    # 5. Combine all valid masks
    for r in results:
        if r is not None:
            valid_mask |= r

    # 6. Relabel filtered detections for compact labeling
    labeled_filtered, n_filtered = label(valid_mask, structure=kernel)

    return labeled_filtered, n_filtered