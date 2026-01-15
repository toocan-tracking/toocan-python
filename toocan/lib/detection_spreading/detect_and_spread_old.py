# =========================================================================================
# File        : detect_and_spread.py
# Author      : Thomas Fiolleau
# Date        : 2025-07-04
# Description : Main detection and spreading logic for TOOCAN cloud segmentation.
#               Performs hierarchical thresholding and temporal-spatial growth.
#
# Project     : TOOCAN - Tracking Of Organized Convection Algorithm using 3D segmentatioN 
# Repository  : https://github.com/yourusername/toocan_project
# License     : 
# ========================================================================================


import time
import numpy as np
from skimage.morphology import remove_small_objects
from scipy.ndimage import label

from toocan.detection_spreading.spreading import spread_labels
from toocan.detection_spreading.kernel import get_custom_kernel_3d
from toocan.detection_spreading.detection import detect_new_objects
from toocan.detection_spreading.detection import detect_objects
from toocan.utils.array_utils import compute_area_per_timestep

import numpy as np
import pandas as pd
import numpy as np
import pandas as pd
import numpy as np
import pandas as pd

def extract_system_properties(global_label_volume, surface_area_2d, lon_array, lat_array):
    """
    Extracts basic properties of labeled cloud systems from a 3D label volume.

    Parameters:
        global_label_volume (np.ndarray): 3D array (time, lat, lon) with labeled cloud systems.
        surface_area_2d (np.ndarray): 2D array of pixel areas (lat, lon) in km².
        lat_array (np.ndarray): 1D array of latitude values (size Y).
        lon_array (np.ndarray): 1D array of longitude values (size X).

    Returns:
        pd.DataFrame: DataFrame with one row per label and associated properties.
    """
    T, Y, X = global_label_volume.shape
    labels = np.unique(global_label_volume)
    labels = labels[labels > 0]  # Exclude background

    results = []

    for label in labels:
        mask = (global_label_volume == label)

        # Temporal extent
        time_presence = np.any(np.any(mask, axis=2), axis=1)
        t_indices = np.where(time_presence)[0]
        t_start = int(t_indices[0])
        t_end = int(t_indices[-1])
        duration = len(t_indices)

        # Spatial bounding box
        indices = np.argwhere(mask)
        min_t, min_y, min_x = np.min(indices, axis=0)
        max_t, max_y, max_x = np.max(indices, axis=0)

        # Area and pixel count
        n_pixels = np.sum(mask)
        total_area = sum(np.sum(surface_area_2d[mask[t]]) for t in t_indices)

        # Geolocation and area of bottom-right bounding pixel
        lat_max = lat_array[max_y, max_x]
        lon_max = lon_array[max_y, max_x]
        pixel_area_max = surface_area_2d[max_y, max_x]

        results.append({
            "label": label,
            "t_start_index": t_start,
            "t_end_index": t_end,
            "duration_timesteps": duration,
            "n_pixels": n_pixels,
            "area_km2": total_area,
            "minY": min_y,
            "maxY": max_y,
            "minX": min_x,
            "maxX": max_x,
            "lat_maxY": lat_max,
            "lon_maxX": lon_max,
            "pixel_area_km2_max": pixel_area_max
        })

    return pd.DataFrame(results)

def detect_and_spread(data_param,volume_bt, surface_area_2d, lon_array_2d,lat_array_2d,global_label_volume,
                                           params_TOOCAN, next_label_id=1):
    """
    Hierarchical labeling with temporal and size filtering of new detections using physical area threshold.
    
    Parameters:
        volume_bt: 3D ndarray of brightness temperature (time, lat, lon)
        surface_area_2d: 2D ndarray (lat, lon)
        global_label_volume: 3D ndarray (same shape as volume_bt), initialized by caller
        next_label_id: int, value to start labeling from (useful for chaining chunks)
    """
    import time
    from skimage.morphology import remove_small_objects
    from scipy.ndimage import label

    start0 = time.time()
    kernel = get_custom_kernel_3d()
    shape = volume_bt.shape
    area_threshold_km2 = float(params_TOOCAN.get("minAreaSeed"))
    minBT_threshold = int(params_TOOCAN.get("minBT_threshold"))
    maxBT_threshold = int(params_TOOCAN.get("maxBT_threshold"))
    stepBT_threshold = int(params_TOOCAN.get("stepBT_threshold"))
    deltaBT_Spread  = int(params_TOOCAN.get("deltaBT_Spread"))
    nCluster = 0
    print("Volume shape:", volume_bt.shape)
    print("Min/Max BT:", np.nanmin(volume_bt), np.nanmax(volume_bt))
     
   
    #print("Bottom latitude:", np.min(volume_bt),np.max(volume_bt))
    print(np.shape(lat_array_2d))
    for threshold in range(minBT_threshold, maxBT_threshold + 1, stepBT_threshold):
        
        threshold_IntermediateCloudshield = threshold + stepBT_threshold
        if(threshold + stepBT_threshold > maxBT_threshold):
            threshold_IntermediateCloudshield = maxBT_threshold
        print(f"\n>>> Threshold {threshold}K - Threshold Intermediate Cloudshield {threshold_IntermediateCloudshield}K")
        start = time.time()

        # Step 0: Re-expand existing labels from previous chunk 
        print("  Re-expanding all existing objects...")
        tmask = volume_bt < (threshold)
        cloud_mask = tmask.astype(np.uint8)

        global_label_volume = spread_labels(
            volume_bt       = volume_bt,
            labeled_volume  = global_label_volume,
            cloud_mask      = cloud_mask,
            delta_spread    = deltaBT_Spread
        )
 
        # Step 0: Re-expand existing labels from previous chunk 
        print("  Re-expanding all existing objects...")

        # Step 1: Detect new objects
        labeled_new_filtered, n_new_filtered = detect_objects(
            volume_bt=volume_bt,
            surface_area_2d=surface_area_2d,
            global_label_volume=global_label_volume,
            threshold=threshold,
            area_threshold_km2=area_threshold_km2,
            kernel=kernel
        )

        # Step 0: Re-expand existing labels from previous chunk 
        print("label OK")
        
        if n_new_filtered > 0:
            labeled_new_filtered[labeled_new_filtered > 0] += next_label_id - 1
            global_label_volume[labeled_new_filtered > 0] = labeled_new_filtered[labeled_new_filtered > 0]
            next_label_id = global_label_volume.max() + 1
            print(f"  New objects kept after filtering: {n_new_filtered}")
        else:
            print("  No new objects met the temporal/area condition.")
            
        nCluster=nCluster+n_new_filtered
        print(f"  nb of clusters: {nCluster}")


        # After detection and labeling is done:
        #df_props = extract_system_properties(
        #    global_label_volume=global_label_volume,
        #    surface_area_2d=surface_area_2d, lon_array=lon_array_2d, lat_array=lat_array_2d
        #)

        #print(df_props)


        # Step 2: Re-expand existing labels
        print("  Re-expanding all existing objects...")
        tmask = volume_bt <= (threshold_IntermediateCloudshield)
        cloud_mask = tmask.astype(np.uint8)
      
        global_label_volume = spread_labels(
            volume_bt=volume_bt,
            labeled_volume=global_label_volume,
            cloud_mask=cloud_mask,
            delta_spread=deltaBT_Spread
        )
  
        # After detection and labeling is done:
        #df_props = extract_system_properties(
        #    global_label_volume=global_label_volume,
        #    surface_area_2d=surface_area_2d, lon_array=lon_array_2d, lat_array=lat_array_2d
        #)
        #print(df_props)

        end = time.time()
        print(f"Runtime: {end - start:.2f} seconds")

    end0 = time.time()
    print(f"Total Runtime: {end0 - start0:.2f} seconds")
    return global_label_volume, next_label_id