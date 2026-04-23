# =========================================================================================
# File        : detect_and_spread.py
# Author      : Thomas Fiolleau
# Date        : 2025-07-04
# Description : Main detection and spreading logic for TOOCAN cloud segmentation.
#               Performs hierarchical thresholding and temporal-spatial growth.
#
# Project     : TOOCAN - Tracking Of Organized Convection Algorithm using 3D segmentatioN
# License     : MIT
# ========================================================================================

import time
import numpy as np
from toocan.lib.detection_spreading.spreading import spread_labels,spread_labels_fast
from toocan.lib.detection_spreading.kernel import get_custom_kernel_3d
from toocan.lib.detection_spreading.detection import detect_objects
from toocan.lib.struct.data_param import DataParam

# (La fonction detect_new_objects reste disponible dans detection.py si besoin)
def list_clusters(global_label_volume):
    unique_labels, counts = np.unique(global_label_volume, return_counts=True)
    
    # Supprime le label 0 (fond / non-cluster)
    valid = unique_labels != 0
    unique_labels = unique_labels[valid]
    counts = counts[valid]
    
    print(f"\n📦 {len(unique_labels)} clusters détectés :\n")
    for label, count in zip(unique_labels, counts):
        print(f" - Label {label} → taille : {count} voxels")
    
    return list(zip(unique_labels, counts))


def detect_and_spread(
    data_param,
    clusters,
    volume_bt,
    surface_area_2d,
    lon_array_2d,
    lat_array_2d,
    global_label_volume,
    params_TOOCAN,
    nb_ConvSeeds,
    labelMin=1
):
    """
    Hierarchical labeling with temporal and size filtering of new detections
    using physical area threshold.
    """
    start_time = time.time()
    kernel = get_custom_kernel_3d()

    area_threshold_km2 = float(params_TOOCAN.get("minAreaSeed"))
    minBT_threshold = int(params_TOOCAN.get("minBT_threshold"))
    maxBT_threshold = int(params_TOOCAN.get("maxBT_threshold"))
    stepBT_threshold = int(params_TOOCAN.get("stepBT_threshold"))
    deltaBT_Spread = int(params_TOOCAN.get("deltaBT_Spread", 1))

    nCluster = 0

    print("Min/Max BT:", np.nanmin(volume_bt), np.nanmax(volume_bt))


    for threshold in range(minBT_threshold, maxBT_threshold + stepBT_threshold, stepBT_threshold):
        threshold_next = min(threshold + stepBT_threshold, maxBT_threshold)

        print(f"\n>>> Threshold {threshold}K - Next {threshold_next}K")
        # param_clusters = list_clusters(global_label_volume)

        # --- Step 1 : Re-expanding existing objects... ---
        mask = (volume_bt < threshold).astype(np.uint8)
        # list_clusters(global_label_volume)

        global_label_volume = spread_labels_fast(
            volume_bt=volume_bt,
            labeled_volume=global_label_volume,
            cloud_mask=mask,
            delta_spread=deltaBT_Spread
        )


        global_label_volume_tmp = global_label_volume.copy()
        # list_clusters(global_label_volume)

        # --- Step 2 : Detection of new conv seeds ---
        global_label_volume,nb_ConvSeeds = detect_objects(data_param,clusters,volume_bt,surface_area_2d,global_label_volume,threshold,area_threshold_km2,kernel,nb_ConvSeeds,labelMin )

        global_label_volume_tmp =  global_label_volume - global_label_volume_tmp

        # param_clusters = list_clusters(global_label_volume)

        # --- step 3 : Dilation of conv seeds up to intermediate boundaries ---
        print("  Expanding with next threshold...",np.shape(global_label_volume),np.nanmax(global_label_volume))
        mask = (volume_bt <= threshold_next).astype(np.uint8)
        global_label_volume = spread_labels_fast(
            volume_bt=volume_bt,
            labeled_volume=global_label_volume,
            cloud_mask=mask,
            delta_spread=deltaBT_Spread
        )

        # list_clusters(global_label_volume)

        #clusters = list_clusters(global_label_volume)


    mask = (volume_bt <= threshold_next).astype(np.uint8)
    global_label_volume = spread_labels_fast(
        volume_bt=volume_bt,
        labeled_volume=global_label_volume,
        cloud_mask=mask,
        delta_spread=10
    )

    print(f"Total clusters detected: {nCluster}")
    print(f"Total runtime: {time.time() - start_time:.2f} s")
    


    return global_label_volume,nb_ConvSeeds