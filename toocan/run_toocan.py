# =============================================================================
# File        : main.py
# Author      : Thomas Fiolleau
# Date        : 2025-07-04
# Description : Main driver for TOOCAN algorithm execution. Handles chunk-based
#               processing of satellite brightness temperature data to detect and
#               label organized deep convection systems.
#
# Usage       : python main.py
# Dependencies: Requires parameter files, navigation files, and a valid file list.
#
# Project     : TOOCAN - Tracking Organized Deep Convection
# License     : 
# =============================================================================
# pip install -e . gcc -O3 -fPIC -shared -o label.so label.c

# =============================================================================
# File        : run_toocan.py
# Author      : Thomas Fiolleau
# Project     : TOOCAN - Tracking Of Organized Convection using 3-dimensional segmentatioN
# Description : Main execution pipeline for the TOOCAN algorithm. Handles:
#                   1) Parameter loading
#                   2) Navigation cropping
#                   3) IR volume extraction
#                   4) Cloud cluster detection & spreading
#                   5) Tracking (INT/LC)
#                   6) FileTracking NetCDF writing
#
# Usage :
#       python run_toocan.py fileparam_TOOCAN.dat fileparam_GEO.dat ...
#
# Or python API :
#       from toocan.run_toocan import run_toocan
#       run_toocan(...)
#
# License     : CNRS / LEGOS / IPSL
# =============================================================================

#import os
#import numpy as np
#import pandas as pd
#import warnings
#from datetime import datetime, timedelta
#from xarray.conventions import SerializationWarning
#import sys 
#import ctypes
#
#import glob
#from ctypes import Structure, c_int, c_float, c_char, c_ulong
#
#from datetime import datetime
#import re
#import os
#import pandas as pd  # Only if you still use df, otherwise skip
#
## Suppress xarray warnings
#warnings.filterwarnings("ignore", category=SerializationWarning)
#
#
#sys.path.insert(0, os.path.expanduser("~/TOOCAN/pyTOOCAN/src"))
#
## === Local TOOCAN imports ===
#from toocan.utils.parse_param_file import parse_param_file
#from toocan.io.open_navigationFile import crop_navigation_file
#from toocan.io.open_IRdata import extract_volume
#from toocan.io.open_IRdata import read_irbt_subset
#from toocan.io.writer_toocanImage import save_labels_slot_by_slot
#from toocan.detection_spreading.detect_and_spread import detect_and_spread
#from toocan.struct.data_param import DataParam
#from toocan.struct.blob import make_Blob_class


import os
import sys
import numpy as np
import pandas as pd
import warnings
import ctypes
from datetime import datetime, timedelta
from xarray.conventions import SerializationWarning
import xarray as xr
#sys.path.insert(0, "/home/fiolleau/TOOCAN/pyTOOCAN")
# Silence certains warnings xarray
warnings.filterwarnings("ignore", category=SerializationWarning)

# =============================================================================
# TOOCAN internal imports
# =============================================================================

# ---- Utils ----
from toocan.lib.utils.parse_param_file import parse_param_file
from toocan.lib.utils.memory_utils import (
    mem_used_mb,
    mem_available_mb,
    estimate_array_mb,
    warn_if_not_enough_memory,
)

# ---- IO ----
from toocan.lib.io.open_navigationFile import crop_navigation_file
from toocan.lib.io.open_IRdata import extract_volume
from toocan.lib.io.writer_toocanImage import save_labels_slot_by_slot
from toocan.lib.io.file_listing import build_ir_filelist,fast_extract_datetime
from toocan.lib.io.writer_navigation import save_navigation_grid
from toocan.lib.io.lauch_resumption import launch_resumption
from toocan.lib.io.compute_VZA_correction import extract_VZARegcoefs, compute_VZA_correction

# ---- Detection / Spreading ----
from toocan.lib.detection_spreading.detect_and_spread import detect_and_spread
from toocan.lib.detection_spreading.cluster_tools import clean_clusters, create_clusters_reprise

# ---- Structures ----
from toocan.lib.struct.data_param import DataParam
from toocan.lib.struct.blob import make_Blob_class

# =============================================================================
# End imports
# =============================================================================



toocan_param_path   = sys.argv[1]
geo_param_path      = sys.argv[2]

yearBEGIN  = int(sys.argv[3])
monthBEGIN = int(sys.argv[4])
dayBEGIN   = int(sys.argv[5])
hourBEGIN  = int(sys.argv[6])
minBEGIN   = int(sys.argv[7])

yearEND    = int(sys.argv[8])
monthEND   = int(sys.argv[9])
dayEND     = int(sys.argv[10])
hourEND    = int(sys.argv[11])
minEND     = int(sys.argv[12])

lonmin     = float(sys.argv[13])
lonmax     = float(sys.argv[14])
latmin     = float(sys.argv[15])
latmax     = float(sys.argv[16])




# === Load TOOCAN Parameters ===
params_TOOCAN = parse_param_file(toocan_param_path)


start_str = (
    f"{yearBEGIN:04d}-"
    f"{monthBEGIN:02d}-"
    f"{dayBEGIN:02d} "
    f"{hourBEGIN:02d}:"
    f"{minBEGIN:02d}"
)
start_time = datetime.strptime(start_str, "%Y-%m-%d %H:%M")

end_str = (
    f"{yearEND:04d}-"
    f"{monthEND:02d}-"
    f"{dayEND:02d} "
    f"{hourEND:02d}:"
    f"{minEND:02d}"
)
end_time = datetime.strptime(end_str, "%Y-%m-%d %H:%M")

print (start_time,end_time)

pathout_TOOCAN = params_TOOCAN.get("pathout_TOOCAN")

file_list = params_TOOCAN.get("file_list")

# Thresholds and constants
overlap_window_size = int(params_TOOCAN.get("overlap_window_size"))
firstlabel = int(params_TOOCAN.get("firstlabel"))

##############
data_param = DataParam()

data_param.reprise = params_TOOCAN.get("reprise")   # if bug encountered
data_param.date_reprise = params_TOOCAN.get("date_reprise")
data_param.hour_reprise = params_TOOCAN.get("hour_reprise")

# --- Time window ---
data_param.yearBegin = int(yearBEGIN)
data_param.monthBegin = int(monthBEGIN)
data_param.dayBegin = int(dayBEGIN)
data_param.hourBegin = int(hourBEGIN)
data_param.minBegin = int(minBEGIN)

data_param.yearEnd = int(yearEND)
data_param.monthEnd = int(monthEND)
data_param.dayEnd = int(dayEND)
data_param.hourEnd = int(hourEND)
data_param.minEnd = int(minEND)

# --- Domain ---
data_param.latmin = float(latmin)
data_param.latmax = float(latmax)
data_param.lonmin = float(lonmin)
data_param.lonmax = float(lonmax)


# --- Metadata ---
data_param.version = bytes(str(params_TOOCAN.get("version", "")), "utf-8")
data_param.path_out = bytes(str(params_TOOCAN.get("pathout_TOOCAN", "")), "utf-8")
data_param.path_fileIN = bytes(str(params_TOOCAN.get("file_list", "")), "utf-8")

# --- BT Thresholds---
data_param.minBT = int(params_TOOCAN["minBT_threshold"])
data_param.maxBT = int(params_TOOCAN["maxBT_threshold"])
data_param.stepBT = int(params_TOOCAN["stepBT_threshold"])

# --- Parametres ---
data_param.deltaDetect = float(params_TOOCAN.get("deltaDetect", 1.0))
data_param.deltaSpread = float(params_TOOCAN.get("deltaSpread", 1.0))

data_param.timin = int(params_TOOCAN["minAreaSeed"])
data_param.lifemin = int(params_TOOCAN["minLifetime"])
data_param.labelFirstMCS = int(params_TOOCAN["firstlabel"])

data_param.ZSIZE = int(params_TOOCAN["VolumeImage"])
data_param.overlap_window_size = int(params_TOOCAN["overlap_window_size"])
data_param.nbMaxCluster = int(params_TOOCAN["nbMaxCluster"])
data_param.maxMising = int(params_TOOCAN["max_missing"])

#############


# clusters class creation
Blob = make_Blob_class()
ClustersArray = Blob * data_param.nbMaxCluster
clusters = ClustersArray()      
# Allocate buffers for EACH blob
for i in range(data_param.nbMaxCluster):
    clusters[i].seed_area_perFrame = (ctypes.c_int * data_param.ZSIZE)()
    clusters[i].labelVoisin       = (ctypes.c_int * 1000)()

# === Load GEO parameters ===
params_GEO = parse_param_file(geo_param_path)
temporalresolution = int(params_GEO.get("temporalresolution"))
file_navigation = params_GEO.get("file_navigation")
file_list_path = params_GEO.get("file_listGEO")
vza_path = params_GEO.get("path_vza")
model_name    = params_GEO.get("GEOplatform")   # ex: "ARPEGENH"
variable_name = params_GEO.get("variable")      # ex: "BT"

# === Crop Navigation File ===
cropped_ds = crop_navigation_file(
    file_navigation, model_name,
    lonmin=lonmin, lonmax=lonmax,
    latmin=latmin, latmax=latmax,
)


# Load regression coefficients alpha and beta for rapid scan
df_dict = {}
if model_name == 'MSGrss':
    regression_path = params_GEO.get("path_coeffs")
    df_file = pd.read_csv(regression_path)
    df_dict = df_file.set_index('sat')[['alpha', 'beta', 'nuc']].to_dict('index')


# Computing VZA correction
# Get VZA reg coeffs
coefVZA_a, coefVZA_b, coefVZA_c, VZAmax, BTmax = extract_VZARegcoefs(vza_path)
mat_coefVZA_ax, mat_coefVZA_bx, mat_coefVZA_cx = compute_VZA_correction(coefVZA_a, coefVZA_b, coefVZA_c, cropped_ds)
VZA_coeffs = [mat_coefVZA_ax, mat_coefVZA_bx, mat_coefVZA_cx, VZAmax, BTmax]

surface_area_2d = cropped_ds['mat_surfacePix'].values
lon_array_2d = cropped_ds['mat_longitude'].values
lat_array_2d = cropped_ds['mat_latitude'].values
lat_size, lon_size = surface_area_2d.shape


nav_file = save_navigation_grid(cropped_ds, params_TOOCAN, params_GEO, "navigation_grid.nc")



# After crop navigation :
data_param.XSIZE = int(lon_array_2d.shape[1])
data_param.YSIZE = int(lon_array_2d.shape[0])

## Case 1: Use static CSV
#if "file_list" in params_GEO and params_GEO["file_list"] and params_GEO["file_list"] != "AUTO":
#    print("Reading file list from CSV...")
#    df = pd.read_csv(params_GEO["file_list"], delim_whitespace=True, engine='python')
#    df['full_path'] = df['path_ir'].str.rstrip('/') + '/' + df['file_ir']
#    df['datetime'] = pd.to_datetime(df[['year', 'month', 'day', 'hour', 'minute']])
#
## Case 2: Use dynamic directory scan
#else:

# read and filter file list netCDF
file_list = xr.open_dataset(file_list_path)

print("Building file list from path_ir...")

dates = pd.to_datetime({
    "year": file_list["year"].values,
    "month": file_list["month"].values,
    "day": file_list["day"].values,
    "hour": file_list["hour"].values,
    "minute": file_list["minute"].values,
})

df = pd.DataFrame({
    "full_path": file_list['path_ir'] + file_list['file_ir'],
    "datetime": dates
})
df = df.sort_values("datetime").reset_index(drop=True)

# === Chunking Logic ===
start_global = pd.to_datetime(start_time)
end_global = pd.to_datetime(end_time)

chunk_length = timedelta(minutes=temporalresolution * data_param.ZSIZE)
timesteps_per_chunk = int((chunk_length.total_seconds() // 60) / temporalresolution)
print('timesteps_per_chunk: ',chunk_length.total_seconds(),temporalresolution, int((chunk_length.total_seconds() // 60) / temporalresolution))


# ===== MEMORY CHECK BEFORE CREATING HUGE ARRAYS =====
# Global label volume shape = (T, Y, X)
T = timesteps_per_chunk
Y = lat_size
X = lon_size

mem_BT   = estimate_array_mb((T, Y, X), dtype=np.float32)   # volume_BT
mem_LAB  = estimate_array_mb((T, Y, X), dtype=np.int32)     # global_label_volume
mem_TEMP = estimate_array_mb((T, Y, X), dtype=np.uint8)     # masks, temporary arrays
mem_lon = estimate_array_mb(( Y, X), dtype=np.float32)     # masks, temporary arrays
mem_lat = estimate_array_mb(( Y, X), dtype=np.float32)     # masks, temporary arrays
mem_area = estimate_array_mb(( Y, X), dtype=np.float32)     # masks, temporary arrays

estimated_total = mem_BT + mem_LAB + mem_TEMP +mem_lon +mem_lat+mem_area

if not warn_if_not_enough_memory(estimated_total):
    sys.exit(1)
# =====================================================


global_label_volume = np.zeros((timesteps_per_chunk, lat_size, lon_size), dtype=np.int32)

current_start = start_global
labelMin      = firstlabel
nb_ConvSeeds  = firstlabel

# If resuming a run, retrieve the last overlap_window_size TOOCAN outputs and update clusters
if data_param.reprise == 1:
    global_label_volume = launch_resumption(current_start, data_param, df, temporalresolution, global_label_volume, nomenclature="ToocanCloudMask_")
    labels_present = np.unique(global_label_volume[global_label_volume > 0])
    nb_ConvSeeds = np.max(labels_present)
    labelMin = np.min(labels_present)
    clusters = create_clusters_reprise(
        data_param,
        clusters,
        labels_present=labels_present,
        labelMin=labelMin,
        nbMax=data_param.nbMaxCluster
    )
    data_param.reprise = 0
    start_global = start_global - 2 * timedelta(minutes=temporalresolution * overlap_window_size)


while current_start < end_global:

    # ---------------------------------------------------------
    # 1) Define temporal window for this chunk
    # ---------------------------------------------------------
    current_end = current_start + chunk_length
    actual_end  = min(current_end, end_global)

    pre_overlap = max(
        current_start - timedelta(minutes=temporalresolution * overlap_window_size),
        start_global
    )
        

    post_overlap = min(
        pre_overlap + chunk_length - timedelta(minutes=temporalresolution),
        end_global
    )

    print(f"\nProcessing window: {pre_overlap} → {post_overlap}")

    # ---------------------------------------------------------
    # 2) Handle label continuity through overlap region
    # ---------------------------------------------------------
    if current_start > start_global:

        # Keep labels from previous overlap
        prev_overlap_labels = global_label_volume[-overlap_window_size:].copy()
        global_label_volume[:] = 0
        global_label_volume[:overlap_window_size] = prev_overlap_labels

        # Active labels inside overlap
        labels_present = np.unique(global_label_volume[global_label_volume > 0])

        if labels_present.size > 0:

            labelMin = np.min(labels_present)

            # Reorganize and compact the cluster table
            clusters = clean_clusters(
                data_param,
                clusters,
                labels_present=labels_present,
                labelMin=labelMin,
                nbMax=data_param.nbMaxCluster
            )

            labelMin = labelMin - 1
            print(f"Minimum label shifted to {labelMin}")


    else:
        # First chunk: no continuity to preserve
        global_label_volume = np.zeros((timesteps_per_chunk, lat_size, lon_size), dtype=np.int32)       # redefine size of volume after cut

        # Allocation des buffers pour CHAQUE blob
        clusters = ClustersArray()
        for i in range(data_param.nbMaxCluster):
            clusters[i].seed_area_perFrame = (ctypes.c_int * data_param.ZSIZE)()
            clusters[i].labelVoisin       = (ctypes.c_int * 1000)()

    # ---------------------------------------------------------
    # 3) Extract 3D IR brightness temperature volume
    # ---------------------------------------------------------
    volume_BT, times, lat, lon, flag_cut, next_date = extract_volume(
        df, pre_overlap, post_overlap, file_list, vza_path, data_param.maxMising, df_dict, VZA_coeffs, nav=cropped_ds,
        model_name=model_name
    )

    T_actual = volume_BT.shape[0]
    T_full   = data_param.ZSIZE
    missing = 0

    print("T_actual", T_actual)
    print("T_full", T_full)
    print("shape global_label_volume", global_label_volume.shape)

    # ---------------------------------------------------------
    # 4) Determine chunk type (normal / last / missing data)
    # ---------------------------------------------------------

    if T_actual == T_full:
        # Normal chunk: nothing to adjust
        pass

    else:
        # volume shorter than expected: could be last window or missing images
        if post_overlap >= end_global:
            # Legitimate last chunk (end of requested period)
            print(f"[INFO] Last chunk is incomplete: {T_actual}/{T_full}")
            global_label_volume = global_label_volume[:T_actual, :, :]

        else:
            # Missing IR images inside the dataset (unexpected)
            missing = T_full - T_actual
            if missing < 0:
                exit()
            print("⚠️ WARNING: Incomplete chunk BEFORE reaching end of period")
            print(f"    → {missing} missing IR image(s)")
            print(f"    pre_overlap  = {pre_overlap}")
            print(f"    post_overlap = {post_overlap}")
            print(f"    end_global   = {end_global}")
            global_label_volume = global_label_volume[:T_actual, :, :]

            # strict mode — stop processing
            # sys.exit("ERROR: Missing IR slot inside a non-terminal chunk.")

    # Final safety check
    # if volume_BT.shape != global_label_volume.shape:
    #     print("❌ ERROR: Shape mismatch after dimensional adjustment.")
    #     print("volume_BT:          ", volume_BT.shape)
    #     print("global_label_volume:", global_label_volume.shape)
    #     sys.exit(1)

    # ---------------------------------------------------------
    # 5) Run detection + spatial-temporal spreading
    # ---------------------------------------------------------
    
    if missing < data_param.ZSIZE - data_param.lifemin:
        global_label_volume,nb_ConvSeeds = detect_and_spread(
            data_param,
            clusters,
            volume_BT,
            surface_area_2d,
            lon_array_2d,
            lat_array_2d,
            global_label_volume,
            params_TOOCAN,
            nb_ConvSeeds,
            labelMin
        )


        # ---------------------------------------------------------
        # 6) Save output cloud labels for this chunk
        # ---------------------------------------------------------
        save_labels_slot_by_slot(
            global_label_volume, times, lat, lon,
            nomenclature="ToocanCloudMask_",
            var_name="DCS_number",
            params_TOOCAN=params_TOOCAN, params_GEO=params_GEO
        )


    
    # ---------------------------------------------------------
    # 7) Move to next temporal window
    # ---------------------------------------------------------
    
    if flag_cut:
        current_start = pd.Timestamp(next_date, unit="ns")
        labelMin = nb_ConvSeeds
        start_global = current_start
        overlap_window_size = 0
        data_param.overlap_window_size = 0
    else:
        current_start = post_overlap + timedelta(minutes=temporalresolution)
        overlap_window_size = int(params_TOOCAN.get("overlap_window_size"))
        data_param.overlap_window_size = int(params_TOOCAN["overlap_window_size"])
    data_param.ZSIZE = int(params_TOOCAN["VolumeImage"])




