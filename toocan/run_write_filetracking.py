import sys
from datetime import datetime
import pandas as pd
import os
sys.path.insert(0, os.path.expanduser("~/TOOCAN/pyTOOCAN/toocan/"))

from toocan.lib.utils.parse_param_file import parse_param_file
from toocan.lib.io.open_navigationFile import crop_navigation_file
from toocan.lib.io.file_listing import (
    build_ir_filelist,
    build_toocan_mask_filelist,
    merge_ir_and_toocan
)
from toocan.lib.postprocessing.tracking.tracking_DCS import track_month_from_df
from toocan.lib.postprocessing.tracking.int_builder import compute_INT_variables
from toocan.lib.postprocessing.tracking.lc_builder import build_LC_variables
from toocan.lib.postprocessing.tracking.nc_writer import write_tracking_nc

# ---- Structures ----
from toocan.lib.struct.data_param import DataParam




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


# --- Fenêtre temporelle ---
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

# --- Domaine ---
data_param.latmin = float(latmin)
data_param.latmax = float(latmax)
data_param.lonmin = float(lonmin)
data_param.lonmax = float(lonmax)


# --- Métadonnées ---
data_param.version = bytes(str(params_TOOCAN.get("version", "")), "utf-8")
data_param.path_out = bytes(str(params_TOOCAN.get("pathout_TOOCAN", "")), "utf-8")
data_param.path_fileIN = bytes(str(params_TOOCAN.get("file_list", "")), "utf-8")

# --- Seuils BT ---
data_param.minBT = int(params_TOOCAN["minBT_threshold"])
data_param.maxBT = int(params_TOOCAN["maxBT_threshold"])
data_param.stepBT = int(params_TOOCAN["stepBT_threshold"])

# --- Paramètres ---
data_param.deltaDetect = float(params_TOOCAN.get("deltaDetect", 1.0))
data_param.deltaSpread = float(params_TOOCAN.get("deltaSpread", 1.0))

data_param.timin = int(params_TOOCAN["minAreaSeed"])
data_param.lifemin = int(params_TOOCAN["minLifetime"])
data_param.labelFirstMCS = int(params_TOOCAN["firstlabel"])

data_param.ZSIZE = int(params_TOOCAN["VolumeImage"])
data_param.overlap_window_size = int(params_TOOCAN["overlap_window_size"])
data_param.nbMaxCluster = int(params_TOOCAN["nbMaxCluster"])
path_masks = data_param.path_out.decode()+"/toocan_2.08/"

#############

# === Load GEO parameters ===
params_GEO = parse_param_file(geo_param_path)
temporalresolution = int(params_GEO.get("temporalresolution"))
file_navigation = params_GEO.get("file_navigation")
model_name    = params_GEO.get("GEOplatform")   # ex: "ARPEGENH"
variable_name = params_GEO.get("variable")      # ex: "BT"


# === Crop Navigation File ===
cropped_ds = crop_navigation_file(
    file_navigation,
    lonmin=lonmin, lonmax=lonmax,
    latmin=latmin, latmax=latmax
)
surface_area_2d = cropped_ds['mat_surfacePix'].values
lon_array_2d = cropped_ds['mat_longitude'].values
lat_array_2d = cropped_ds['mat_latitude'].values
lat_size, lon_size = surface_area_2d.shape

fileIR_time = build_ir_filelist(params_GEO["path_ir"], start_time, end_time)
df_ir = pd.DataFrame(fileIR_time, columns=["ir_path", "datetime"])
maskTOOCAN_time = build_toocan_mask_filelist(path_masks, start_time, end_time)
df_mask = pd.DataFrame(maskTOOCAN_time, columns=["mask_path", "datetime"])
df_all = merge_ir_and_toocan(df_ir, df_mask)

print("lat_array_2d shape =", lon_array_2d.shape)
print("lon_array_2d shape =", lat_array_2d.shape)
print("Surface grid shape =", surface_area_2d.shape)

# 1) tracking complet
tracks = track_month_from_df(df_all, cropped_ds,surface_area_2d,lon_array_2d,lat_array_2d, lonmin, lonmax, latmin, latmax,model_name=model_name)

INT = compute_INT_variables(tracks,dt_minutes=temporalresolution)
global_times = sorted({t for tr in tracks.values() for t in tr["times"]})
LC = build_LC_variables(tracks, global_times)
# 5) écriture NetCDF
outfile = "TOOCAN-AFRICA-20120801-20120831.nc"  # à construire depuis start/end
write_tracking_nc(outfile, INT, LC, global_times,
                  region="AFRICA",
                  platform=params_GEO["GEOplatform"],
                  version=params_TOOCAN["version"])
