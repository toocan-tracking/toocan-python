import numpy as np
import xarray as xr
from datetime import datetime
from scipy.ndimage import center_of_mass
import os
import numpy as np
import glob

#def compute_dcs_statistics(global_label_volume, volume_bt, surface_area_2d, time_array, lat_array, lon_array):
#    dcs_ids = np.unique(global_label_volume)
#    dcs_ids = dcs_ids[dcs_ids > 0]  # Exclude background
#
#    time_seconds = np.array([(t - datetime(1970, 1, 1)).total_seconds() for t in time_array])
#    dt_seconds = (time_array[1] - time_array[0]).total_seconds()
#    num_times = len(time_array)
#
#    dcs_stats_LC = []
#
#    for label_id in dcs_ids:
#        mask = (global_label_volume == label_id)
#        times_present = np.where(np.any(mask, axis=(1, 2)))[0]
#        if len(times_present) == 0:
#            continue
#
#        t0, t1 = times_present[0], times_present[-1]
#        lifetime = (t1 - t0 + 1) * dt_seconds / 3600.0  # in hours
#
#        # Position (center of mass)
#        y0, x0 = center_of_mass(mask[t0])
#        y1, x1 = center_of_mass(mask[t1])
#        lat0 = float(lat_array[int(round(y0))])
#        lon0 = float(lon_array[int(round(x0))])
#        lat1 = float(lat_array[int(round(y1))])
#        lon1 = float(lon_array[int(round(x1))])
#
#        # Distance and velocity
#        dist = haversine_km(lat0, lon0, lat1, lon1)
#        velocity = dist * 1000.0 / ((t1 - t0 + 1) * dt_seconds)  # m/s
#
#        # Tb stats
#        tb_values = volume_bt[mask]
#        tbmin = tb_values.min()
#
#        # Surface area stats at 235K
#        max_area = 0
#        cum_area = 0
#
#        lc_lat = np.full(num_times, -999.0, dtype=np.float32)
#        lc_lon = np.full(num_times, -999.0, dtype=np.float32)
#        lc_tbmin = np.full(num_times, -999.0, dtype=np.float32)
#        lc_surf235 = np.full(num_times, -999.0, dtype=np.float32)
#        lc_time = np.full(num_times, -999, dtype=np.int32)
#
#        for t in times_present:
#            cold_mask = (volume_bt[t] <= 235) & (global_label_volume[t] == label_id)
#            area = surface_area_2d[cold_mask].sum()
#            cum_area += area
#            max_area = max(max_area, area)
#
#            lc_time[t] = int(time_seconds[t])
#            cy, cx = center_of_mass(mask[t])
#            lc_lat[t] = float(lat_array[int(round(cy))])
#            lc_lon[t] = float(lon_array[int(round(cx))])
#            lc_tbmin[t] = float(volume_bt[t][mask[t]].min())
#            lc_surf235[t] = float(area)
#
#        dcs_stats_LC.append({
#            "label": int(label_id),
#            "t0": int(time_seconds[t0]),
#            "t1": int(time_seconds[t1]),
#            "duration": lifetime,
#            "lat_init": lat0,
#            "lon_init": lon0,
#            "lat_end": lat1,
#            "lon_end": lon1,
#            "distance": dist,
#            "velocity": velocity,
#            "tbmin": float(tbmin),
#            "surfmax_235K": float(max_area),
#            "surfcum_235K": float(cum_area),
#            "LC_lat": lc_lat,
#            "LC_lon": lc_lon,
#            "LC_tbmin": lc_tbmin,
#            "LC_surfkm2_235K": lc_surf235,
#            "LC_UTC_time": lc_time
#        })
#
#    return dcs_stats_LC
import pandas as pd
import numpy as np
from datetime import datetime
import xarray as xr
import numpy as np
import pandas as pd
import warnings

def read_irbt_subset(filepath, lonmin, lonmax, latmin, latmax):
    """
    Reads a GEO L1C-MSG NetCDF file and extracts a spatial subset of Harmonized IRBT data.

    Parameters:
        filepath (str): Path to the NetCDF file.
        lonmin (float): Minimum longitude for subset.
        lonmax (float): Maximum longitude for subset.
        latmin (float): Minimum latitude for subset.
        latmax (float): Maximum latitude for subset.

    Returns:
        dict: Dictionary containing:
            - 'irbt': 2D array of subsetted brightness temperature (in K)
            - 'lat': 1D array of latitudes (subset)
            - 'lon': 1D array of longitudes (subset)
            - 'time': observation time (as pd.Timestamp)
    """
    # Open dataset
    ds = xr.open_dataset(filepath)

    # Subset the domain
#    ds_sub = ds.sel(
#        latitude=slice(latmax, latmin),  # latitude usually decreasing
#        longitude=slice(lonmin, lonmax)
#    )
    if ds['latitude'][0] > ds['latitude'][-1]:  # decreasing
        ds_sub = ds.sel(
           latitude=slice(latmax, latmin),
           longitude=slice(lonmin, lonmax)
        )
    else:  # increasing
        ds_sub = ds.sel(
           latitude=slice(latmin, latmax),
           longitude=slice(lonmin, lonmax)
        )
    # Extract and scale IRBT
    raw_irbt = ds_sub['Harmonized_irBT'].isel(time=0).astype(np.float32)
    scale = ds_sub['Harmonized_irBT'].attrs.get('scale_factor', 1.0)
    fill_value = ds_sub['Harmonized_irBT'].attrs.get('_FillValue', -99800)

    irbt_scaled = np.where(raw_irbt == fill_value, np.nan, raw_irbt * scale)

    # Extract lat/lon grid
    lat = ds_sub['latitude'].values
    lon = ds_sub['longitude'].values
    
    # Time
    time = pd.to_datetime(ds['time'].values[0])

    return {
        'irbt': irbt_scaled,
        'lat': lat,
        'lon': lon,
        'time': time
    }


import xarray as xr
import pandas as pd

def read_toocan_label(filepath):
    """
    Reads a TOOCAN label NetCDF file and extracts the label mask and time.

    Parameters:
        filepath (str): Path to the TOOCAN NetCDF file.

    Returns:
        dict: Dictionary containing:
            - 'label': 2D array of labeled objects
            - 'lat': 1D or 2D latitude grid (if available)
            - 'lon': 1D or 2D longitude grid (if available)
            - 'time': timestamp (as pd.Timestamp)
    """
    ds = xr.open_dataset(filepath)

    label = ds['ToocanCloudMask_'].isel(time=0).values.astype(np.int32)

    # Handle optional lat/lon
    lat = ds['latitude'].values if 'latitude' in ds else None
    lon = ds['longitude'].values if 'longitude' in ds else None

    #print(ds['latitude'][0],ds['latitude'][-1])
    #print(ds['longitude'][0],ds['longitude'][-1])

    # Time
    if 'time' in ds:
        time = pd.to_datetime(ds['time'].values[0])
    elif 'UTC_time' in ds.attrs:
        time = pd.to_datetime(ds.attrs['UTC_time'])
    else:
        raise ValueError("Time information not found in TOOCAN file.")

    return {
        'label': label,
        'lat': lat,
        'lon': lon,
        'time': time
    }

def create_dcs_netcdf_from_dicts(dcs_stats_LC, output_path):
    num_dcs = len(dcs_stats_LC)
    num_times = len(dcs_stats_LC[0]['LC_lat'])

    def stack(name):
        return np.stack([d[name] for d in dcs_stats_LC], axis=0)

    ds = xr.Dataset(
        {
            "INT_DCSnumber": ("DCS", [d["label"] for d in dcs_stats_LC]),
            "INT_duration": ("DCS", [d["duration"] for d in dcs_stats_LC]),
            "INT_UTC_timeInit": ("DCS", [d["t0"] for d in dcs_stats_LC]),
            "INT_UTC_timeEnd": ("DCS", [d["t1"] for d in dcs_stats_LC]),
            "INT_latInit": ("DCS", [d["lat_init"] for d in dcs_stats_LC]),
            "INT_lonInit": ("DCS", [d["lon_init"] for d in dcs_stats_LC]),
            "INT_latEnd": ("DCS", [d["lat_end"] for d in dcs_stats_LC]),
            "INT_lonEnd": ("DCS", [d["lon_end"] for d in dcs_stats_LC]),
            "INT_distance": ("DCS", [d["distance"] for d in dcs_stats_LC]),
            "INT_velocityAvg": ("DCS", [d["velocity"] for d in dcs_stats_LC]),
            "INT_tbmin": ("DCS", [d["tbmin"] for d in dcs_stats_LC]),
            "INT_surfmaxkm2_235K": ("DCS", [d["surfmax_235K"] for d in dcs_stats_LC]),
            "INT_surfcumkm2_235K": ("DCS", [d["surfcum_235K"] for d in dcs_stats_LC]),

            "LC_lat": (("DCS", "time"), stack("LC_lat")),
            "LC_lon": (("DCS", "time"), stack("LC_lon")),
            "LC_tbmin": (("DCS", "time"), stack("LC_tbmin")),
            "LC_surfkm2_235K": (("DCS", "time"), stack("LC_surfkm2_235K")),
            "LC_UTC_time": (("DCS", "time"), stack("LC_UTC_time"))
        },
        coords={
            "DCS": np.arange(1, num_dcs + 1, dtype=np.int32),
            "time": np.arange(num_times, dtype=np.int32)
        },
        attrs={
            "title": "TOOCAN-derived DCS tracking dataset with life cycle variables",
            "institution": "CNRS/LEGOS/IPSL",
            "conventions": "CF-1.6",
            "creator_name": "Thomas Fiolleau",
            "contact": "thomas.fiolleau@cnrs.fr",
        }
    )

    ds.to_netcdf(output_path, format="NETCDF4")
    print(f"NetCDF written to: {output_path}")


####################################################################################







import numpy as np
from datetime import datetime
from scipy.ndimage import center_of_mass
from collections import defaultdict

def haversine_km(lat1, lon1, lat2, lon2):
    R = 6371.0
    phi1, phi2 = np.radians(lat1), np.radians(lat2)
    dphi = phi2 - phi1
    dlambda = np.radians(lon2 - lon1)
    a = np.sin(dphi / 2.0) ** 2 + np.cos(phi1) * np.cos(phi2) * np.sin(dlambda / 2.0) ** 2
    return R * 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))

class DCSLifeCycleTracker:
    def __init__(self, surface_area_2d, lat_array, lon_array):
        self.surface_area_2d = surface_area_2d
        self.lat = lat_array
        self.lon = lon_array
        self.track_dict = defaultdict(lambda: {
            "times": [],
            "lat": [],
            "lon": [],
            "tbmin": [],
            "surf_235": [],
            "cum_area": 0.0,
            "max_area": 0.0,
            "t0": None,
            "t1": None
        })

    def update(self, label_array, irbt_array, timestamp):
        """
        Update the tracker with a single timestep's label and IR image.
        """
        for label_id in np.unique(label_array):
            if label_id == 0:
                continue

            mask = (label_array == label_id)
            if not np.any(mask):
                continue

            cy, cx = center_of_mass(mask)
            #lat = float(self.lat[int(round(cy))])
            #lon = float(self.lon[int(round(cx))])
            if self.lat.ndim == 2:
                lat = float(self.lat[int(round(cy)), int(round(cx))])
                lon = float(self.lon[int(round(cy)), int(round(cx))])
            else:
                lat = float(self.lat[int(round(cy))])
                lon = float(self.lon[int(round(cx))])

            tbmin = float(irbt_array[mask].min())
            area_235 = float(self.surface_area_2d[mask & (irbt_array <= 235)].sum())
            print(tbmin,area_235)
            record = self.track_dict[label_id]
            record["times"].append(timestamp)
            record["lat"].append(lat)
            record["lon"].append(lon)
            record["tbmin"].append(tbmin)
            record["surf_235"].append(area_235)
            record["cum_area"] += area_235
            record["max_area"] = max(record["max_area"], area_235)

            if record["t0"] is None:
                record["t0"] = timestamp
            record["t1"] = timestamp

    def finalize(self):
        """
        Convert track_dict into a list of DCS summary dicts.
        """
        stats = []
        for label_id, rec in self.track_dict.items():
            if len(rec["times"]) < 2:
                continue

            t0 = rec["t0"]
            t1 = rec["t1"]
            dt = (t1 - t0).total_seconds()
            duration_hr = dt / 3600.0

            lat0, lon0 = rec["lat"][0], rec["lon"][0]
            lat1, lon1 = rec["lat"][-1], rec["lon"][-1]
            dist = haversine_km(lat0, lon0, lat1, lon1)
            velocity = dist * 1000.0 / dt if dt > 0 else 0.0

            time_seconds = [int((t - datetime(1970, 1, 1)).total_seconds()) for t in rec["times"]]
            n_time = len(time_seconds)
            lc_pad = lambda arr: np.pad(arr, (0, self.max_length - n_time), constant_values=-999)

            stats.append({
                "label": label_id,
                "t0": time_seconds[0],
                "t1": time_seconds[-1],
                "duration": duration_hr,
                "lat_init": lat0,
                "lon_init": lon0,
                "lat_end": lat1,
                "lon_end": lon1,
                "distance": dist,
                "velocity": velocity,
                "tbmin": float(np.min(rec["tbmin"])),
                "surfmax_235K": rec["max_area"],
                "surfcum_235K": rec["cum_area"],
                "LC_lat": np.array(rec["lat"], dtype=np.float32),
                "LC_lon": np.array(rec["lon"], dtype=np.float32),
                "LC_tbmin": np.array(rec["tbmin"], dtype=np.float32),
                "LC_surfkm2_235K": np.array(rec["surf_235"], dtype=np.float32),
                "LC_UTC_time": np.array(time_seconds, dtype=np.int32)
            })

        return stats

import xarray as xr
import numpy as np
import pandas as pd
from datetime import timedelta

def crop_navigation_file(file_path, lonmin=0, lonmax=20, latmin=-5, latmax=15, save_path=None):
    """
    Robustly reads and crops a NetCDF navigation file to a lat/lon bounding box.
    Supports multiple variable naming conventions and both 1D or 2D lat/lon.
    """
    ds = xr.open_dataset(file_path)

    # Candidate names for lat/lon
    lat_candidates = ['latitude', 'lat', 'mat_latitude', 'matlat']
    lon_candidates = ['longitude', 'lon', 'mat_longitude', 'matlon']

    def find_coord_var(ds, candidates, axis='lat'):
        # Try name-based search
        for name in candidates:
            if name in ds:
                return ds[name]
        # Try attribute-based search
        for var in ds.variables:
            da = ds[var]
            if hasattr(da, 'standard_name') and (
                (axis == 'lat' and 'latitude' in da.standard_name.lower()) or
                (axis == 'lon' and 'longitude' in da.standard_name.lower())
            ):
                return da
            if hasattr(da, 'units'):
                if axis == 'lat' and 'degrees_north' in da.units:
                    return da
                if axis == 'lon' and 'degrees_east' in da.units:
                    return da
        return None

    lat_var = find_coord_var(ds, lat_candidates, axis='lat')
    lon_var = find_coord_var(ds, lon_candidates, axis='lon')

    if lat_var is None or lon_var is None:
        raise ValueError("No recognizable latitude/longitude variables found.")

    # Determine whether coordinates are 1D or 2D
    if lat_var.ndim == 1 and lon_var.ndim == 1:
        # Use sel on axis-based coordinates
        lat_mask = (lat_var >= latmin) & (lat_var <= latmax)
        lon_mask = (lon_var >= lonmin) & (lon_var <= lonmax)

        cropped_ds = ds.sel(
            {lat_var.dims[0]: lat_var[lat_mask],
             lon_var.dims[0]: lon_var[lon_mask]}
        )

    elif lat_var.ndim == 2 and lon_var.ndim == 2:
        # Use 2D mask for spatial filtering
        mask = ((lat_var >= latmin) & (lat_var <= latmax) &
                (lon_var >= lonmin) & (lon_var <= lonmax))
        cropped_ds = ds.where(mask, drop=True)

    else:
        raise ValueError("Unsupported coordinate dimensionality.")

    if save_path:
        cropped_ds.to_netcdf(save_path)
        print(f"Cropped file saved to: {save_path}")

    return cropped_ds

import re
from collections import defaultdict

def parse_param_file(filepath):
    params = defaultdict(list)

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()

            # Skip blank lines and comments
            if not line or line.startswith(';') or '=' not in line:
                continue

            # Parse key and value
            key, value = map(str.strip, line.split('=', 1))

            # Try converting value to int, then float, then leave as string
            try:
                value = int(value)
            except ValueError:
                try:
                    value = float(value)
                except ValueError:
                    value = value.strip('"').strip("'")  # remove quotes if any

            params[key].append(value)

    # Simplify lists that only have one item
    clean_params = {k: v[0] if len(v) == 1 else v for k, v in params.items()}
    return clean_params


def get_toocan_file(base_path, dt):
    """
    Constructs the full path to a TOOCAN file based on datetime.

    Parameters:
        base_path (str): Base directory, e.g. /bdd/MT_WORKSPACE/MCS/TOOCAN/pyTOOCAN/toocan_2.08/AFRICA/
        dt (datetime): Desired timestamp

    Returns:
        str: Full path to the TOOCAN file
    """
    ymd = dt.strftime("%Y%m%d")
    hm = dt.strftime("%H%M")
    subdir = dt.strftime("%Y/%Y_%m_%d/")
    filename = f"ToocanCloudMask_{ymd}_{hm}.nc"
    full_path = os.path.join(base_path, subdir, filename)
    return full_path

def get_GEO_IR_file(base_path, dt):
    """
    Finds the GEO IR file for a given datetime by globbing the directory.

    Parameters:
        base_path (str): e.g., "/bdd/GEOgrid_coldcloud/MSG+0000/"
        dt (datetime): Desired timestamp

    Returns:
        str or None: Full path to the matching IR file, or None if not found.
    """
    subdir = dt.strftime("%Y/%Y_%m_%d/")
    pattern_time = dt.strftime("%Y-%m-%dT%H-%M-%S")
    search_dir = os.path.join(base_path, subdir)
    pattern = os.path.join(search_dir, f"*{pattern_time}*.nc")

    matches = glob.glob(pattern)
    if matches:
        return matches[0]  # return the first match
    else:
        print(f"No IR file found for {dt} in {search_dir}")
        return None





from datetime import datetime
import xarray as xr

toocan_param_path="../../../config/fileparam_TOOCAN.dat"
geo_param_path="../../../config/fileparam_GEO.dat"

# === Load TOOCAN Parameters ===
params_TOOCAN = parse_param_file(toocan_param_path)

latmin = float(params_TOOCAN.get("latmin"))
latmax = float(params_TOOCAN.get("latmax"))
lonmin = float(params_TOOCAN.get("lonmin"))
lonmax = float(params_TOOCAN.get("lonmax"))

start_str = (
    f"{params_TOOCAN['yearBEGIN']:04d}-"
    f"{params_TOOCAN['monthBEGIN']:02d}-"
    f"{params_TOOCAN['dayBEGIN']:02d} "
    f"{params_TOOCAN['hourBEGIN']:02d}:"
    f"{params_TOOCAN['minBEGIN']:02d}"
)
start_time = datetime.strptime(start_str, "%Y-%m-%d %H:%M")

end_str = (
    f"{params_TOOCAN['yearEND']:04d}-"
    f"{params_TOOCAN['monthEND']:02d}-"
    f"{params_TOOCAN['dayEND']:02d} "
    f"{params_TOOCAN['hourEND']:02d}:"
    f"{params_TOOCAN['minEND']:02d}"
)
end_time = datetime.strptime(end_str, "%Y-%m-%d %H:%M")

version = str(params_TOOCAN.get("version"))

pathout_TOOCAN = params_TOOCAN.get("pathout_TOOCAN")
file_list = params_TOOCAN.get("file_list")


# === Load GEO parameters ===
params_GEO = parse_param_file(geo_param_path)
temporalresolution = int(params_GEO.get("temporalresolution"))
file_navigation = params_GEO.get("file_navigation")
region = params_GEO.get("REGION")
path_ir= params_GEO.get("path_ir")
# === Crop Navigation File ===
cropped_ds = crop_navigation_file(
    file_navigation,
    lonmin=lonmin, lonmax=lonmax,
    latmin=latmin, latmax=latmax
)
surface_area_2d = cropped_ds['mat_surfacePix'].values
lon_array_2d = cropped_ds['mat_longitude'].values
lat_array_2d = cropped_ds['mat_latitude'].values

# Initialize tracker
tracker = DCSLifeCycleTracker(surface_area_2d, lat_array_2d, lon_array_2d)


# === Chunking Logic ===
start_global = pd.to_datetime(start_time)
end_global = pd.to_datetime(end_time)

current_start = start_global

# === Define 30-minute increment ===
step = timedelta(minutes=30)
output_pathTOOCAN = f"{pathout_TOOCAN}toocan_{version}/{region}/"


while current_start < end_global:
    print(current_start,end_global)


    file_toocan = get_toocan_file(output_pathTOOCAN, current_start)
    file_ir   = get_GEO_IR_file(path_ir, current_start)
    print(file_ir,file_toocan)

    
    #with xr.open_dataset(file_toocan) as ds_label, xr.open_dataset(file_ir) as ds_ir:
    #    label_img = ds_label["ToocanCloudMask_"].values
        #irbt_img = ds_ir["Harmonized_irBT"].values


    result_toocan = read_toocan_label(file_toocan)
    label_img = result_toocan['label']
    result = read_irbt_subset(file_ir, lonmin, lonmax, latmin, latmax)
    irbt_img = result['irbt']

    #sys.exit()
#   
    tracker.update(label_img, irbt_img, current_start)
    # Example components

    current_start = current_start + step

# Finalize
#dcs_stats_LC = tracker.finalize()

create_dcs_netcdf_from_dicts(dcs_stats_LC, "TOOCAN_DCS_output.nc")

