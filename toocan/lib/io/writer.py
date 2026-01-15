# =============================================================================
# File        : writer_toocanImage.py
# Author      : Thomas Fiolleau
# Date        : 2025-07-04
# Description : Save TOOCAN segmented label volumes as CF-compliant NetCDF images,
#               one per time slot, with appropriate metadata and attributes.
#
# Functions   :
#   - save_labels_slot_by_slot(): Export each time slice of a 3D label volume as
#                                 a separate NetCDF file with TOOCAN-compliant metadata.
#
# Project     : TOOCAN - Tracking Organized Deep Convection
# License     : MIT (or specify your license)
# =============================================================================

import os
import xarray as xr
import numpy as np
import pandas as pd
from datetime import datetime

def save_labels_slot_by_slot(label_volume, time_array, lat_array, lon_array,
                              var_name="DCS_number",
                              params_TOOCAN=None, params_GEO=None):
    """
    Save 3D label volume as individual NetCDF files, one per time step,
    including CF-compliant time metadata and TOOCAN-style attributes.

    Parameters:
        label_volume (np.ndarray): 3D array of shape (time, lat, lon)
        time_array (np.ndarray): 1D array of datetime objects
        lat_array (np.ndarray): 1D array of latitude values
        lon_array (np.ndarray): 1D array of longitude values
        output_dir (str): Directory where NetCDF files will be saved
        var_name (str): Name of the variable to save
        params_TOOCAN (dict): Dictionary of TOOCAN parameters (optional)
        params_GEO (dict): Dictionary of GEO platform parameters (optional)
    """

    n_time = label_volume.shape[0]

    # Get metadata from param files (with defaults)
    platform = (params_GEO or {}).get("GEOplatform", "MSG")
    channel = (params_GEO or {}).get("channel", "IR108")
    nadir = str((params_GEO or {}).get("nadir", 0.0))
    resolution_deg = float((params_GEO or {}).get("spatialresolution", 0.04))
    version = str((params_TOOCAN or {}).get("version", "2.08"))
    scan_duration_min = int((params_GEO or {}).get("temporalresolution", 30))

    output_dir = params_TOOCAN.get("pathout_TOOCAN")
    os.makedirs(output_dir, exist_ok=True)

    output_dir =output_dir+'/toocan_'+version+'/'
    os.makedirs(output_dir, exist_ok=True)

    for i in range(n_time):
        tstamp = pd.to_datetime(time_array[i]) if not isinstance(time_array[i], datetime) else time_array[i]
        filename = f"{var_name}_{tstamp:%Y%m%d_%H%M}.nc"
        path = os.path.join(output_dir, filename)

        # Prepare data
        data = label_volume[i].astype(np.int32)
        data[data == 0] = -998  # Replace zeros with fill value

        seconds_since_epoch = (tstamp - datetime(1970, 1, 1)).total_seconds()

        da = xr.DataArray(
            data[np.newaxis, :, :],
            coords={
                "time": [seconds_since_epoch],
                "latitude": lat_array,
                "longitude": lon_array
            },
            dims=["time", "latitude", "longitude"],
            name=var_name,
            attrs={
                "units": " ",
                "standard_name": "number_of_Deep_Cloud_Systems",
                "long_name": "number of Deep Cloud Systems",
                "_FillValue": -998,
                "missing_value": -999
            }
        )

        ds = xr.Dataset({var_name: da})

        ds["time"].attrs.update({
            "standard_name": "time",
            "long_name": "scan start time",
            "units": "seconds since 1970-01-01 00:00:00 UTC"
        })

        ds["latitude"].attrs.update({
            "units": "degrees_north",
            "long_name": "latitude"
        })

        ds["longitude"].attrs.update({
            "units": "degrees_east",
            "long_name": "longitude",
            "FlagDegEst": 0
        })

        ds.attrs.update({
            "title": f"TOOCAN segmented images - Gridded data {resolution_deg:.3f}°",
            "conventions": "CF-1.6",
            "institution": "CNRS/LEGOS/IPSL",
            "creator_name": "Thomas Fiolleau",
            "image_time": tstamp.strftime("%Y-%m-%d-T%H-%M-%S UTC"),
            "scan_start_time": tstamp.strftime("%H-%M-%S UTC"),
            "scan_end_time": (tstamp + pd.Timedelta(minutes=scan_duration_min)).strftime("%H-%M-%S UTC"),
            "platform": platform,
            "channel": f"Infrared {channel}",
            "nadir": f"sub satellite point longitude: {nadir}°",
            "grid_resolution_in_degrees": resolution_deg,
            "version": version
        })

        ds.to_netcdf(path, unlimited_dims=["time"])
        print(f"✅ Saved {path}")