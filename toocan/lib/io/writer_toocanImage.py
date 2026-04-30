import os
import xarray as xr
import numpy as np
import pandas as pd
from datetime import datetime

def save_labels_slot_by_slot_v0(label_volume, time_array, lat_array, lon_array,
                              nomenclature="ToocanCloudMask_",var_name="DCS_number",
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
    region  = str(params_GEO.get("REGION"))
    scan_duration_min = int((params_GEO or {}).get("temporalresolution", 30))

    output_dir = params_TOOCAN.get("pathout_TOOCAN")
    os.makedirs(output_dir, exist_ok=True)

    output_dir =output_dir+'/toocan_'+version+'/'
    os.makedirs(output_dir, exist_ok=True)

    output_dir =output_dir+'/'+region+'/'
    #output_dir = os.path.join(output_dir, region)
    os.makedirs(output_dir, exist_ok=True)
    base_output_dir = output_dir  # keep original

    for i in range(n_time):
        tstamp = pd.to_datetime(time_array[i]) if not isinstance(time_array[i], datetime) else time_array[i]

        # Extract year, month, day
        year = str(tstamp.year)  # convert to string
        month = f"{tstamp.month:02d}"
        day = f"{tstamp.day:02d}"

        # Construct year-level path
        year_dir = os.path.join(base_output_dir, year)
        if not os.path.exists(year_dir):
            os.makedirs(year_dir)
    
        # Construct full date-level path (e.g., 2023/2023_07_16)
        day_dir = os.path.join(year_dir, f"{year}_{month}_{day}")
        if not os.path.exists(day_dir):
            os.makedirs(day_dir)

        filename = f"{nomenclature}{tstamp:%Y%m%d_%H%M}.nc"
        path = os.path.join(day_dir, filename)

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




def save_labels_slot_by_slot_v1(label_volume, time_array, lat_array, lon_array,
                              nomenclature="ToocanCloudMask_",var_name="DCS_number",
                              params_TOOCAN=None, params_GEO=None):
    """
    Grid-aware TOOCAN mask writer:
      - MSG_native  (HRIT)  → save only (ny,nx) coords
      - Regular grid        → save lat(y) / lon(x)
      - 2D grid (navigation)-> save lat(y,x) / lon(y,x)
    """

    import os
    import xarray as xr
    import numpy as np

    outdir = params_TOOCAN["pathout_TOOCAN"]
    os.makedirs(outdir, exist_ok=True)

    T, ny, nx = label_volume.shape

    for i in range(T):

        slot = label_volume[i]
        t = times[i]

        fname = f"{nomenclature}{t:%Y%m%d_%H%M}.nc"
        fpath = os.path.join(outdir, fname)

        # ======================================================
        # CASE 1 — MSG HRIT native → DO NOT SAVE LAT/LON
        # ======================================================
        if model_name == "MSG_native":

            da = xr.DataArray(
                slot.astype(np.int32),
                dims=("y", "x"),
                coords={
                    "y": np.arange(ny),
                    "x": np.arange(nx)
                },
                name="DCS_number",
            )

            ds = xr.Dataset({"DCS_number": da})
            ds["time"] = xr.DataArray([np.int64(t.timestamp())], dims=("time",))

            ds.to_netcdf(fpath)
            print(f"✔ MSG-native TOOCAN mask saved: {fpath}")
            continue

        # ======================================================
        # CASE 2 — Regular grid (1D lat, 1D lon)
        # ======================================================
        if lat_array.ndim == 1 and lon_array.ndim == 1:

            da = xr.DataArray(
                slot.astype(np.int32),
                dims=("y", "x"),
                coords={
                    "y": np.arange(ny),
                    "x": np.arange(nx),
                    "latitude":  ("y", lat_array),
                    "longitude": ("x", lon_array),
                },
                name="DCS_number",
            )

            ds = xr.Dataset({"DCS_number": da})
            ds["time"] = xr.DataArray([np.int64(t.timestamp())], dims=("time",))
            ds.to_netcdf(fpath)
            print(f"✔ Regular-grid TOOCAN mask saved: {fpath}")
            continue

        # ======================================================
        # CASE 3 — Navigation 2D grid
        # ======================================================
        if lat_array.ndim == 2 and lon_array.ndim == 2:

            da = xr.DataArray(
                slot.astype(np.int32),
                dims=("y", "x"),
                coords={
                    "y": np.arange(ny),
                    "x": np.arange(nx),
                    "latitude":  (("y", "x"), lat_array),
                    "longitude": (("y", "x"), lon_array),
                },
                name="DCS_number",
            )

            ds = xr.Dataset({"DCS_number": da})
            ds["time"] = xr.DataArray([np.int64(t.timestamp())], dims=("time",))
            ds.to_netcdf(fpath)
            print(f"✔ 2D-grid TOOCAN mask saved: {fpath}")
            continue

        raise ValueError("Unrecognized grid format")        


import os
import numpy as np
import xarray as xr
import pandas as pd
from datetime import datetime


def save_labels_slot_by_slot(
    label_volume, time_array,
    lat_array, lon_array,
    nomenclature="ToocanCloudMask_",
    var_name="DCS_number",
    params_TOOCAN=None, params_GEO=None
):
    """
    Universal TOOCAN mask writer supporting:
    - MSG_native  → save (time,y,x) only
    - Regular grid → lat(y), lon(x)
    - Navigation 2D grid → lat(y,x), lon(y,x)
    And including rich CF metadata.
    """

    n_time, ny, nx = label_volume.shape

    # ---------------------------------------------------------------------
    # PARAMETERS AND OUTPUT DIRECTORY CONSTRUCTION
    # ---------------------------------------------------------------------
    params_TOOCAN = params_TOOCAN or {}
    params_GEO = params_GEO or {}

    platform   = params_GEO.get("GEOplatform", "MSG")
    region     = params_GEO.get("REGION", "UNKNOWN")
    channel    = params_GEO.get("channel", "IR108")
    nadir      = params_GEO.get("nadir", 0.0)
    res_deg    = float(params_GEO.get("spatialresolution", 0.04))
    version    = str(params_TOOCAN.get("version", "2.08"))
    scan_dur   = int(params_GEO.get("temporalresolution", 30))

    root = params_TOOCAN.get("pathout_TOOCAN")
    if root is None:
        raise ValueError("Missing params_TOOCAN['pathout_TOOCAN']")

    outdir = os.path.join(root, f"toocan_{version}", region)
    os.makedirs(outdir, exist_ok=True)

    # ---------------------------------------------------------------------
    # DETERMINE GRID TYPE FOR WRITING
    # ---------------------------------------------------------------------
    if platform == "MSG_native":
        grid_type = "MSG_NATIVE"
    elif platform == "MSGrss":
        grid_type = "MSG_RSS"
    elif lat_array.ndim == 1 and lon_array.ndim == 1:
        grid_type = "REGULAR"
    elif lat_array.ndim == 2 and lon_array.ndim == 2:
        grid_type = "GRID2D"
    else:
        raise ValueError("Unsupported lat/lon format")

    print(f"🛰  Writer grid type: {grid_type}")

    # ---------------------------------------------------------------------
    # LOOP ON EACH TIME SLOT
    # ---------------------------------------------------------------------
    for i in range(n_time):

        tstamp = pd.to_datetime(time_array[i])
        year  = f"{tstamp.year}"
        month = f"{tstamp.month:02d}"
        day   = f"{tstamp.day:02d}"

        # Directory structure .../region/YYYY/YYYY_MM_DD/
        ydir  = os.path.join(outdir, year)
        ddir  = os.path.join(ydir, f"{year}_{month}_{day}")
        os.makedirs(ddir, exist_ok=True)

        fname = f"{nomenclature}{tstamp:%Y%m%d_%H%M}.nc"
        fpath = os.path.join(ddir, fname)

        # Replace zeros with fill value
        data = label_volume[i].astype(np.int32)
        data[data == 0] = -998

        # CF-compliant time
        epoch_sec = float((tstamp - datetime(1970,1,1)).total_seconds())

        # =================================================================
        # CASE 1 — MSG_NATIVE → Only (time,y,x)
        # =================================================================
        if grid_type == "MSG_NATIVE" or grid_type == "MSG_RSS":
            da = xr.DataArray(
                data[np.newaxis, :, :],
                dims=("time", "y", "x"),
                coords={
                    "time": [epoch_sec],
                    "y": np.arange(ny),
                    "x": np.arange(nx),
                },
                name=var_name,
                attrs={
                    "long_name": "number of Deep Cloud Systems",
                    "standard_name": "number_of_Deep_Cloud_Systems",
                    "_FillValue": -998,
                    "missing_value": -999,
                }
            )

            ds = xr.Dataset({var_name: da})

        # =================================================================
        # CASE 2 — REGULAR GRID (lat(y), lon(x))
        # =================================================================
        elif grid_type == "REGULAR":

            da = xr.DataArray(
                data[np.newaxis, :, :],
                dims=("time", "latitude", "longitude"),
                coords={
                    "time": [epoch_sec],
                    "latitude": ("latitude", lat_array),
                    "longitude": ("longitude", lon_array),
                },
                name=var_name,
                attrs={
                    "long_name": "number of Deep Cloud Systems",
                    "standard_name": "number_of_Deep_Cloud_Systems",
                    "_FillValue": -998,
                    "missing_value": -999,
                }
            )

            ds = xr.Dataset({var_name: da})

            ds["latitude"].attrs.update({
                "units": "degrees_north",
                "long_name": "latitude"
            })
            ds["longitude"].attrs.update({
                "units": "degrees_east",
                "long_name": "longitude"
            })

        # =================================================================
        # CASE 3 — 2D GRID (lat(y,x), lon(y,x))
        # =================================================================
        elif grid_type == "GRID2D":

            da = xr.DataArray(
                data[np.newaxis, :, :],
                dims=("time", "y", "x"),
                coords={
                    "time": [epoch_sec],
                    "latitude":  (("y","x"), lat_array.astype(np.float32)),
                    "longitude": (("y","x"), lon_array.astype(np.float32)),
                },
                name=var_name,
                attrs={
                    "long_name": "number of Deep Cloud Systems",
                    "standard_name": "number_of_Deep_Cloud_Systems",
                    "_FillValue": -998,
                    "missing_value": -999,
                }
            )

            ds = xr.Dataset({var_name: da})

            ds["latitude"].attrs.update({
                "units": "degrees_north",
                "long_name": "latitude"
            })
            ds["longitude"].attrs.update({
                "units": "degrees_east",
                "long_name": "longitude"
            })

        # -----------------------------------------------------------------
        # Global attributes
        # -----------------------------------------------------------------
        ds["time"].attrs.update({
            "standard_name": "time",
            "long_name": "scan start time",
            "units": "seconds since 1970-01-01 00:00:00 UTC"
        })

        ds.attrs.update({
            "title": f"TOOCAN segmented images - Gridded data {res_deg:.3f}°",
            "conventions": "CF-1.6",
            "institution": "CNRS/LEGOS/IPSL",
            "creator_name": "Thomas Fiolleau",
            "image_time": tstamp.strftime("%Y-%m-%d-T%H-%M-%S UTC"),
            "scan_start_time": tstamp.strftime("%H-%M-%S UTC"),
            "scan_end_time": (tstamp + pd.Timedelta(minutes=scan_dur)).strftime("%H-%M-%S UTC"),
            "platform": platform,
            "channel": f"Infrared {channel}",
            "nadir": f"sub satellite point longitude: {nadir}°",
            "grid_resolution_in_degrees": res_deg,
            "version": version
        })


        encoding = {
            var: {
                "zlib": True,
                "complevel": 4,
                "shuffle": True
            }
            for var in ds.data_vars
        }

        ds.to_netcdf(fpath, encoding=encoding, unlimited_dims=["time"])
        print(f"✔ Saved {fpath}")