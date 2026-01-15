import os
import numpy as np
import xarray as xr
from datetime import datetime


def save_navigation_grid(nav_ds, params_TOOCAN, params_GEO, filename="navigation_grid.nc"):
    """
    Save the cropped navigation grid (lat, lon, surface pixel area)
    into the TOOCAN output directory:

        <pathout_TOOCAN>/navigation/<REGION>/navigation_grid.nc

    Parameters
    ----------
    nav_ds : xr.Dataset
        Output of crop_navigation_file(), containing:
        - mat_latitude (y,x)
        - mat_longitude (y,x)
        - mat_surfacePix (y,x)
        attrs:
        - xmin, xmax, ymin, ymax  (crop indices)

    params_TOOCAN : dict
        TOOCAN parameter file contents (parsed)

    params_GEO : dict
        GEO parameter file contents (parsed)

    filename : str
        Name of the navigation file to write (default: navigation_grid.nc)
    """

    # ----------------------------
    # 1. Build output directory
    # ----------------------------
    base_out = params_TOOCAN.get("pathout_TOOCAN")
    if isinstance(base_out, bytes):   # handle DataParam bytes
        base_out = base_out.decode()

    version = str(params_TOOCAN.get("version", "2.08"))
    region  = str(params_GEO.get("REGION", "GLOBAL"))

    # <pathout_TOOCAN>/navigation/<REGION>/
    outdir = os.path.join(base_out,"toocan_"+version,region, "navigation")
    os.makedirs(outdir, exist_ok=True)

    outfile = os.path.join(outdir, filename)

    print(f"🛰 Saving navigation file → {outfile}")

    # ----------------------------
    # 2. Extract arrays
    # ----------------------------
    LAT  = nav_ds["mat_latitude"].values.astype(np.float32)
    LON  = nav_ds["mat_longitude"].values.astype(np.float32)
    AREA = nav_ds["mat_surfacePix"].values.astype(np.float32)

    ny, nx = LAT.shape

    # Crop indices
    xmin = int(nav_ds.coords["xmin"].values)
    xmax = int(nav_ds.coords["xmax"].values)
    ymin = int(nav_ds.coords["ymin"].values)
    ymax = int(nav_ds.coords["ymax"].values)

    # Spatial resolution (if available)
    spatialres = float(params_GEO.get("spatialresolution", 0.05))

    # ----------------------------
    # 3. Build Dataset
    # ----------------------------
    ds_out = xr.Dataset(
        data_vars={
            "mat_latitude":    (("y", "x"), LAT),
            "mat_longitude":   (("y", "x"), LON),
            "mat_surfacePix":  (("y", "x"), AREA),
        },
        coords={
            "y": np.arange(ny, dtype=np.int32),
            "x": np.arange(nx, dtype=np.int32)
        },
        attrs={
            "title": "Cropped Navigation Grid for TOOCAN",
            "creation_time": datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ"),
            "platform": params_GEO.get("GEOplatform", "UNKNOWN"),
            "region": region,
            "version": version,
            "spatial_resolution_deg": spatialres,
            "xmin": xmin,
            "xmax": xmax,
            "ymin": ymin,
            "ymax": ymax
        }
    )

    # ----------------------------
    # 4. Save NetCDF
    # ----------------------------
    print(outfile)
    ds_out.to_netcdf(outfile)
    print(f"   ✔ Navigation saved ({nx}×{ny})")
    return outfile