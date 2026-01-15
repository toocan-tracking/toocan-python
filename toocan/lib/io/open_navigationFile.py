# =============================================================================
# File        : open_navigationFile.py
# Author      : Thomas Fiolleau
# Date        : 2025-07-04
# Description : Provides functionality to crop a navigation NetCDF file
#               to a user-defined latitude/longitude bounding box.
#
# Functions   :
#   - crop_navigation_file(): Reads, crops, and optionally saves a navigation file.
#
# Project     : TOOCAN - Tracking Organized Deep Convection
# License     : 
# =============================================================================
import xarray as xr
import numpy as np

def crop_navigation_file_v0(file_path, lonmin=0, lonmax=20, latmin=-5, latmax=15, save_path=None):
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
    print(lonmin,lonmax,latmin,latmax)
    print(np.shape(lat_var))
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


import xarray as xr
import numpy as np

def crop_navigation_file(navfile, lonmin, lonmax, latmin, latmax):
    """
    Crop navigation using ONLY min/max indices.
    Returns an xr.Dataset with:
      - mat_latitude(y,x)
      - mat_longitude(y,x)
      - mat_surfacePix(y,x)
      - xmin, xmax, ymin, ymax  (scalar coords)
    """

    ds = xr.open_dataset(navfile)

    # ------------------------------------------
    # 1) Identify latitude / longitude variables
    # ------------------------------------------
    lat_candidates = ['lat','latitude','mat_latitude','matlat']
    lon_candidates = ['lon','longitude','mat_longitude','matlon']
    area_candidates = ['mat_surfacePix','surface','area']

    lat_var = None
    lon_var = None
    area = None

    for key in lat_candidates:
        if key in ds:
            lat_var = ds[key].values
            break
    if lat_var is None:
        raise ValueError("No latitude field in navigation file")

    for key in lon_candidates:
        if key in ds:
            lon_var = ds[key].values
            break
    if lon_var is None:
        raise ValueError("No longitude field in navigation file")

    for key in area_candidates:
        if key in ds:
            area = ds[key].values
            break

    # ------------------------------------------
    # 2) Case A : regular grid  lat(y), lon(x)
    # ------------------------------------------
    if lat_var.ndim == 1 and lon_var.ndim == 1:

        y_idx = np.where((lat_var >= latmin) & (lat_var <= latmax))[0]
        x_idx = np.where((lon_var >= lonmin) & (lon_var <= lonmax))[0]

        if y_idx.size == 0 or x_idx.size == 0:
            raise ValueError("No pixel found inside requested domain")

        ymin, ymax = y_idx[0], y_idx[-1]
        xmin, xmax = x_idx[0], x_idx[-1]

        LAT = lat_var[ymin:ymax+1][:, None] * np.ones((1, xmax-xmin+1))
        LON = lon_var[xmin:xmax+1][None, :] * np.ones((ymax-ymin+1, 1))

        if area is None:
            area = np.ones_like(LAT)
        else:
            area = area[ymin:ymax+1, xmin:xmax+1]

    # ------------------------------------------
    # 3) Case B : 2D grid lat(y,x), lon(y,x)
    # ------------------------------------------
    elif lat_var.ndim == 2 and lon_var.ndim == 2:

        LAT = lat_var
        LON = lon_var

        mask = ((LAT >= latmin) & (LAT <= latmax) &
                (LON >= lonmin) & (LON <= lonmax))

        if not np.any(mask):
            raise ValueError("No pixel found inside requested domain")

        y_idx, x_idx = np.where(mask)
        ymin, ymax = y_idx.min(), y_idx.max()
        xmin, xmax = x_idx.min(), x_idx.max()

        LAT = LAT[ymin:ymax+1, xmin:xmax+1]
        LON = LON[ymin:ymax+1, xmin:xmax+1]

        if area is None:
            area = np.ones_like(LAT)
        else:
            area = area[ymin:ymax+1, xmin:xmax+1]

    else:
        raise ValueError("Unsupported navigation format")

    # ------------------------------------------
    # 4) Build output dataset with crop indices
    # ------------------------------------------
    ds_out = xr.Dataset(
        data_vars={
            "mat_latitude":   (("y", "x"), LAT),
            "mat_longitude":  (("y", "x"), LON),
            "mat_surfacePix": (("y", "x"), area),
        },
        coords={
            "xmin": int(xmin),
            "xmax": int(xmax),
            "ymin": int(ymin),
            "ymax": int(ymax),
        }
    )

    return ds_out