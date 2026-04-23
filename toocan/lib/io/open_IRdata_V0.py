# =============================================================================
# File        : open_IRdata.py
# Description : Fast and simple IR/BT/OLR readers for multiple satellite/models.
#               Uses a ROUTER + fast dedicated readers.
# =============================================================================

import os
import numpy as np
import pandas as pd
import xarray as xr

SIGMA = 5.670374419e-8  # Stefan-Boltzmann constant

# ============================================================================
# FAST HELPERS
# ============================================================================

def crop_latlon_fast(da, lat_name, lon_name, latmin, latmax, lonmin, lonmax):
    """Fast lat/lon cropping without metadata overhead."""
    lat = da[lat_name].values
    if lat[0] > lat[-1]:
        return da.sel(
            **{
                lat_name: slice(latmax, latmin),
                lon_name: slice(lonmin, lonmax),
            }
        )
    else:
        return da.sel(
            **{
                lat_name: slice(latmin, latmax),
                lon_name: slice(lonmin, lonmax),
            }
        )


def convert_olr_to_irbt(da):
    """Convert OLR (W/m²) to IRBT (K)."""
    olr = da.astype(np.float32)
    Tb = (olr / SIGMA)**0.25
    return Tb**2 * 0.00443828 - Tb * 0.543846 + 129.544


# ============================================================================
# FAST MODEL-SPECIFIC READERS
# ============================================================================

# ---------------------- MSG SEVIRI -------------------------
# ---------------------- MSG SEVIRI -------------------------
def read_msg_irbt(filepath, lonmin, lonmax, latmin, latmax):
    """Very fast MSG reader (Harmonized_irBT)."""

    ds = xr.open_dataset(filepath)

    da = ds["brighness_temperature"].isel(time=0).astype(np.float32)       # or Harmonized_irBT with correction or irBT without correction

    da = crop_latlon_fast(da, "latitude", "longitude",
                          latmin, latmax, lonmin, lonmax)

    scale = da.attrs.get("scale_factor", 1.0)
    fill = da.attrs.get("_FillValue", None)

    arr = da.values * scale
    if fill is not None:
        arr = np.where(da.values == fill, np.nan, arr)

    return dict(
        irbt=arr,
        lat=da["latitude"].values,
        lon=da["longitude"].values,
        time=pd.to_datetime(ds["time"].values[0])
    )



# ---------------------- MSG SEVIRI -------------------------
# ---------------------- MSG SEVIRI -------------------------
def read_msg_irbt_native(filepath, lonmin, lonmax, latmin, latmax):
    """Very fast MSG reader (Harmonized_irBT)."""

    ds = xr.open_dataset(filepath)

    # -------------------------------
    # Extract IR108 brightness temp
    # -------------------------------
    da = ds["IR_108"].astype(np.float32)

    scale = da.attrs.get("scale_factor", 1.0)
    offset = da.attrs.get("add_offset", 0.0)
    fill = da.attrs.get("_FillValue", None)

    irbt_raw = da.values * scale + offset
    if fill is not None:
        irbt_raw = np.where(da.values == fill, np.nan, irbt_raw)
        
    da = crop_latlon_fast(da, "latitude", "longitude",
                          latmin, latmax, lonmin, lonmax)

    scale = da.attrs.get("scale_factor", 1.0)
    add_offset = da.attrs.get("add_offset", 0.0)
    fill = da.attrs.get("_FillValue", None)

    arr = da.values * scale
    if fill is not None:
        arr = np.where(da.values == fill, np.nan, arr)

    # -------------------------------
    # Read time (safe)
    # -------------------------------
    try:
        unix = float(ds["time"].values)
        time = datetime.utcfromtimestamp(unix)
    except:
        # fallback: parse filename
        import re
        m = re.search(r"(\d{8})_(\d{6})", filepath)
        if m:
            time = datetime.strptime(m.group(1)+m.group(2), "%Y%m%d%H%M%S")
        else:
            time = None

    return dict(
        irbt=arr,
        lat=da["ny"].values,
        lon=da["nx"].values,
        time=time
    )


# ---------------------- ARPEGE-NH -------------------------
def read_arpege_bt(filepath, lonmin, lonmax, latmin, latmax):
    """Reader for ARPEGE NH: BT(time, spec, latitude, longitude)."""

    ds = xr.open_dataset(filepath)

    da = ds["BT"].isel(time=0, spec=0).astype(np.float32)

    da = crop_latlon_fast(da, "latitude", "longitude",
                          latmin, latmax, lonmin, lonmax)

    return dict(
        irbt=da.values,
        lat=da["latitude"].values,
        lon=da["longitude"].values,
        time=pd.to_datetime(ds["time"].values[0])
    )


# ---------------------- ICON MODEL -------------------------
def read_icon_bt(filepath, lonmin, lonmax, latmin, latmax):
    """Example ICON reader (adjust variable names as needed)."""
    ds = xr.open_dataset(filepath)

    da = ds["BT"].isel(time=0).astype(np.float32)

    da = crop_latlon_fast(da, "lat", "lon",
                          latmin, latmax, lonmin, lonmax)

    return dict(
        irbt=da.values,
        lat=da["lat"].values,
        lon=da["lon"].values,
        time=pd.to_datetime(ds["time"].values[0])
    )


# ---------------------- RLUT / OLR-GENERIC -------------------------
def read_sam_olr(filepath, lonmin, lonmax, latmin, latmax):
    """Reader for OLR fields converting to IRBT."""

    ds = xr.open_dataset(filepath)

    if "rlut" in ds:
        da = ds["rlut"].isel(time=0)
        lat_name = "lat"
        lon_name = "lon"
    elif "LWNTA" in ds:
        da = ds["LWNTA"].isel(time=0)
        lat_name = "lat"
        lon_name = "lon"
    else:
        raise ValueError(f"No RLUT/LWNTA in file: {filepath}")

    da = da.astype(np.float32)

    da = crop_latlon_fast(da, lat_name, lon_name,
                          latmin, latmax, lonmin, lonmax)

    irbt = convert_olr_to_irbt(da)

    return dict(
        irbt=irbt.values.astype(np.float32),
        lat=da[lat_name].values,
        lon=da[lon_name].values,
        time=pd.to_datetime(ds["time"].values[0])
    )


# ============================================================================
# ROUTER — CHOOSE MODEL READER
# ============================================================================

def read_ir_file(filepath, lonmin, lonmax, latmin, latmax, model_name):
    """
    Routeur principal : envoie au lecteur rapide selon model_name.
    """
    if model_name == "MSG":
        return read_msg_irbt(filepath, lonmin, lonmax, latmin, latmax)

    if model_name == "MSG_native":
        return read_msg_irbt_native(filepath, lonmin, lonmax, latmin, latmax)

    elif model_name == "ARPEGENH":
        return read_arpege_bt(filepath, lonmin, lonmax, latmin, latmax)

    elif model_name == "ICON":
        return read_icon_olr(filepath, lonmin, lonmax, latmin, latmax)

    elif model_name == "SAM":
        return read_sam_olr(filepath, lonmin, lonmax, latmin, latmax)

    else:
        raise ValueError(f"Unknown model_name: {model_name}")


# ============================================================================
# VOLUMETRIC EXTRACTOR (called by main.py)
# ============================================================================

def extract_volume(df, start_time, end_time,
                   lonmin, lonmax, latmin, latmax,
                   model_name):
    """
    Build (T, lat, lon) volume from files and timestamps.
    """

    start_dt = pd.to_datetime(start_time)
    end_dt = pd.to_datetime(end_time)

    df_sel = df[(df["datetime"] >= start_dt) &
                (df["datetime"] <= end_dt)].copy()

    volume_list = []
    time_list = []
    lat_ref, lon_ref = None, None

    print(f" Time window: {start_dt} → {end_dt}")
    print(f" Files in range: {len(df_sel)}\n")

    for _, row in df_sel.iterrows():
        full_path = row["full_path"]
        print("→ Reading:", full_path)

        try:
            result = read_ir_file(
                full_path, lonmin, lonmax, latmin, latmax,
                model_name=model_name
            )

            volume_list.append(result["irbt"])
            time_list.append(result["time"])

            if lat_ref is None:
                lat_ref = result["lat"]
                lon_ref = result["lon"]

        except Exception as e:
            print(f"❌ Could not read {full_path}: {e}")

    volume = np.stack(volume_list, axis=0) if volume_list else np.array([])

    return volume, time_list, lat_ref, lon_ref


# ============================================================================
# LEGACY FUNCTION (optional)
# ============================================================================
def read_irbt_subset(filepath, lonmin, lonmax, latmin, latmax):
    """Compatibility wrapper kept for old code (MSG only)."""
    return read_msg_irbt(filepath, lonmin, lonmax, latmin, latmax)