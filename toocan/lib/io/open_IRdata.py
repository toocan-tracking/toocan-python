# =============================================================================
# File        : open_IRdata.py
# Description : IR/BT/OLR readers using navigation-crop indices (pixel-accurate)
# =============================================================================

import os
import numpy as np
import pandas as pd
import xarray as xr
from datetime import datetime

SIGMA = 5.670374419e-8  # Stefan–Boltzmann constant

def get_xy_dims(da):
    """
    Return the (y_dim, x_dim) from the DataArray.
    Accepts dim names among:
      ('y','x'), ('ny','nx'), ('latitude','longitude'), ('lat','lon')
    """
    Y_CAND = ["y", "ny", "latitude", "lat"]
    X_CAND = ["x", "nx", "longitude", "lon"]

    dims = da.dims

    y_dim = next((d for d in dims if d.lower() in Y_CAND), None)
    x_dim = next((d for d in dims if d.lower() in X_CAND), None)

    if y_dim is None or x_dim is None:
        raise ValueError(f"Cannot identify 2D dims in {dims}")

    return y_dim, x_dim

# ======================================================================
# COLD CLOUD / OLR CONVERSION
# ======================================================================

def convert_olr_to_irbt(da):
    """Convert TOA OLR (W/m²) to IR brightness temperature (K)."""
    olr = da.astype(np.float32)
    Tb = (olr / SIGMA)**0.25
    return Tb**2 * 0.00443828 - Tb * 0.543846 + 129.544


# ======================================================================
# NEW: Cropping IR using navigation-output
# ======================================================================
def crop_with_nav(da, nav):
    """
    Crop DataArray da using navigation pixel indices.
    Works with dim names y/x or ny/nx.
    """

    y_dim, x_dim = get_xy_dims(da)

    ymin = int(nav["ymin"])
    ymax = int(nav["ymax"])
    xmin = int(nav["xmin"])
    xmax = int(nav["xmax"])

    return da.isel({y_dim: slice(ymin, ymax+1),
                    x_dim: slice(xmin, xmax+1)})
    
# ======================================================================
# =====================     MSG  SEVIRI  ===============================
# ======================================================================

def read_msg_irbt(filepath, nav, model_name=None):
    """
    Read Harmonized_irBT from MSG HRIT–>NETCDF product.
    Cropping is done using navigation indices (pixel exact).
    """

    ds = xr.open_dataset(filepath)

    da = ds["Harmonized_irBT"].astype(np.float32)
    if da.ndim == 3:
        da = da.isel(time=0)

    # Pixel-accurate crop
    da = crop_with_nav(da, nav)

    scale = da.attrs.get("scale_factor", 1.0)
    fill  = da.attrs.get("_FillValue", None)

    arr_raw = da.values
    arr = arr_raw * scale

    if fill is not None:
        arr = np.where(arr_raw == fill, np.nan, arr)

    # Time
    t = pd.to_datetime(ds["time"].values[0])

    return dict(
        irbt=arr.astype(np.float32),
        lat=nav["mat_latitude"].values,
        lon=nav["mat_longitude"].values,
        time=t
    )


# ======================================================================
# ================= MSG HRIT native (IR_108) ===========================
# ======================================================================
import os
import re
import numpy as np
import xarray as xr
from datetime import datetime
import pandas as pd


def _extract_hrit_time(filepath):
    """
    Extract time from HRIT native filename.
    Examples:
       *_20230319_2130.nc
       *_202303192130.nc
    """

    name = os.path.basename(filepath)

    # Pattern 1: YYYYMMDD_HHMM
    m = re.search(r"(\d{8})_(\d{4})", name)
    if m:
        return datetime.strptime(m.group(1) + m.group(2), "%Y%m%d%H%M")

    # Pattern 2: YYYYMMDDHHMM
    m = re.search(r"(\d{8})(\d{4})", name)
    if m:
        return datetime.strptime(m.group(1) + m.group(2), "%Y%m%d%H%M")

    # Pattern 3: YYYYMMDDHHMMSS
    m = re.search(r"(\d{8})(\d{6})", name)
    if m:
        return datetime.strptime(m.group(1) + m.group(2), "%Y%m%d%H%M%S")

    print(f"⚠ WARNING: cannot read time from filename: {name}")
    return None


def read_msg_irbt_native(filepath, nav, model_name=None):
    """
    Read MF HRIT-native (converted) IR_108 and crop using a pre-cut navigation.
    
    Parameters
    ----------
    filepath : str
         HRIT-native IR file
    nav : xr.Dataset
         Contains 'mat_latitude', 'mat_longitude', 'mat_surfacePix'
    """

    ds = xr.open_dataset(filepath)

    # --- select IR108 ---
    da = ds["IR_108"].astype(np.float32)
    if da.ndim == 3:
        da = da.isel(time=0)

    scale  = da.attrs.get("scale_factor", 1.0)
    offset = da.attrs.get("add_offset", 0.0)
    fill   = da.attrs.get("_FillValue", None)

    # --- crop using navigation ---
    # nav already contains the reduced domain
    da_crop = crop_with_nav(da, nav)

    raw = da_crop.values
    irbt = raw * scale + offset
    if fill is not None:
        irbt = np.where(raw == fill, np.nan, irbt)

    # --- get time ---
    t = None

    # Try CF time first
    if "time" in ds:
        try:
            t = pd.to_datetime(ds["time"].values[0])
        except:
            t = None

    # Fallback to filename
    if t is None:
        t = _extract_hrit_time(filepath)

    if t is None:
        raise ValueError(f"❌ Could not determine time for HRIT file {filepath}")

    return dict(
        irbt=irbt.astype(np.float32),
        lat=nav["mat_latitude"].values,
        lon=nav["mat_longitude"].values,
        time=t
    )

# ======================================================================
# ==================== ARPEGE NH MODEL ================================
# ======================================================================

def read_arpege_bt(filepath, nav, model_name=None):
    """ARPEGE NH BT(time, spec, lat, lon)."""

    ds = xr.open_dataset(filepath)

    da = ds["BT"].isel(time=0, spec=0).astype(np.float32)

    da = crop_with_nav(da, nav)

    return dict(
        irbt=da.values.astype(np.float32),
        lat=nav["mat_latitude"].values,
        lon=nav["mat_longitude"].values,
        time=pd.to_datetime(ds["time"].values[0])
    )


# ======================================================================
# ========================== ICON MODEL ================================
# ======================================================================

def read_icon_bt(filepath, nav, model_name=None):
    """ICON BT(time, lat, lon)."""

    ds = xr.open_dataset(filepath)

    da = ds["BT"].isel(time=0).astype(np.float32)
    da = crop_with_nav(da, nav)

    return dict(
        irbt=da.values.astype(np.float32),
        lat=nav["mat_latitude"].values,
        lon=nav["mat_longitude"].values,
        time=pd.to_datetime(ds["time"].values[0])
    )


# ======================================================================
# ====================== OLR → IRBT PRODUCTS ===========================
# ======================================================================

def read_sam_olr(filepath, nav, model_name=None):
    """Generic OLR reader converting to IRBT."""

    ds = xr.open_dataset(filepath)

    if "rlut" in ds:
        da = ds["rlut"].isel(time=0)
    elif "LWNTA" in ds:
        da = ds["LWNTA"].isel(time=0)
    else:
        raise ValueError(f"No RLUT/LWNTA in {filepath}")

    da = crop_with_nav(da, nav)
    irbt = convert_olr_to_irbt(da)

    return dict(
        irbt=irbt.values.astype(np.float32),
        lat=nav["mat_latitude"].values,
        lon=nav["mat_longitude"].values,
        time=pd.to_datetime(ds["time"].values[0])
    )


# ======================================================================
# ROUTER
# ======================================================================

def read_ir_file(filepath, nav, model_name):
    """
    Universal IR reader: delegates to correct backend.
    """

    if model_name == "MSG":
        return read_msg_irbt(filepath, nav)

    elif model_name == "MSG_native":
        return read_msg_irbt_native(filepath, nav)

    elif model_name == "ARPEGENH":
        return read_arpege_bt(filepath, nav)

    elif model_name == "ICON":
        return read_icon_bt(filepath, nav)

    elif model_name == "SAM":
        return read_sam_olr(filepath, nav)

    else:
        raise ValueError(f"Unknown model_name {model_name}")


# ======================================================================
# VOLUME EXTRACTOR
# ======================================================================

def extract_volume(df, start_time, end_time, nav, model_name):
    """
    Extract (T, y, x) volume cropped according to nav.
    """

    df_sel = df[(df["datetime"] >= start_time) &
                (df["datetime"] <= end_time)].copy()

    volume_list = []
    time_list = []

    for _, row in df_sel.iterrows():
        fp = row["full_path"]

        try:
            r = read_ir_file(fp, nav, model_name)
            volume_list.append(r["irbt"])
            time_list.append(r["time"])

        except Exception as e:
            print(f"❌ Failed reading {fp}: {e}")

    if len(volume_list) == 0:
        return np.array([]), [], nav["mat_latitude"].values, nav["mat_longitude"].values

    volume = np.stack(volume_list, axis=0)

    return volume, time_list, nav["mat_latitude"].values, nav["mat_longitude"].values