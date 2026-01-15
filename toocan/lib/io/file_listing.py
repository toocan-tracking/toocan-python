# =============================================================================
# File        : datetime_utils.py
# Author      : Thomas Fiolleau
# Project     : TOOCAN - Tracking Organized Deep Convection
# Description : Helper functions to extract datetimes from complex filenames
#               used by meteorological satellites (MSG, Himawari, GOES, INSAT,
#               TOOCAN output files, etc).
#               File system scanning functions dedicated to IR satellite
#               products. Builds chronological file lists with timestamp
#               extraction based solely on filenames (no NetCDF opening)
#
# Provides    :
#   - extract_datetime_from_filename()
#   - fast_extract_datetime()
#   - scan_ir_files()
#   - build_ir_filelist()
#
# These functions do not open files; they rely solely on filename patterns,
# offering very high speed for large directory scans.
#
# License     : CNRS / LEGOS 
# =============================================================================


import re
from datetime import datetime
import os
import glob
import pandas as pd

def scan_ir_files(path_ir, extension=".nc"):
    filepaths = []
    for root, dirs, files in os.walk(path_ir):
        for f in files:
            if f.endswith(extension):
                filepaths.append(os.path.join(root, f))
    return sorted(filepaths)


def extract_datetime_from_filename(fname):
    """
    Extracts datetime from filenames using multiple patterns.
    Supports:
      - YYYYMMDDHHMM
      - YYYYMMDDHHMMSS
      - YYYY-MM-DDTHH-MM-SS
      - YYYYMMDD_HHMM
      - YYYYMMDDTHHMM
      - Anything containing those sequences
    """
    # 1) Compact format YYYYMMDDHHMM or YYYYMMDDHHMMSS
    match = re.search(r"(20\d{10})(\d{2})?", fname)
    if match:
        s = match.group(0)
        fmt = "%Y%m%d%H%M" if len(s) == 12 else "%Y%m%d%H%M%S"
        return datetime.strptime(s, fmt)

    # 2) MSG/EUMETSAT format YYYY-MM-DDTHH-MM-SS
    match = re.search(r"\d{4}-\d{2}-\d{2}T\d{2}-\d{2}-\d{2}", fname)
    if match:
        return datetime.strptime(match.group(0), "%Y-%m-%dT%H-%M-%S")

    # 3) Format YYYYMMDD_HHMM
    match = re.search(r"(20\d{6})_(\d{4})", fname)
    if match:
        s = match.group(1) + match.group(2)
        return datetime.strptime(s, "%Y%m%d%H%M")

    # 4) Format YYYYMMDDTHHMM or YYYYMMDDTHHMMSS
    match = re.search(r"20\d{6}T\d{4}(\d{2})?", fname)
    if match:
        s = match.group(0).replace("T", "")
        fmt = "%Y%m%d%H%M" if len(s) == 12 else "%Y%m%d%H%M%S"
        return datetime.strptime(s, fmt)

    raise ValueError(f"Cannot parse datetime from filename: {fname}")

import os
import glob
from datetime import datetime

def fast_extract_datetime_old(fname):
    """
    Efficient multi-pattern date extraction.
    Tries:
      - YYYYMMDDHHMM
      - YYYYMMDDHHMMSS
      - YYYYMMDDTHHMM
      - YYYYMMDDTHHMMSS
      - YYYY-MM-DDTHH-MM-SS
    (no regex: faster)
    """
    # Remove extension
    base = fname.split('.')[0]

    # 1) Look for YYYYMMDDHHMMSS (14 digits)
    for i in range(len(base)-13):
        chunk = base[i:i+14]
        if chunk.isdigit():
            try:
                return datetime.strptime(chunk, "%Y%m%d%H%M%S")
            except:
                pass

    # 2) Look for YYYYMMDDHHMM (12 digits)
    for i in range(len(base)-11):
        chunk = base[i:i+12]
        if chunk.isdigit():
            try:
                return datetime.strptime(chunk, "%Y%m%d%H%M")
            except:
                pass

    # 3) Look for YYYY-MM-DDTHH-MM-SS
    for i in range(len(base)-18):
        chunk = base[i:i+19]
        try:
            return datetime.strptime(chunk, "%Y-%m-%dT%H-%M-%S")
        except:
            pass

    # 4) Look for YYYYMMDDTHHMM or THHMMSS
    cleaned = base.replace("T", "")
    for i in range(len(cleaned)-11):
        chunk = cleaned[i:i+12]
        if chunk.isdigit():
            try:
                return datetime.strptime(chunk, "%Y%m%d%H%M")
            except:
                pass

    return None

def fast_extract_datetime(fname, dt_minutes=30):
    """
    Efficient multi-pattern date extraction.
    Tries:
      - YYYYMMDDHHMM
      - YYYYMMDDHHMMSS
      - YYYYMMDD_HHMM
      - YYYYMMDDTHHMM
      - YYYYMMDDTHHMMSS
      - YYYY-MM-DDTHH-MM-SS
      - TOOCAN slot format : YYYYMMDD-SLOT
    """
    base = fname.split('.')[0]

    # --------------------------------------------------------------
    # 0) NEW : TOOCAN format YYYYMMDD_HHMM
    # Example : ToocanCloudMask_20120801_1430.nc
    # --------------------------------------------------------------
    m = re.search(r"(\d{8})_(\d{4})$", base)
    if m:
        ymd = m.group(1)
        hm  = m.group(2)
        try:
            return datetime.strptime(ymd + hm, "%Y%m%d%H%M")
        except:
            pass

    # --------------------------------------------------------------
    # 1) Look for YYYYMMDDHHMMSS (14 digits)
    # --------------------------------------------------------------
    for i in range(len(base)-13):
        chunk = base[i:i+14]
        if chunk.isdigit():
            try:
                return datetime.strptime(chunk, "%Y%m%d%H%M%S")
            except:
                pass

    # --------------------------------------------------------------
    # 2) Look for YYYYMMDDHHMM (12 digits)
    # --------------------------------------------------------------
    for i in range(len(base)-11):
        chunk = base[i:i+12]
        if chunk.isdigit():
            try:
                return datetime.strptime(chunk, "%Y%m%d%H%M")
            except:
                pass

    # --------------------------------------------------------------
    # 3) Look for YYYY-MM-DDTHH-MM-SS
    # --------------------------------------------------------------
    for i in range(len(base)-18):
        chunk = base[i:i+19]
        try:
            return datetime.strptime(chunk, "%Y-%m-%dT%H-%M-%S")
        except:
            pass

    # --------------------------------------------------------------
    # 4) Look for YYYYMMDDTHHMM (T removed)
    # --------------------------------------------------------------
    cleaned = base.replace("T", "")
    for i in range(len(cleaned)-11):
        chunk = cleaned[i:i+12]
        if chunk.isdigit():
            try:
                return datetime.strptime(chunk, "%Y%m%d%H%M")
            except:
                pass

    # --------------------------------------------------------------
    # 5) TOOCAN slot format :  YYYYMMDD-SLOT
    # --------------------------------------------------------------
    #m = re.search(r"(\d{8})-(\d{1,3})$", base)
    #if m:
    #    ymd = m.group(1)
    ##    slot = int(m.group(2))
    ##    date0 = datetime.strptime(ymd, "%Y%m%d")
    #    return date0 + timedelta(minutes=(slot - 1) * dt_minutes)

    # --------------------------------------------------------------
    # Nothing matched
    # --------------------------------------------------------------
    return None


def build_ir_filelist(path_ir, start_time, end_time):
    """
    Fast file lister:
      - NO xarray (100x faster)
      - NO netcdf opening
      - date extracted from filename only
      - filters by period instantly
    """
    file_time_pairs = []

    # 1) Gather all candidate NetCDF files
    patterns = ["Mmultic3kmNC4_msg04_*.nc", "*.nc4", "*.nc.gz", "*.nc4.gz"]
    files = []
    for pat in patterns:
        files.extend(glob.glob(os.path.join(path_ir, "**", pat), recursive=True))

    print(f"Found {len(files)} files under {path_ir}")

    # 2) Loop over filenames ONLY
    for f in files:
        fname = os.path.basename(f)

        dt = fast_extract_datetime(fname)
        if dt is None:
            # Skip if no datetime found
            continue

        # 3) Period filtering (FAST)
        if dt < start_time or dt > end_time:
            continue

        # 4) Keep file
        file_time_pairs.append((f, dt))

    # 5) Sort once
    file_time_pairs.sort(key=lambda x: x[1])
    return file_time_pairs


def build_toocan_mask_filelist(path_toocan_mask, start_time, end_time):
    """
    Fast file lister:
      - NO xarray (100x faster)
      - NO netcdf opening
      - date extracted from filename only
      - filters by period instantly
    """
    file_time_pairs = []

    # 1) Gather all candidate NetCDF files
    patterns = ["*.nc", "*.nc4", "*.nc.gz", "*.nc4.gz"]
    files = []
    for pat in patterns:
        files.extend(glob.glob(os.path.join(path_toocan_mask, "**", pat), recursive=True))

    print(f"Found {len(files)} files under {path_toocan_mask}")

    # 2) Loop over filenames ONLY
    for f in files:
        fname = os.path.basename(f)
        

        dt = fast_extract_datetime(fname)
        if dt is None:
            # Skip if no datetime found
            continue

        # 3) Period filtering (FAST)
        if dt < start_time or dt > end_time:
            continue

        # 4) Keep file
        file_time_pairs.append((f, dt))

    # 5) Sort once
    file_time_pairs.sort(key=lambda x: x[1])
    #print(file_time_pairs)
    return file_time_pairs


def merge_ir_and_toocan(df_ir, df_mask, tolerance_minutes=1):
    """
    Merge IR filelist and TOOCAN mask filelist on datetime.

    Parameters
    ----------
    df_ir : DataFrame with columns ["datetime", "ir_path"]
    df_mask : DataFrame with columns ["datetime", "mask_path"]
    tolerance_minutes : allowed time difference in minutes (default = 1)

    Returns
    -------
    df_all : merged dataframe with columns :
         ["datetime", "ir_path", "mask_path"]
    """

    # Ensure datetime is datetime64
    df_ir = df_ir.copy()
    df_mask = df_mask.copy()

    df_ir["datetime"] = pd.to_datetime(df_ir["datetime"])
    df_mask["datetime"] = pd.to_datetime(df_mask["datetime"])

    # Sort both
    df_ir = df_ir.sort_values("datetime")
    df_mask = df_mask.sort_values("datetime")

    # As-of merge = match nearest timestamp
    df_all = pd.merge_asof(
        df_ir,
        df_mask,
        on="datetime",
        direction="nearest",
        tolerance=pd.Timedelta(minutes=tolerance_minutes)
    )

    # Identify rows where no mask matched (time gap)
    missing_mask = df_all["mask_path"].isna()
    if missing_mask.any():
        print(f"⚠ WARNING: {missing_mask} IR timestamps have no matching TOOCAN mask")

    # Identify rows where TOOCAN mask timestamps exist but no IR
    missing_ir = ~df_mask["datetime"].isin(df_all["datetime"])
    if missing_ir.any():
        print(f"⚠ WARNING: {missing_ir.sum()} TOOCAN masks not used in merged table")

    return df_all

