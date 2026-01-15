import xarray as xr
import numpy as np
import pandas as pd
import glob
from datetime import datetime
import numpy as np
from datetime import datetime
import re
import sys 
import os

sys.path.insert(0, os.path.expanduser("~/TOOCAN/pyTOOCAN/src"))

from toocan.utils.parse_param_file import parse_param_file
from toocan.io.open_navigationFile import crop_navigation_file
from toocan.io.open_IRdata import read_ir_file
from toocan.struct.data_param import DataParam

import warnings
from xarray.conventions import SerializationWarning
warnings.filterwarnings("ignore", category=SerializationWarning)



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


def compute_cold_surfaces(bt, mask, lab, surface_area_2d):
    """
    Renvoie surfaces froides pour 235/220/210/200 K (en pixels et km²)
    pour un DCS donné.
    """
    m = (mask == lab)

    out = {}

    for thr in [235, 220, 210, 200]:
        cold = m & (bt <= thr)

        pix = int(cold.sum())
        km2 = float(surface_area_2d[cold].sum())

        out[f"pix_{thr}K"] = pix
        out[f"km2_{thr}K"] = km2

    return out

def compute_temperature_stats(bt, mask, lab):
    m = (mask == lab)
    if not np.any(m):
        return {"tbmin": np.nan, "tb90": np.nan, "avg235": np.nan, "avg208": np.nan, "avg200": np.nan}

    vals = bt[m]

    tb_min  = float(np.nanmin(vals))
    tb_90   = float(np.nanpercentile(vals, 90))

    avg235 = float(vals[vals <= 235].mean()) if np.any(vals <= 235) else np.nan
    avg208 = float(vals[vals <= 208].mean()) if np.any(vals <= 208) else np.nan
    avg200 = float(vals[vals <= 200].mean()) if np.any(vals <= 200) else np.nan

    return {
        "tbmin": tb_min,
        "tb90":  tb_90,
        "avg235": avg235,
        "avg208": avg208,
        "avg200": avg200,
    }



def compute_equivalent_ellipse(bt, mask, lab, thr, lat_array_2d, lon_array_2d):
    """
    Calcule ellipse équivalente à un seuil BT.
    thr : 235 ou 220 K

    Paramètres :
      bt  = champ IR (2D)
      mask = masque TOOCAN
      lab = label du DCS
      thr = seuil froid
      lat_array_2d, lon_array_2d : grilles de lat/lon

    Retourne (a_km, b_km, eccentricite, orientation_deg)
    """

    m = (mask == lab)
    cold = m & (bt <= thr)

    yy, xx = np.nonzero(cold)
    if yy.size < 5:
        return (np.nan, np.nan, np.nan, np.nan)

    # -------- Centre de masse --------
    y0 = yy.mean()
    x0 = xx.mean()

    # -------- Coordonnées en km --------
    lat0 = lat_array_2d[int(y0), int(x0)]
    lon0 = lon_array_2d[int(y0), int(x0)]

    # Conversion deg → km approx
    R = 6371.0
    km_per_deg_lat = 111.0
    km_per_deg_lon = 111.0 * np.cos(np.deg2rad(lat0))

    x_km = (lon_array_2d[yy, xx] - lon0) * km_per_deg_lon
    y_km = (lat_array_2d[yy, xx] - lat0) * km_per_deg_lat

    # -------- Covariance --------
    coords = np.vstack([x_km, y_km])
    cov = np.cov(coords)

    eigvals, eigvecs = np.linalg.eig(cov)

    # Semi-axes
    a = 2 * np.sqrt(eigvals.max())
    b = 2 * np.sqrt(eigvals.min())

    vec = eigvecs[:, np.argmax(eigvals)]
    angle = np.degrees(np.arctan2(vec[1], vec[0]))

    ecc = b / a

    return float(a), float(b), float(ecc), float(angle)


def parse_time_from_filename(filepath):
    import re
    name = filepath.split("/")[-1]
    match = re.search(r"(\d{8})_(\d{4})", name)
    if not match:
        raise ValueError(f"Bad filename format: {name}")
    ymd, hm = match.groups()
    return datetime.strptime(ymd + hm, "%Y%m%d%H%M")


def list_files_with_times(path_pattern):
    files = sorted(glob(path_pattern))
    ft = [(f, parse_time_from_filename(f)) for f in files]
    ft.sort(key=lambda x: x[1])
    return ft


def iter_masks(file_times, expected_dt_minutes=30):
    prev_time = None
    dt_expected = pd.Timedelta(minutes=expected_dt_minutes)

    for filepath, t in file_times:

        if prev_time is not None:
            gap = t - prev_time
            if abs(gap - dt_expected) > pd.Timedelta(minutes=1):
                print(f"⚠ WARNING time gap: {prev_time} → {t} (gap={gap})")

        prev_time = t

        with xr.open_dataset(filepath) as ds:
            mask = ds["DCS_number"].isel(time=0).values.astype(np.int32)

        yield mask, t


import xarray as xr
import numpy as np
import pandas as pd

def track_month_from_df(df_all, surface_area_2d, lon_array_2d, lat_array_2d,
                        lonmin, lonmax, latmin, latmax,
                        model_name):

    tracks = {}

    LAT = lat_array_2d
    LON = lon_array_2d

    for _, row in df_all.iterrows():

        t = row["datetime"]
        print(f"[TRACK] {t}")

        # =============================
        # 1) LOAD MASK (skip if missing)
        # =============================
        if pd.isna(row["mask_path"]):
            print(f"⚠ skip {t} : NO MASK")
            continue

        try:
            with xr.open_dataset(row["mask_path"]) as ds_mask:
                mask = ds_mask["DCS_number"].isel(time=0).values.astype(np.int32)
        except Exception as e:
            print(f"⚠ skip {t} : ERROR reading mask → {e}")
            continue

        # =============================
        # 2) LOAD IR BT (skip if missing)
        # =============================
        if pd.isna(row["ir_path"]):
            print(f"⚠ skip {t} : NO IR data")
            continue

        try:
            ir = read_ir_file(
                row["ir_path"],
                lonmin, lonmax,
                latmin, latmax,
                model_name=model_name
            )
            bt = ir["irbt"]
            bt2d = bt if bt.ndim == 2 else bt[0]
        except Exception as e:
            print(f"⚠ skip {t} : ERROR reading IR BT → {e}")
            continue

        # =============================
        # 3) PROCESS ALL LABELS
        # =============================
        labels = np.unique(mask)
        labels = labels[labels > 0]

        for lab in labels:

            m = (mask == lab)
            if not np.any(m):
                continue

            yy, xx = np.nonzero(m)
            y_cm = int(yy.mean())
            x_cm = int(xx.mean())

            lat_cm = float(LAT[y_cm, x_cm])
            lon_cm = float(LON[y_cm, x_cm])

            # Create track if new
            if lab not in tracks:
                tracks[lab] = {
                    "times": [],
                    "area": [],
                    "lon_cm": [],
                    "lat_cm": [],
                    "tbmin": [],
                    "tb90": [],
                    "avg235": [],
                    "avg208": [],
                    "avg200": [],
                    "pix_235": [],
                    "pix_220": [],
                    "pix_210": [],
                    "pix_200": [],
                    "km2_235": [],
                    "km2_220": [],
                    "km2_210": [],
                    "km2_200": [],
                    "ellipse_235": [],
                    "ellipse_220": [],
                }

            tr = tracks[lab]

            # ---- Geometry
            tr["times"].append(t)
            tr["area"].append(int(m.sum()))
            tr["lon_cm"].append(lon_cm)
            tr["lat_cm"].append(lat_cm)

            # ---- Cold surfaces
            cold = compute_cold_surfaces(bt2d, mask, lab, surface_area_2d)

            tr["pix_235"].append(cold["pix_235K"])
            tr["pix_220"].append(cold["pix_220K"])
            tr["pix_210"].append(cold["pix_210K"])
            tr["pix_200"].append(cold["pix_200K"])

            tr["km2_235"].append(cold["km2_235K"])
            tr["km2_220"].append(cold["km2_220K"])
            tr["km2_210"].append(cold["km2_210K"])
            tr["km2_200"].append(cold["km2_200K"])

            # ---- Temperature stats
            st = compute_temperature_stats(bt2d, mask, lab)
            tr["tbmin"].append(st["tbmin"])
            tr["tb90"].append(st["tb90"])
            tr["avg235"].append(st["avg235"])
            tr["avg208"].append(st["avg208"])
            tr["avg200"].append(st["avg200"])

            # ---- Ellipses
            tr["ellipse_235"].append(
                compute_equivalent_ellipse(
                    bt2d, mask, lab, thr=235,
                    lat_array_2d=lat_array_2d,
                    lon_array_2d=lon_array_2d
                )
            )

            tr["ellipse_220"].append(
                compute_equivalent_ellipse(
                    bt2d, mask, lab, thr=220,
                    lat_array_2d=lat_array_2d,
                    lon_array_2d=lon_array_2d
                )
            )

    return tracks

def compute_INT_variables(tracks, dt_minutes=30):
    """
    Calcule toutes les INT_* principales :
      - durée
      - UTC_init / UTC_end
      - lon/lat init & end
      - lon/lat min/max
      - distance, vitesse moyenne
      - tbmin globale
      - surfmax km2 aux différents seuils
      - surfcum km2 à 235K
    """
    INT = {}
    dt_hours = dt_minutes / 60.0
    R = 6370.0  # rayon Terre (km)

    for lab, tr in tracks.items():
        times   = tr["times"]
        lon_cm  = np.array(tr["lon_cm"],  float)
        lat_cm  = np.array(tr["lat_cm"],  float)
        tb_min  = np.array(tr["tb_min"],  float)
        s235    = np.array(tr["surfKm2_235"], float)
        s220    = np.array(tr["surfKm2_220"], float)
        s210    = np.array(tr["surfKm2_210"], float)
        s200    = np.array(tr["surfKm2_200"], float)

        t_init = times[0]
        t_end  = times[-1]

        lon_init, lat_init = lon_cm[0],  lat_cm[0]
        lon_end,  lat_end  = lon_cm[-1], lat_cm[-1]

        duration_hours = len(times) * dt_hours

        # distance (Haversine entre init et fin)
        dlon = np.deg2rad(lon_end - lon_init)
        dlat = np.deg2rad(lat_end - lat_init)
        a = (np.sin(dlat/2)**2 +
             np.cos(np.deg2rad(lat_init)) *
             np.cos(np.deg2rad(lat_end)) *
             np.sin(dlon/2)**2)
        dist_km = 2 * R * np.sqrt(a)
        vel_ms = (dist_km * 1000.0) / (duration_hours*3600.0) if duration_hours > 0 else np.nan

        INT[lab] = {
            # ID de base
            "DCS_number": int(lab),

            # durée
            "duration_hours": float(duration_hours),

            # temps UTC
            "UTC_init": int(t_init.timestamp()),
            "UTC_end":  int(t_end.timestamp()),

            # coordonnées init / fin
            "lon_init": float(lon_init),
            "lat_init": float(lat_init),
            "lon_end":  float(lon_end),
            "lat_end":  float(lat_end),

            # min/max lon/lat
            "lon_min": float(np.nanmin(lon_cm)),
            "lon_max": float(np.nanmax(lon_cm)),
            "lat_min": float(np.nanmin(lat_cm)),
            "lat_max": float(np.nanmax(lat_cm)),

            # TB globale
            "tb_min_global": float(np.nanmin(tb_min)),

            # surfaces max
            "surfmaxkm2_235": float(np.nanmax(s235)),
            "surfmaxkm2_220": float(np.nanmax(s220)),
            "surfmaxkm2_210": float(np.nanmax(s210)),
            "surfmaxkm2_200": float(np.nanmax(s200)),

            # surface cumulée 235K
            "surfcumkm2_235": float(np.nansum(s235)),

            # cinématique
            "distance_km": float(dist_km),
            "velocity_ms": float(vel_ms),

            # place-holders pour QC et classifications (à remplir plus tard)
            "qltyDCS": 11100,   # ex : tout OK, 0 image interpolée
            "classif_Roca": -999,
            "classif_JIRAK": -999,
            "classif_MADDOX": -999,
            "TS_number_IBTRACS": -999,
            "TS_nature_IBTRACS": -999,
            "TS_mindistance_IBTRACS": -999.0,
        }

    return INT

def build_LC_variables(tracks, global_times):
    """
    Construit tous les LC_* principaux / temps.

    global_times : liste triée de tous les timestamps (pandas.Timestamp)
    """

    dcs_ids = sorted(tracks.keys())
    n_dcs   = len(dcs_ids)
    n_time  = len(global_times)

    # index temporel
    index_time = {t: i for i, t in enumerate(global_times)}

    # Création des tableaux
    fill_int  = -999
    fill_float = -999.0

    LC_tbmin      = np.full((n_dcs, n_time), fill_float, dtype=np.float32)
    LC_tb90th     = np.full((n_dcs, n_time), fill_float, dtype=np.float32)
    LC_tbavg_235K = np.full((n_dcs, n_time), fill_float, dtype=np.float32)
    LC_tbavg_208K = np.full((n_dcs, n_time), fill_float, dtype=np.float32)
    LC_tbavg_200K = np.full((n_dcs, n_time), fill_float, dtype=np.float32)

    LC_lon  = np.full((n_dcs, n_time), fill_float, dtype=np.float32)
    LC_lat  = np.full((n_dcs, n_time), fill_float, dtype=np.float32)
    LC_x    = np.full((n_dcs, n_time), fill_int,   dtype=np.int32)
    LC_y    = np.full((n_dcs, n_time), fill_int,   dtype=np.int32)
    LC_UTC  = np.full((n_dcs, n_time), fill_int,   dtype=np.int32)
    LC_vel  = np.full((n_dcs, n_time), fill_float, dtype=np.float32)

    LC_surfPix_235K = np.full((n_dcs, n_time), fill_int,   dtype=np.int32)
    LC_surfPix_210K = np.full((n_dcs, n_time), fill_int,   dtype=np.int32)
    LC_surfkm2_235K = np.full((n_dcs, n_time), fill_float, dtype=np.float32)
    LC_surfkm2_220K = np.full((n_dcs, n_time), fill_float, dtype=np.float32)
    LC_surfkm2_210K = np.full((n_dcs, n_time), fill_float, dtype=np.float32)
    LC_surfkm2_200K = np.full((n_dcs, n_time), fill_float, dtype=np.float32)

    # Ellipses (on ne les calcule pas vraiment ici, mais on prépare les matrices)
    LC_semiminor_235K   = np.full((n_dcs, n_time), fill_float, dtype=np.float32)
    LC_semimajor_235K   = np.full((n_dcs, n_time), fill_float, dtype=np.float32)
    LC_ecc_235K         = np.full((n_dcs, n_time), fill_float, dtype=np.float32)
    LC_orientation_235K = np.full((n_dcs, n_time), fill_float, dtype=np.float32)

    LC_semiminor_220K   = np.full((n_dcs, n_time), fill_float, dtype=np.float32)
    LC_semimajor_220K   = np.full((n_dcs, n_time), fill_float, dtype=np.float32)
    LC_ecc_220K         = np.full((n_dcs, n_time), fill_float, dtype=np.float32)
    LC_orientation_220K = np.full((n_dcs, n_time), fill_float, dtype=np.float32)

    # Remplissage
    for i, lab in enumerate(dcs_ids):
        tr = tracks[lab]

        # pour vitesse instantanée : on gardera (lon,lat,t) pour ce DCS
        prev_lon = prev_lat = prev_t = None

        for t, area_pix, lon_cm, lat_cm, x_cm, y_cm, tb_min, tb90, ta235, ta208, ta200, \
            s235pix, s235km2, s220km2, s210km2, s200km2 in zip(
                tr["times"],
                tr["area_pix"],
                tr["lon_cm"], tr["lat_cm"],
                tr["x_cm"], tr["y_cm"],
                tr["tb_min"], tr["tb90"],
                tr["tbavg_235"], tr["tbavg_208"], tr["tbavg_200"],
                tr["surfPix_235"], tr["surfKm2_235"],
                tr["surfKm2_220"], tr["surfKm2_210"], tr["surfKm2_200"]
            ):

            j = index_time[t]

            LC_tbmin[i, j]      = tb_min
            LC_tb90th[i, j]     = tb90
            LC_tbavg_235K[i, j] = ta235
            LC_tbavg_208K[i, j] = ta208
            LC_tbavg_200K[i, j] = ta200

            LC_lon[i, j] = lon_cm
            LC_lat[i, j] = lat_cm
            LC_x[i, j]   = x_cm
            LC_y[i, j]   = y_cm
            LC_UTC[i, j] = int(t.timestamp())

            LC_surfPix_235K[i, j] = s235pix
            # NB : pas de LC_surfPix_220K dans l’original, seulement 235 & 210
            LC_surfPix_210K[i, j] = fill_int  # placeholder si tu ne le calcules pas
            LC_surfkm2_235K[i, j] = s235km2
            LC_surfkm2_220K[i, j] = s220km2
            LC_surfkm2_210K[i, j] = s210km2
            LC_surfkm2_200K[i, j] = s200km2

            # vitesse instantanée (m/s) le long de la trajectoire
            if prev_lon is not None:
                R = 6370.0
                dlon = np.deg2rad(lon_cm - prev_lon)
                dlat = np.deg2rad(lat_cm - prev_lat)
                a = (np.sin(dlat/2)**2 +
                     np.cos(np.deg2rad(prev_lat)) *
                     np.cos(np.deg2rad(lat_cm)) *
                     np.sin(dlon/2)**2)
                dist_km = 2 * R * np.sqrt(a)
                dt_sec  = (t - prev_t).total_seconds()
                LC_vel[i, j] = (dist_km*1000.0)/dt_sec if dt_sec > 0 else fill_float

            prev_lon, prev_lat, prev_t = lon_cm, lat_cm, t

        # ⚠ les ellipses LC_*_235K / 220K restent à remplir si tu veux
        #   calculer ellipse équivalente à chaque pas de temps. Ici on les laisse à -999.

    return {
        "DCS_ids": dcs_ids,
        "LC_tbmin": LC_tbmin,
        "LC_tb90th": LC_tb90th,
        "LC_tbavg_235K": LC_tbavg_235K,
        "LC_tbavg_208K": LC_tbavg_208K,
        "LC_tbavg_200K": LC_tbavg_200K,
        "LC_lon": LC_lon,
        "LC_lat": LC_lat,
        "LC_x": LC_x,
        "LC_y": LC_y,
        "LC_UTC": LC_UTC,
        "LC_velocity": LC_vel,
        "LC_surfPix_235K": LC_surfPix_235K,
        "LC_surfPix_210K": LC_surfPix_210K,
        "LC_surfkm2_235K": LC_surfkm2_235K,
        "LC_surfkm2_220K": LC_surfkm2_220K,
        "LC_surfkm2_210K": LC_surfkm2_210K,
        "LC_surfkm2_200K": LC_surfkm2_200K,
        "LC_semiminor_235K": LC_semiminor_235K,
        "LC_semimajor_235K": LC_semimajor_235K,
        "LC_ecc_235K": LC_ecc_235K,
        "LC_orientation_235K": LC_orientation_235K,
        "LC_semiminor_220K": LC_semiminor_220K,
        "LC_semimajor_220K": LC_semimajor_220K,
        "LC_ecc_220K": LC_ecc_220K,
        "LC_orientation_220K": LC_orientation_220K,
    }

import xarray as xr
import numpy as np

def write_tracking_nc(outfile, INT, LC, global_times,
                      title="TOOCAN - Morphological characteristics of the Deep Convective Systems",
                      region="AFRICA",
                      platform="MSG",
                      version="2.08"):
    """
    Écrit un NetCDF proche de TOOCAN v2.08.
    """

    dcs_ids = LC["DCS_ids"]
    n_dcs   = len(dcs_ids)
    n_time  = len(global_times)

    # vecteurs INT
    INT_duration     = [INT[lab]["duration_hours"]      for lab in dcs_ids]
    INT_UTC_timeInit = [INT[lab]["UTC_init"]            for lab in dcs_ids]
    INT_UTC_timeEnd  = [INT[lab]["UTC_end"]             for lab in dcs_ids]
    INT_lonInit      = [INT[lab]["lon_init"]            for lab in dcs_ids]
    INT_latInit      = [INT[lab]["lat_init"]            for lab in dcs_ids]
    INT_lonEnd       = [INT[lab]["lon_end"]             for lab in dcs_ids]
    INT_latEnd       = [INT[lab]["lat_end"]             for lab in dcs_ids]
    INT_velocityAvg  = [INT[lab]["velocity_ms"]         for lab in dcs_ids]
    INT_distance     = [INT[lab]["distance_km"]         for lab in dcs_ids]
    INT_lonmin       = [INT[lab]["lon_min"]             for lab in dcs_ids]
    INT_lonmax       = [INT[lab]["lon_max"]             for lab in dcs_ids]
    INT_latmin       = [INT[lab]["lat_min"]             for lab in dcs_ids]
    INT_latmax       = [INT[lab]["lat_max"]             for lab in dcs_ids]
    INT_tbmin        = [INT[lab]["tb_min_global"]       for lab in dcs_ids]
    INT_surfmaxkm2_235 = [INT[lab]["surfmaxkm2_235"]    for lab in dcs_ids]
    INT_surfmaxkm2_220 = [INT[lab]["surfmaxkm2_220"]    for lab in dcs_ids]
    INT_surfmaxkm2_210 = [INT[lab]["surfmaxkm2_210"]    for lab in dcs_ids]
    INT_surfmaxkm2_200 = [INT[lab]["surfmaxkm2_200"]    for lab in dcs_ids]
    INT_surfcumkm2_235 = [INT[lab]["surfcumkm2_235"]    for lab in dcs_ids]
    INT_qltyDCS        = [INT[lab]["qltyDCS"]           for lab in dcs_ids]
    INT_classif        = [INT[lab]["classif_Roca"]      for lab in dcs_ids]
    INT_classif_JIRAK  = [INT[lab]["classif_JIRAK"]     for lab in dcs_ids]
    INT_classif_MADDOX = [INT[lab]["classif_MADDOX"]    for lab in dcs_ids]
    INT_TS_number_IBTRACS   = [INT[lab]["TS_number_IBTRACS"]   for lab in dcs_ids]
    INT_TS_nature_IBTRACS   = [INT[lab]["TS_nature_IBTRACS"]   for lab in dcs_ids]
    INT_TS_mindistance_IBTRACS = [INT[lab]["TS_mindistance_IBTRACS"] for lab in dcs_ids]

    # temps en int
    time_seconds = np.array([np.int32(t.timestamp()) for t in global_times], dtype=np.int32)

    # QC GEO IR image par pas de temps -> ici 1 (= full data OK)
    QCgeo_IRimage = np.ones(n_time, dtype=np.int32)

    ds = xr.Dataset(
        coords = {
            "DCS": ("DCS", np.array(dcs_ids, dtype=np.int32)),
            "time": ("time", time_seconds),
        },
        data_vars = {
            # === Intégral (DCS) ===
            "DCS": (("DCS",), np.array(dcs_ids, dtype=np.int32)),
            "INT_DCSnumber": (("DCS",), np.array(dcs_ids, dtype=np.int32)),
            "INT_DCS_qualitycontrol": (("DCS",), np.array(INT_qltyDCS, dtype=np.int32)),
            "INT_classif": (("DCS",), np.array(INT_classif, dtype=np.int32)),
            "INT_duration": (("DCS",), np.array(INT_duration, dtype=np.float32)),
            "INT_UTC_timeInit": (("DCS",), np.array(INT_UTC_timeInit, dtype=np.int32)),
            "INT_UTC_timeEnd":  (("DCS",), np.array(INT_UTC_timeEnd,  dtype=np.int32)),
            "INT_lonInit": (("DCS",), np.array(INT_lonInit, dtype=np.float32)),
            "INT_latInit": (("DCS",), np.array(INT_latInit, dtype=np.float32)),
            "INT_lonEnd":  (("DCS",), np.array(INT_lonEnd,  dtype=np.float32)),
            "INT_latEnd":  (("DCS",), np.array(INT_latEnd,  dtype=np.float32)),
            "INT_velocityAvg": (("DCS",), np.array(INT_velocityAvg, dtype=np.float32)),
            "INT_distance":    (("DCS",), np.array(INT_distance,    dtype=np.float32)),
            "INT_lonmin": (("DCS",), np.array(INT_lonmin, dtype=np.float32)),
            "INT_lonmax": (("DCS",), np.array(INT_lonmax, dtype=np.float32)),
            "INT_latmin": (("DCS",), np.array(INT_latmin, dtype=np.float32)),
            "INT_latmax": (("DCS",), np.array(INT_latmax, dtype=np.float32)),
            "INT_tbmin":  (("DCS",), np.array(INT_tbmin,  dtype=np.float32)),
            "INT_surfmaxkm2_235K": (("DCS",), np.array(INT_surfmaxkm2_235, dtype=np.float32)),
            "INT_surfmaxkm2_220K": (("DCS",), np.array(INT_surfmaxkm2_220, dtype=np.float32)),
            "INT_surfmaxkm2_210K": (("DCS",), np.array(INT_surfmaxkm2_210, dtype=np.float32)),
            "INT_surfmaxkm2_200K": (("DCS",), np.array(INT_surfmaxkm2_200, dtype=np.float32)),
            "INT_surfcumkm2_235K": (("DCS",), np.array(INT_surfcumkm2_235, dtype=np.float32)),
            "INT_classif_JIRAK":   (("DCS",), np.array(INT_classif_JIRAK,  dtype=np.int32)),
            "INT_classif_MADDOX":  (("DCS",), np.array(INT_classif_MADDOX, dtype=np.int32)),
            "INT_TS_number_IBTRACS":   (("DCS",), np.array(INT_TS_number_IBTRACS, dtype=np.int32)),
            "INT_TS_nature_IBTRACS":   (("DCS",), np.array(INT_TS_nature_IBTRACS, dtype=np.int32)),
            "INT_TS_mindistance_IBTRACS": (("DCS",), np.array(INT_TS_mindistance_IBTRACS, dtype=np.float32)),

            # === QC GEO ===
            "QCgeo_IRimage": (("time",), QCgeo_IRimage),

            # === LC : températures ===
            "LC_tbmin":      (("DCS","time"), LC["LC_tbmin"]),
            "LC_tbavg_235K": (("DCS","time"), LC["LC_tbavg_235K"]),
            "LC_tbavg_208K": (("DCS","time"), LC["LC_tbavg_208K"]),
            "LC_tbavg_200K": (("DCS","time"), LC["LC_tbavg_200K"]),
            "LC_tb90th":     (("DCS","time"), LC["LC_tb90th"]),

            # === LC : temps, positions ===
            "LC_UTC_time": (("DCS","time"), LC["LC_UTC"]),
            "LC_lon":      (("DCS","time"), LC["LC_lon"]),
            "LC_lat":      (("DCS","time"), LC["LC_lat"]),
            "LC_x":        (("DCS","time"), LC["LC_x"]),
            "LC_y":        (("DCS","time"), LC["LC_y"]),
            "LC_velocity": (("DCS","time"), LC["LC_velocity"]),

            # === LC : surfaces froides ===
            "LC_surfPix_235K": (("DCS","time"), LC["LC_surfPix_235K"]),
            "LC_surfPix_210K": (("DCS","time"), LC["LC_surfPix_210K"]),
            "LC_surfkm2_235K": (("DCS","time"), LC["LC_surfkm2_235K"]),
            "LC_surfkm2_220K": (("DCS","time"), LC["LC_surfkm2_220K"]),
            "LC_surfkm2_210K": (("DCS","time"), LC["LC_surfkm2_210K"]),
            "LC_surfkm2_200K": (("DCS","time"), LC["LC_surfkm2_200K"]),

            # === LC : ellipses équivalentes (placeholders pour le moment) ===
            "LC_semiminor_235K":   (("DCS","time"), LC["LC_semiminor_235K"]),
            "LC_semimajor_235K":   (("DCS","time"), LC["LC_semimajor_235K"]),
            "LC_ecc_235K":         (("DCS","time"), LC["LC_ecc_235K"]),
            "LC_orientation_235K": (("DCS","time"), LC["LC_orientation_235K"]),
            "LC_semiminor_220K":   (("DCS","time"), LC["LC_semiminor_220K"]),
            "LC_semimajor_220K":   (("DCS","time"), LC["LC_semimajor_220K"]),
            "LC_ecc_220K":         (("DCS","time"), LC["LC_ecc_220K"]),
            "LC_orientation_220K": (("DCS","time"), LC["LC_orientation_220K"]),
        }
    )

    # Attributs minimaux
    ds["time"].attrs.update({
        "units": "seconds since 1970-01-01",
        "long_name": "time",
    })

    ds.attrs.update({
        "title": title,
        "creator_name": "Thomas Fiolleau",
        "institution": "CNRS/LEGOS/IPSL",
        "conventions": "CF-1.6, ACDD-1.3",
        "tracker": "TOOCAN",
        "version": version,
        "Geostationary_platform": platform,
        "region": region,
    })

    ds.to_netcdf(outfile)
    print("✔ Tracking saved →", outfile)


import re
from datetime import datetime, timedelta
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
    m = re.search(r"(\d{8})-(\d{1,3})$", base)
    if m:
        ymd = m.group(1)
        slot = int(m.group(2))
        date0 = datetime.strptime(ymd, "%Y%m%d")
        return date0 + timedelta(minutes=(slot - 1) * dt_minutes)

    # --------------------------------------------------------------
    # Nothing matched
    # --------------------------------------------------------------
    return None


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
    patterns = ["*.nc", "*.nc4", "*.nc.gz", "*.nc4.gz"]
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

import pandas as pd
import numpy as np

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



path_masks = data_param.path_out.decode()+"/toocan_2.08/"
path_ir =params_GEO["path_ir"]+ "/2012/*/*"



fileIR_time = build_ir_filelist(params_GEO["path_ir"], start_time, end_time)
df_ir = pd.DataFrame(fileIR_time, columns=["ir_path", "datetime"])
maskTOOCAN_time = build_toocan_mask_filelist(path_masks, start_time, end_time)
df_mask = pd.DataFrame(maskTOOCAN_time, columns=["mask_path", "datetime"])
df_all = merge_ir_and_toocan(df_ir, df_mask)

# 1) tracking complet
tracks = track_month_from_df(df_all, surface_area_2d,lon_array_2d,lat_array_2d, lonmin, lonmax, latmin, latmax,model_name=model_name)

# 2) INT
INT = compute_INT_variables(tracks, dt_minutes=temporalresolution)

# 3) grille temporelle
global_times = sorted({t for tr in tracks.values() for t in tr["times"]})

# 4) LC
LC = build_LC_variables(tracks, global_times)

# 5) écriture NetCDF
outfile = "TOOCAN-AFRICA-20120801-20120831.nc"  # à construire depuis start/end
write_tracking_nc(outfile, INT, LC, global_times,
                  region="AFRICA",
                  platform=params_GEO["GEOplatform"],
                  version=params_TOOCAN["version"])