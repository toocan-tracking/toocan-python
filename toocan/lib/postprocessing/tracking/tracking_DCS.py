import pandas as pd
import numpy as np
import xarray as xr

import warnings
from xarray.conventions import SerializationWarning
warnings.filterwarnings("ignore", category=SerializationWarning)


from .tracking_utils import (
    compute_cold_surfaces,
    compute_temperature_stats,
    compute_equivalent_ellipse
)
from toocan.lib.io.open_IRdata import read_ir_file

def track_month_from_df(df_all,nav, surface_area_2d, lon_array_2d, lat_array_2d,
                        lonmin, lonmax, latmin, latmax,
                        model_name):

    tracks = {}

    LAT = lat_array_2d
    LON = lon_array_2d

    for _, row in df_all.iterrows():

        t = row["datetime"]
        print(f"[TRACK] {t}")
        print(row["ir_path"])
        # =============================
        # 1) LOAD MASK (skip if missing)
        # =============================
        if pd.isna(row["mask_path"]):
            print(f"⚠ skip {t} : NO MASK")
            continue

        try:
            with xr.open_dataset(row["mask_path"], decode_cf=False) as ds:
                raw = ds["DCS_number"].isel(time=0).values
                mask = raw.astype(np.int32)

            print("RAW unique:", np.unique(raw)[:10])
            print("MASK unique:", np.unique(mask)[:10])

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
            #ir = read_ir_file(
            #    row["ir_path"],
            #    lonmin, lonmax,
            #    latmin, latmax,
            #    model_name=model_name
            #)

            ir = read_ir_file(row["ir_path"], nav, model_name)


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
        print("Labels détectés :", labels)
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
            #print(st["tbmin"],st["tb90"],st["avg235"],st["avg208"])

            tr["tbmin"].append(st["tbmin"])
            tr["tb90"].append(st["tb90"])
            tr["avg235"].append(st["avg235"])
            tr["avg208"].append(st["avg208"])
            tr["avg200"].append(st["avg200"])

            # ---- Ellipses
            tr["ellipse_235"].append(
                compute_equivalent_ellipse(
                    bt2d, mask, lab, 235, lat_array_2d, lon_array_2d
                )
            )


            tr["ellipse_220"].append(
                compute_equivalent_ellipse(
                    bt2d, mask, lab, 220, lat_array_2d, lon_array_2d
                )
            )

    return tracks
