import numpy as np

def compute_INT_variables(tracks, dt_minutes=30):

    INT = {}
    R = 6371.0
    dt_hours = dt_minutes / 60.0

    for lab, tr in tracks.items():
        times = tr["times"]
        lon = np.array(tr["lon_cm"])
        lat = np.array(tr["lat_cm"])

        # Duration
        duration_hours = len(times) * dt_hours

        # Distance (init → end)
        lon0, lat0 = lon[0], lat[0]
        lon1, lat1 = lon[-1], lat[-1]

        dlon = np.deg2rad(lon1 - lon0)
        dlat = np.deg2rad(lat1 - lat0)

        a = (np.sin(dlat/2)**2 +
             np.cos(np.deg2rad(lat0))*np.cos(np.deg2rad(lat1))*np.sin(dlon/2)**2)

        dist_km = 2 * R * np.sqrt(a)
        vel_ms = (dist_km*1000)/(duration_hours*3600) if duration_hours>0 else np.nan

        INT[lab] = {
            "DCS_number": int(lab),
            "duration_hours": duration_hours,
            "UTC_init": int(times[0].timestamp()),
            "UTC_end":  int(times[-1].timestamp()),
            "lon_init": float(lon0),
            "lat_init": float(lat0),
            "lon_end":  float(lon1),
            "lat_end":  float(lat1),
            "lon_min": float(np.nanmin(lon)),
            "lon_max": float(np.nanmax(lon)),
            "lat_min": float(np.nanmin(lat)),
            "lat_max": float(np.nanmax(lat)),
            "tb_min_global": float(np.nanmin(tr["tbmin"])),
            "surfmaxkm2_235": float(np.nanmax(tr["km2_235"])),
            "surfmaxkm2_220": float(np.nanmax(tr["km2_220"])),
            "surfmaxkm2_210": float(np.nanmax(tr["km2_210"])),
            "surfmaxkm2_200": float(np.nanmax(tr["km2_200"])),
            "surfcumkm2_235": float(np.nansum(tr["km2_235"])),
            "distance_km": float(dist_km),
            "velocity_ms": float(vel_ms),
            "qltyDCS": 11100,
            "classif_Roca": -999,
            "classif_JIRAK": -999,
            "classif_MADDOX": -999,
            "TS_number_IBTRACS": -999,
            "TS_nature_IBTRACS": -999,
            "TS_mindistance_IBTRACS": -999.0,
        }

    return INT