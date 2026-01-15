import numpy as np

def build_LC_variables(tracks, global_times):

    dcs_ids = sorted(tracks.keys())
    n_dcs = len(dcs_ids)
    n_time = len(global_times)

    # Index time
    index_time = {t: i for i, t in enumerate(global_times)}

    fill_int = -999
    fill_float = -999.0

    LC = {
        "DCS_ids": dcs_ids,
        "LC_tbmin": np.full((n_dcs, n_time), fill_float, dtype=np.float32),
        "LC_tb90th": np.full((n_dcs, n_time), fill_float, dtype=np.float32),
        "LC_tbavg_235K": np.full((n_dcs, n_time), fill_float, dtype=np.float32),
        "LC_tbavg_208K": np.full((n_dcs, n_time), fill_float, dtype=np.float32),
        "LC_tbavg_200K": np.full((n_dcs, n_time), fill_float, dtype=np.float32),
        "LC_lon": np.full((n_dcs, n_time), fill_float, dtype=np.float32),
        "LC_lat": np.full((n_dcs, n_time), fill_float, dtype=np.float32),
        "LC_UTC": np.full((n_dcs, n_time), fill_int, dtype=np.int32),
        "LC_velocity": np.full((n_dcs, n_time), fill_float, dtype=np.float32),
        "LC_surfPix_235K": np.full((n_dcs, n_time), fill_int, dtype=np.int32),
        "LC_surfkm2_235K": np.full((n_dcs, n_time), fill_float, dtype=np.float32),
        "LC_surfkm2_220K": np.full((n_dcs, n_time), fill_float, dtype=np.float32),
        "LC_surfkm2_210K": np.full((n_dcs, n_time), fill_float, dtype=np.float32),
        "LC_surfkm2_200K": np.full((n_dcs, n_time), fill_float, dtype=np.float32),
        # Ellipse placeholders
        "LC_semiminor_235K": np.full((n_dcs, n_time), fill_float),
        "LC_semimajor_235K": np.full((n_dcs, n_time), fill_float),
        "LC_ecc_235K": np.full((n_dcs, n_time), fill_float),
        "LC_orientation_235K": np.full((n_dcs, n_time), fill_float),
        "LC_semiminor_220K": np.full((n_dcs, n_time), fill_float),
        "LC_semimajor_220K": np.full((n_dcs, n_time), fill_float),
        "LC_ecc_220K": np.full((n_dcs, n_time), fill_float),
        "LC_orientation_220K": np.full((n_dcs, n_time), fill_float),
    }

    # Fill each DCS
    for i, lab in enumerate(dcs_ids):
        tr = tracks[lab]
        times = tr["times"]
        lon = tr["lon_cm"]
        lat = tr["lat_cm"]

        prev_lon = prev_lat = prev_t = None

        for idx in range(len(times)):
            t = times[idx]
            j = index_time[t]

            LC["LC_lon"][i, j] = lon[idx]
            LC["LC_lat"][i, j] = lat[idx]
            LC["LC_tbmin"][i, j] = tr["tbmin"][idx]
            LC["LC_tb90th"][i, j] = tr["tb90"][idx]
            LC["LC_tbavg_235K"][i, j] = tr["avg235"][idx]
            LC["LC_tbavg_208K"][i, j] = tr["avg208"][idx]
            LC["LC_tbavg_200K"][i, j] = tr["avg200"][idx]

            LC["LC_surfPix_235K"][i, j] = tr["pix_235"][idx]
            LC["LC_surfkm2_235K"][i, j] = tr["km2_235"][idx]
            LC["LC_surfkm2_220K"][i, j] = tr["km2_220"][idx]
            LC["LC_surfkm2_210K"][i, j] = tr["km2_210"][idx]
            LC["LC_surfkm2_200K"][i, j] = tr["km2_200"][idx]

            LC["LC_UTC"][i, j] = int(t.timestamp())

            # instantaneous speed
            if prev_lon is not None:
                R = 6371.0
                dlon = np.deg2rad(lon[idx] - prev_lon)
                dlat = np.deg2rad(lat[idx] - prev_lat)
                a = (np.sin(dlat/2)**2 +
                     np.cos(np.deg2rad(prev_lat)) *
                     np.cos(np.deg2rad(lat[idx])) *
                     np.sin(dlon/2)**2)
                dist_km = 2 * R * np.sqrt(a)
                dt_sec = (t - prev_t).total_seconds()
                LC["LC_velocity"][i, j] = (dist_km*1000)/dt_sec if dt_sec>0 else -999

            prev_lon, prev_lat, prev_t = lon[idx], lat[idx], t

    return LC