import xarray as xr
import numpy as np

def write_tracking_nc(outfile, INT, LC, global_times,
                      region="AFRICA",
                      platform="MSG",
                      title="TOOCAN Tracking",
                      version="3.0.0"):

    dcs_ids = LC["DCS_ids"]
    time_seconds = np.array([int(t.timestamp()) for t in global_times], dtype=np.int32)

    ds = xr.Dataset(
        coords={
            "DCS": ("DCS", np.array(dcs_ids, dtype=np.int32)),
            "time": ("time", time_seconds)
        }
    )

    # Fill INT
    for key in INT[dcs_ids[0]]:
        ds["INT_" + key] = ("DCS", np.array([INT[lab][key] for lab in dcs_ids]))

    # Fill LC
    for key in LC:
        if key == "DCS_ids":
            continue
        ds[key] = (("DCS", "time"), LC[key])

    ds.attrs.update({
        "title": title,
        "tracker": "TOOCAN",
        "version": version,
        "region": region,
        "platform": platform
    })

    ds.to_netcdf(outfile)
    print("✔ Saved:", outfile)