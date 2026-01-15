import numpy as np

def compute_cold_surfaces(bt, mask, lab, surface_area_2d):
    m = (mask == lab)
    out = {}

    for thr in [235, 220, 210, 200]:
        cold = m & (bt <= thr)
        out[f"pix_{thr}K"] = int(cold.sum())
        out[f"km2_{thr}K"] = float(surface_area_2d[cold].sum())

    return out


def compute_temperature_stats(bt, mask, lab):
    m = (mask == lab)
    if not np.any(m):
        return dict(tbmin=np.nan, tb90=np.nan, avg235=np.nan, avg208=np.nan, avg200=np.nan)

    vals = bt[m]

    tb_min = float(np.nanmin(vals))
    tb_90 = float(np.nanpercentile(vals, 90))

    avg235 = float(vals[vals <= 235].mean()) if np.any(vals <= 235) else np.nan
    avg208 = float(vals[vals <= 208].mean()) if np.any(vals <= 208) else np.nan
    avg200 = float(vals[vals <= 200].mean()) if np.any(vals <= 200) else np.nan

    return dict(tbmin=tb_min, tb90=tb_90, avg235=avg235, avg208=avg208, avg200=avg200)


def compute_equivalent_ellipse(bt, mask, lab, thr, lat_2d, lon_2d):
    m = (mask == lab)
    cold = m & (bt <= thr)

    yy, xx = np.nonzero(cold)
    if yy.size < 5:
        return (np.nan, np.nan, np.nan, np.nan)

    y0, x0 = yy.mean(), xx.mean()

    lat0 = float(lat_2d[int(y0), int(x0)])
    lon0 = float(lon_2d[int(y0), int(x0)])

    # Convert to km
    km_lat = 111.0
    km_lon = 111.0 * np.cos(np.deg2rad(lat0))

    x_km = (lon_2d[yy, xx] - lon0) * km_lon
    y_km = (lat_2d[yy, xx] - lat0) * km_lat

    coords = np.vstack([x_km, y_km])
    cov = np.cov(coords)
    eigvals, eigvecs = np.linalg.eig(cov)

    a = 2 * np.sqrt(eigvals.max())
    b = 2 * np.sqrt(eigvals.min())

    vec = eigvecs[:, np.argmax(eigvals)]
    angle = np.degrees(np.arctan2(vec[1], vec[0]))

    return float(a), float(b), float(b/a), float(angle)