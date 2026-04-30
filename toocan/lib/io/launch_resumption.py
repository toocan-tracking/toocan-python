import os
import h5py
import numpy as np
from datetime import datetime, timedelta


def launch_resumption(current_start, data_param, df, temporalresolution, global_label_volume, nomenclature="ToocanCloudMask_"):
    """
   Resume a processing run using overlap data.

    Load previous time steps within the overlap window, read corresponding NetCDF
    files, and update the global label volume accordingly.

    Parameters
    ----------
    current_start : datetime
        Start time of the current run.
    data_param : object
        Configuration object containing paths and processing parameters.
    df : pandas.DataFrame
        DataFrame containing at least 'datetime' and 'full_path' columns.
    temporalresolution : int
        Temporal resolution in minutes.
    global_label_volume : ndarray
        Array to be filled with retrieved labels.
    nomenclature : str, optional
        File name prefix, by default "ToocanCloudMask_"

    Returns
    -------
    ndarray
        Updated global label volume.
    """
    dates = df[(df["datetime"] >= current_start - timedelta(minutes=data_param.overlap_window_size * temporalresolution)) & (df["datetime"] < current_start)]
    for i in range(data_param.overlap_window_size):
        date = dates["datetime"].values[i]
        path = dates["full_path"].values[i]
        dt = date.astype('datetime64[s]').astype(object)
        day_str = str(date)[:10].replace('-', '')
        minute_str = str(dt.hour).zfill(2) + str(dt.minute).zfill(2)
        path = os.path.join(data_param.path_out.decode() + 'toocan_' + data_param.version.decode(), 'EUROPE', str(dt.year), str(dt.year) + '_' + str(dt.month).zfill(2) + '_' + str(dt.day).zfill(2), nomenclature + day_str + '_' + minute_str + '.nc')
        with h5py.File(path,'r') as f:
            print(path)
            id_toocan = f['/DCS_number'][:]
            id_toocan[id_toocan == -998] = 0
        global_label_volume[data_param.ZSIZE - data_param.overlap_window_size + i] = id_toocan

    return global_label_volume
