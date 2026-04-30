import xarray as xr
import numpy as np

def extract_VZARegcoefs(filename):
    """
    Extract the a, b, c coefficients for VZA correction from a NetCDF file.

    Args:
        filename (str): Path to the NetCDF file.
        nb_VZAcoefs (int): Number of coefficients to extract.

    Returns:
        tuple: (coefVZA_a, coefVZA_b, coefVZA_c)
    """
    dataset = xr.open_dataset(filename)

    # Variable retrieval
    coefVZA_a = dataset.variables["coef_a"]
    coefVZA_b = dataset.variables["coef_b"]
    coefVZA_c = dataset.variables["coef_c"]
    VZAmax = int(dataset.attrs["VZAmax (degree)"])
    BTmax = int(dataset.attrs["BTmax (K)"])

    dataset.close()

    return coefVZA_a, coefVZA_b, coefVZA_c, VZAmax, BTmax

def compute_VZA_correction(coef_a, coef_b, coef_c, nav):

    # Regression coefficients
    toocan_imVZA = nav.mat_ZenithalAngle.values
    mat_coefVZA_ax = np.zeros_like(toocan_imVZA)
    mat_coefVZA_bx = np.zeros_like(toocan_imVZA)
    mat_coefVZA_cx = np.zeros_like(toocan_imVZA)

    # VZA indices
    ivza = toocan_imVZA.astype(int)

    # Mask
    mask = (ivza >= 0) & (ivza < coef_a.size)  

    mat_coefVZA_ax[mask] = coef_a[ivza[mask]]
    mat_coefVZA_bx[mask] = coef_b[ivza[mask]]
    mat_coefVZA_cx[mask] = coef_c[ivza[mask]]

    return mat_coefVZA_ax, mat_coefVZA_bx, mat_coefVZA_cx