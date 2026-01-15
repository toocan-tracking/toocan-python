# ============================================================================================
# File        : spread.py
# Author      : Thomas Fiolleau
# Date        : 2025-07-04
# Description : Implements the label spreading logic for TOOCAN.
#               Reproduces temporal and spatial growth of detected cloud clusters
#               using brightness temperature constraints.
#
# Functions   : 
#   - spread_labels(): Iteratively propagates cluster labels within valid cloudy regions.
#
# Project     : TOOCAN - Tracking Of Organized Convection Algorithm using 3D segmentatioN 
# License     : 
# ============================================================================================

import numpy as np
from toocan.lib.utils.array_utils import compute_index, valid_coords  # optional, if moved out

def spread_labels_V0(volume_bt, labeled_volume, cloud_mask, delta_spread=1.0):
    """
    Reimplementation of the C Spread function in Python.

    Parameters:
        volume_bt (np.ndarray): 3D brightness temperature array (t, y, x)
        labeled_volume (np.ndarray): 3D label array (same shape as volume_bt)
        cloud_mask (np.ndarray): 3D binary mask (1=allowed to grow, 0=forbidden)
        delta_spread (float): Maximum temperature difference (in Kelvin)

    Returns:
        np.ndarray: Updated labeled_volume after growth
    """
    shape = volume_bt.shape
    T, Y, X = shape
    volume_bt = volume_bt.astype(np.float32)
    labeled_volume = labeled_volume.copy()

    # Flatten arrays for easier indexing
    flat_bt = volume_bt.ravel()
    flat_label = labeled_volume.ravel()
    flat_mask = cloud_mask.ravel()
    size = flat_bt.size

    # Prepare index list of initially labeled pixels
    initial_indices = np.flatnonzero(flat_label > 0).tolist()
    growing_indices = initial_indices.copy()

    # Neighborhood offsets (temporal ±1, spatial 8-connectivity)
    offsets = []
    for dt in [-1, 0, 1]:
        for dy in [-1, 0, 1]:
            for dx in [-1, 0, 1]:
                if dt == dx == dy == 0 or abs(dt) + abs(dx) + abs(dy) > 2:
                    continue
                offsets.append((dt, dy, dx))

    def compute_index(t, y, x):
        return t * (Y * X) + y * X + x

    def valid_coords(t, y, x):
        return 0 <= t < T and 0 <= y < Y and 0 <= x < X

    iteration = 0
    while True:
        nPIX_old = len(growing_indices)
        new_indices = []

        for idx in growing_indices:
            t = idx // (Y * X)
            y = (idx % (Y * X)) // X
            x = (idx % (Y * X)) % X
            label = flat_label[idx]
            center_tb = flat_bt[idx]

            for dt, dy, dx in offsets:
                nt, ny, nx = t + dt, y + dy, x + dx
                if not valid_coords(nt, ny, nx):
                    continue

                nidx = compute_index(nt, ny, nx)
                if (
                    flat_mask[nidx] > 0 and
                    flat_label[nidx] == 0 and
                    flat_bt[nidx] + delta_spread >= center_tb
                ):
                    flat_label[nidx] = label
                    new_indices.append(nidx)

        growing_indices = new_indices
        #print(f"Iteration {iteration}: {len(new_indices)} pixels added")
        iteration += 1

        if not new_indices:
            break

    return flat_label.reshape(shape)




#def spread_labels(volume_bt, labeled_volume, cloud_mask, delta_spread=1.0):
#    """
#    Recode du Spread C en Python : propagation conditionnelle des labels dans un volume 3D.
#
#    Parameters:
#        volume_bt (np.ndarray): 3D array (SIZE, YSIZE, ZSIZE) = (lon, lat, time)
#        labeled_volume (np.ndarray): 3D array (SIZE, YSIZE, ZSIZE)
#        cloud_mask (np.ndarray): 3D array (SIZE, YSIZE, ZSIZE), 1 where propagation is allowed
#        delta_spread (float): max allowed temperature diff for propagation (K)
#
#    Returns:
#        np.ndarray: updated labeled_volume
#    """
#
#    SIZE, YSIZE, ZSIZE = volume_bt.shape  # X, Y, T
#    imLabel = labeled_volume.copy()
#    imCloudShield = cloud_mask
#    im3D_IR = (volume_bt * 100).astype(np.int16)  # Mimic short* in C
#
#    def compute_index(x, y, z):
#        return z * YSIZE * SIZE + y * SIZE + x
#
#    flat_size = SIZE * YSIZE * ZSIZE
#    indice_CloudyPix = np.flatnonzero(imLabel.ravel() > 0).tolist()
#    nPIX = len(indice_CloudyPix)
#    nPIXprec = -1
#
#    while nPIX != nPIXprec:
#        nPIXprec = nPIX
#        new_indices = []
#
#        for k in range(nPIXprec):
#            i = indice_CloudyPix[k]
#            if i == -999:
#                continue
#            z = i // (YSIZE * SIZE)
#            y = (i % (YSIZE * SIZE)) // SIZE
#            x = (i % (YSIZE * SIZE)) % SIZE
#            label = imLabel[x, y, z]
#
#            neighbors = []
#            for dz, dy, dx in [
#                (1, 0, 0), (-1, 0, 0),  # temporal
#                (0, 0, -1), (0, 0, 1),
#                (0, 1, -1), (0, 1, 0), (0, 1, 1),
#                (0, -1, -1), (0, -1, 0), (0, -1, 1)
#            ]:
#                zz, yy, xx = z + dz, y + dy, x + dx
#                if 0 <= zz < ZSIZE and 0 <= yy < YSIZE and 0 <= xx < SIZE:
#                    neighbors.append((xx, yy, zz))
#
#            for xx, yy, zz in neighbors:
#                if imLabel[xx, yy, zz] == 0 and imCloudShield[xx, yy, zz] > 0:
#                    tb_center = im3D_IR[x, y, z] / 100.0
#                    tb_neighb = im3D_IR[xx, yy, zz] / 100.0
#                    if tb_neighb + delta_spread >= tb_center:
#                        imLabel[xx, yy, zz] = label
#                        new_indices.append(compute_index(xx, yy, zz))
#
#        indice_CloudyPix.extend(new_indices)
#        nPIX = len(indice_CloudyPix)
#
#    return imLabel



def spread_labels(volume_bt, labeled_volume, cloud_mask, delta_spread=1.0):

    """
    Recode du Spread C en Python : propagation conditionnelle des labels dans un volume 3D.

    Parameters:
        volume_bt (np.ndarray): 3D array (T, Y, X) of brightness temperature
        labeled_volume (np.ndarray): 3D array (T, Y, X) of labels
        cloud_mask (np.ndarray): 3D array (T, Y, X), 1 where propagation is allowed
        delta_spread (float): max allowed temperature diff for propagation (K)

    Returns:
        np.ndarray: updated labeled_volume
    """
    ZSIZE, YSIZE, XSIZE = volume_bt.shape
    imLabel = labeled_volume.copy()
    imCloudShield = cloud_mask
    im3D_IR = (volume_bt * 100).astype(np.int16)  # Mimic short* in C
    flat_size = ZSIZE * YSIZE * XSIZE

    def compute_index(z, y, x):
        return z * YSIZE * XSIZE + y * XSIZE + x

    # Initial list of labeled pixels
    indice_CloudyPix = np.flatnonzero(imLabel.ravel() > 0).tolist()
    nPIX = len(indice_CloudyPix)
    nPIXprec = -1

    while nPIX != nPIXprec:
        nPIXprec = nPIX
        n = nPIX
        new_indices = []

        for k in range(nPIXprec):
            i = indice_CloudyPix[k]
            if i == -999:
                continue
            z = i // (YSIZE * XSIZE)
            y = (i % (YSIZE * XSIZE)) // XSIZE
            x = (i % (YSIZE * XSIZE)) % XSIZE
            label = imLabel[z, y, x]

            neighbors = []
            for dz, dy, dx in [
                (1, 0, 0), (-1, 0, 0),  # temporal
                (0, 0, -1), (0, 0, 1),
                (0, 1, -1), (0, 1, 0), (0, 1, 1),
                (0, -1, -1), (0, -1, 0), (0, -1, 1)
            ]:
                zz, yy, xx = z + dz, y + dy, x + dx
                if 0 <= zz < ZSIZE and 0 <= yy < YSIZE and 0 <= xx < XSIZE:
                    neighbors.append((zz, yy, xx))

            for zz, yy, xx in neighbors:
                if imLabel[zz, yy, xx] == 0 and imCloudShield[zz, yy, xx] > 0:
                    tb_center = im3D_IR[z, y, x] / 100.0
                    tb_neighb = im3D_IR[zz, yy, xx] / 100.0
                    if tb_neighb + delta_spread >= tb_center:
                        imLabel[zz, yy, xx] = label
                        new_indices.append(compute_index(zz, yy, xx))

        indice_CloudyPix.extend(new_indices)
        nPIX = len(indice_CloudyPix)

    return imLabel



import numpy as np
from numba import njit
from collections import deque

@njit
def spread_labels_fast(volume_bt, labeled_volume, cloud_mask, delta_spread):
    Z, Y, X = volume_bt.shape
    im3D_IR = (volume_bt * 100).astype(np.int16)

    # Copy label volume
    labels = labeled_volume.copy()

    # Build initial queue of labeled pixels
    qz = []
    qy = []
    qx = []
    for z in range(Z):
        for y in range(Y):
            for x in range(X):
                if labels[z, y, x] > 0:
                    qz.append(z)
                    qy.append(y)
                    qx.append(x)

    # Convert lists to deques manually (Numba-friendly)
    head = 0

    # 10-neighbor connectivity, including temporal neighbors
    neigh = np.array([
        [ 1, 0, 0], [-1, 0, 0],
        [ 0, 0,-1], [ 0, 0, 1],
        [ 0, 1,-1], [ 0, 1, 0], [ 0, 1, 1],
        [ 0,-1,-1], [ 0,-1, 0], [ 0,-1, 1]
    ], dtype=np.int32)

    while head < len(qz):
        z = qz[head]
        y = qy[head]
        x = qx[head]
        head += 1

        lbl = labels[z, y, x]
        tb_center = im3D_IR[z, y, x] / 100.0

        for k in range(10):
            zz = z + neigh[k, 0]
            yy = y + neigh[k, 1]
            xx = x + neigh[k, 2]

            # Bounds check
            if zz < 0 or zz >= Z or yy < 0 or yy >= Y or xx < 0 or xx >= X:
                continue

            # Skip already-labeled or not cloud-shielded
            if labels[zz, yy, xx] != 0 or cloud_mask[zz, yy, xx] == 0:
                continue

            # Temperature check
            tb_neigh = im3D_IR[zz, yy, xx] / 100.0
            if tb_neigh + delta_spread >= tb_center:
                labels[zz, yy, xx] = lbl
                qz.append(zz)
                qy.append(yy)
                qx.append(xx)

    return labels