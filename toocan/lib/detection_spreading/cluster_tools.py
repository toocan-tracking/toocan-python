# =============================================================================
# File        : cluster_tools.py
# Author      : Thomas Fiolleau
# Project     : TOOCAN - Tracking Organized Deep Convection
# Description : Cluster management routines used during the detection and
#               spreading phases. Reorders cluster structures after temporal
#               overlap handling, reinitializes internal buffers, and aligns
#               Python-side Blob arrays with labels maintained by C routines.
#
# Provides    :
#   - clean_clusters()
#
# Notes :
#   This is tightly linked to the Blob C-struct definition. It must remain
#   consistent with blob.py and liblabel.so for memory integrity.
#
# License     : CNRS / LEGOS 
# =============================================================================

import ctypes

def clean_clusters(data_param, clusters, labels_present, labelMin, nbMax):
    """
    Reorganizes clusters so that:
        new[label - labelMin] = cluster
    Empty slots are initialized with empty Blob(), equivalent to calloc().
    """
    Blob = type(clusters[0])

    # Allocate new table (all zeros – like calloc)
    new = (Blob * nbMax)()
    
    # Allocation des buffers pour CHAQUE blob
    for i in range(data_param.nbMaxCluster):
        new[i].seed_area_perFrame = (ctypes.c_int * data_param.ZSIZE)()
        new[i].labelVoisin       = (ctypes.c_int * 1000)()

    # Place each cluster at correct index
    for cl in clusters:
        if cl.label in labels_present:
            index = cl.label - labelMin
            if 0 <= index < nbMax:
                new[index] = cl   # copy full struct


    return new
