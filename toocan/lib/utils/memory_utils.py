# =============================================================================
# File        : memory_utils.py
# Author      : Thomas Fiolleau
# Project     : TOOCAN - Tracking Organized Deep Convection
# Description : Utility functions to evaluate RAM usage, memory availability,
#               memory footprint of arrays, and detect potential overloads
#               before allocating large 3D/4D volumes.
#
# Provides    :
#   - mem_used_mb()          : Current RAM usage of the Python process
#   - mem_available_mb()     : Available system memory in MB
#   - estimate_array_mb()    : Estimate memory footprint of NumPy arrays
#   - warn_if_not_enough_memory() : Safety check before large allocations
#
# License     : CNRS / LEGOS 
# =============================================================================

import psutil
import numpy as np
import os

def mem_used_mb():
    """Current memory usage of this Python process in MB."""
    return psutil.Process(os.getpid()).memory_info().rss / 1024**2

def mem_available_mb():
    """Available system RAM in MB."""
    return psutil.virtual_memory().available / 1024**2

def estimate_array_mb(shape, dtype=np.float32):
    """Memory footprint of a NumPy array."""
    itemsize = np.dtype(dtype).itemsize
    nbytes = np.prod(shape) * itemsize
    return nbytes / 1024**2

def warn_if_not_enough_memory(required_mb, safety_factor=1.5):
    """
    Warn or abort if estimated memory is too high.
    safety_factor = recommended margin (1.5 or 2)
    """
    avail = mem_available_mb()

    print(f"\n🔍 Estimated required memory: {required_mb:.1f} MB")
    print(f"💾 Available memory         : {avail:.1f} MB")

    if required_mb * safety_factor > avail:
        print("\n❌ ERROR: Insufficient memory for the requested computation!")
        print("   Increase RAM or reduce VolumeImage / spatial domain.\n")
        return False

    print("✅ Memory check OK.")
    return True
