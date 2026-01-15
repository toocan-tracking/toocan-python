from ctypes import Structure, c_int, c_float, c_char, c_ulong

class DataParam(Structure):
    _fields_ = [
        # --- Période temporelle ---
        ("yearBegin", c_int),
        ("monthBegin", c_int),
        ("dayBegin", c_int),
        ("hourBegin", c_int),
        ("minBegin", c_int),

        ("yearEnd", c_int),
        ("monthEnd", c_int),
        ("dayEnd", c_int),
        ("hourEnd", c_int),
        ("minEnd", c_int),

        # --- Domaine géographique ---
        ("latmin", c_float),
        ("latmax", c_float),
        ("lonmin", c_float),
        ("lonmax", c_float),

        # --- Taille des grilles ---
        ("XSIZE", c_ulong),
        ("YSIZE", c_ulong),
        ("ZSIZE", c_ulong),

        # --- Paramètres TOOCAN ---
        ("deltaDetect", c_float),
        ("deltaSpread", c_float),

        ("timin", c_int),
        ("lifemin", c_int),
        ("labelFirstMCS", c_int),
        ("nbMaxCluster", c_int),
        ("overlap_window_size", c_int),

        # --- Seuils BT ---
        ("minBT", c_int),
        ("maxBT", c_int),
        ("stepBT", c_int),

        # --- Métadonnées ---
        ("version", c_char * 30),
        ("path_out", c_char * 250),
        ("path_fileIN", c_char * 250),
    ]
