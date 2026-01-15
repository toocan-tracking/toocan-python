import ctypes

#class Blob(ctypes.Structure):
#    _fields_ = [
#        ("label", ctypes.c_int),                      # int label;
#        ("seed_area", ctypes.c_int),                 # int seed_area;
#        ("imin", ctypes.c_ulong),                    # unsigned long imin;
#        ("imax", ctypes.c_ulong),                    # unsigned long imax;
#        ("seed_ThresholdDetection", ctypes.c_int),   # int seed_ThresholdDetection;
#        ("slotBegin", ctypes.c_int),                 # int slotBegin;
#        ("slotEnd", ctypes.c_int),                   # int slotEnd;
#        ("flagFIX", ctypes.c_int),                   # int flagFIX;
#        ("flagDilate", ctypes.c_int),                # int flagDilate;
#        ("seed_duration", ctypes.c_int),             # int seed_duration;
#        ("seed_area_perFrame", ctypes.POINTER(ctypes.c_int)),  # int*
#        ("seed_npixels", ctypes.c_int),              # int seed_npixels;
#        ("Flag_obj", ctypes.c_int),                  # int Flag_obj;
#        ("nb_neighbours", ctypes.c_int),             # int nb_neighbours;
#        ("labelVoisin", ctypes.POINTER(ctypes.c_int)), # int*
#        ("NbMCS_alreadyidentified", ctypes.c_int),   # int NbMCS_alreadyidentified;
#        ("labelMCS_alreadyidentified", ctypes.c_int),# int labelMCS_alreadyidentified;
#        ("flagPrint", ctypes.c_int),                 # int flagPrint;
#        ("flagRelabel", ctypes.c_int)                # int flagRelabel;
#    ]

def make_Blob_class():

    class Blob(ctypes.Structure):
        _fields_ = [
            ("label", ctypes.c_int),
            ("seed_area", ctypes.c_int),
            ("imin", ctypes.c_ulong),
            ("imax", ctypes.c_ulong),
            ("seed_ThresholdDetection", ctypes.c_int),
            ("slotBegin", ctypes.c_int),
            ("slotEnd", ctypes.c_int),
            ("flagFIX", ctypes.c_int),
            ("flagDilate", ctypes.c_int),
            ("seed_duration", ctypes.c_int),

            # POINTEUR, pas tableau inline :
            ("seed_area_perFrame", ctypes.POINTER(ctypes.c_int)),

            ("seed_npixels", ctypes.c_int),
            ("Flag_obj", ctypes.c_int),
            ("nb_neighbours", ctypes.c_int),

            # POINTEUR aussi :
            ("labelVoisin", ctypes.POINTER(ctypes.c_int)),

            ("NbMCS_alreadyidentified", ctypes.c_int),
            ("labelMCS_alreadyidentified", ctypes.c_int),
            ("flagPrint", ctypes.c_int),
            ("flagRelabel", ctypes.c_int),
        ]
    return Blob