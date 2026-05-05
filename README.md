# 🌩️ TOOCAN

**Tracking Of Organized Convective Algorithm using 3-dimensional segmentatioN**

TOOCAN (Tracking Of Organized Convective Algorithm using 3-dimensional segmentatioN) is a Python implementation of the algorithm used to detect, segment, and track Deep Convective Systems (DCS) using infrared brightness temperature data from geostationary satellites. This branch introduces new tracking capabilities and preprocessing while remaining scientifically consistent with the reference TOOCAN methodology.

This repository contains the open-source, modular implementation of TOOCAN.

---

## 🌐 Website and References

- **Website**: [https://toocan.ipsl.fr](https://toocan.ipsl.fr)

- **Main Reference**:  
  Fiolleau, T. & Roca, R. (2024). A database of deep convective systems derived from the intercalibrated meteorological geostationary satellite fleet and the TOOCAN algorithm (2012–2020). Earth Syst. Sci. Data, 16, 4021–4050.  
  [https://doi.org/10.5194/essd-16-4021-2024](https://doi.org/10.5194/essd-16-4021-2024)

  Fiolleau, T., R. Roca, S. Cloché, D. Bouniol, P. Raberanto, 2020: Homogenization of geostationary infrared imager channels for cold cloud studies using Megha-Tropiques/ScaRaB. IEEE Trans. Geosci. Remote Sens., vol 58, no. 9, pp. 6609-6622. doi: 10.1109/TGRS.2020.2978171

  Fiolleau, T. and R. Roca, 2013: An Algorithm for the Detection and Tracking of Tropical Mesoscale Convective Systems Using Infrared Images From Geostationary Satellite, IEEE Trans. Geosci. Remote Sens., vol. 51, no. 7, pp. 4302–4315. doi: 10.1109/TGRS.2012.2227762


## 🚀 New features

### 🛰️ 1. MCS Tracking

Compatibility with:
Standard geostationary data (30 min resolution)
Rapid-scan service observations on Europe (high temporal resolution)

### 📂 2. File Inventory System (Filelist.nc)

A new preprocessing step generates:

`FileList.nc`

This file provides a complete inventory of input satellite data, including:
- Available files
- File paths
- Missing time steps

Allowing detection of missing data gaps, preventing artificial continuity in tracking results.

### 🔁 3. Restart / Resume Capability

This branch introduces fault tolerance with a possibility to resume after crash or interruption.
It preserves the cluster IDs and temporal continuity.

## 📦 Structure
pyTOOCAN/
├── config/   # Parameter files
├── toocan/
│ └── config/
│ ├── lib/
│ ├── preprocessing/
│ ├── detection_spreading/
│ ├── preprocessing/
│ └── utils/
├── tests/
└── README.md

## ▶️ How to Run
Preprocessing stage to produce a NetCDF file that indicates which GEO files exist, which do not, and includes their respective file paths:
```bash
PYTHONPATH=src python scripts/run_toocan.py config/params_launch_toocan.dat yearEnd monthEnd dayEnd minuteEnd yearBegin monthBegin dayBegin minuteBegin outputPath
```

Launch TOOCAN
```bash
PYTHONPATH=src python config/fileparam_TOOCAN.dat config/fileparam_GEO.dat yearBegin monthBegin dayBegin yearEnd monthEnd dayEnd lonMin lonMax latMin latMax
```

## 🔁 Restart Capability

If execution stops:

Switch parameter `reprise` in `fileparam_TOOCAN.dat` to 1.
➡️ Simply relaunch the same command with the new date.


## License

Copyright (C) 2026 Thomas Fiolleau

This project is licensed under the GNU General Public License v3.0 or later (GPL-3.0-or-later).
See the LICENSE file for details.
