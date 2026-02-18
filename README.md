# 🌩️ TOOCAN  
**Tracking Of Organized Convective Algorithm using 3-dimensional segmentatioN**

TOOCAN (Tracking Of Organized Convective Algorithm using 3-dimensional segmentatioN) is a Python implementation of the algorithm used to detect, segment, and track Deep Convective Systems (DCS) using infrared brightness temperature data from geostationary satellites.

This repository contains the open-source, modular implementation of TOOCAN.

---

## 🌐 Official Website

🔗 https://toocan.ipsl.fr

---

## 📚 Scientific References

If you use TOOCAN in scientific work, please cite:

📄 Fiolleau, T. and R. Roca, 2013: An Algorithm for the Detection and Tracking of Tropical Mesoscale Convective Systems Using Infrared Images From Geostationary Satellite, IEEE Trans. Geosci. Remote Sens., vol. 51, no. 7, pp. 4302–4315.
🔗 https://doi.org/10.1109/TGRS.2012.2227762

📄 Fiolleau et al. (2020)
Fiolleau, T., R. Roca, S. Cloché, D. Bouniol, P. Raberanto, 2020: Homogenization of geostationary infrared imager channels for cold cloud studies using Megha-Tropiques/ScaRaB. IEEE Trans. Geosci. Remote Sens., vol 58, no. 9, pp. 6609-6622. 
🔗 https://doi.org/10.1109/TGRS.2020.2978171

📄 Fiolleau, T. and R. Roca, 2024: A database of deep convective systems derived from the intercalibrated meteorological geostationary satellite fleet and the TOOCAN algorithm (2012–2020), Earth Syst. Sci. Data, 16, 4021–4050.
🔗 https://doi.org/10.5194/essd-16-4021-2024

---

## 📦 Repository Structure

```text
toocan-python/
├── config/        # Parameter files
├── toocan/
│   └── detection_spreading/
│   └── io/
│   └── postprocessing/
│   └── preprocessing/ 
│   └── struct/
│   └── utils/
├── pyproject.toml
└── README.md
```
---

## ⚙️ Installation
1️⃣ Clone the repository
git clone https://gitlab.in2p3.fr/thomas.fiolleau1/toocan-python.git
cd toocan-python

2️⃣ Install in editable mode
pip install -e .

3️⃣ Compile the C extension (required)
cd src/toocan/detection_spreading
gcc -O3 -fPIC -shared -o label.so label.c

⚠️ Make sure gcc is installed on your system.

---

## ⚙️ Configuration Files

TOOCAN requires two configuration files:

- `fileparam_GEOnative.dat`
- `fileparam_TOOCAN.dat`

These files define:

- Satellite / model parameters  
- Data paths  
- Algorithm thresholds  
- Runtime configuration  

---

### 🛰️ GEO Configuration (`fileparam_GEOnative.dat`)

This file defines the geostationary satellite characteristics and input data paths.

#### Main parameters

- `REGION` — Geographic domain (e.g. AFRICA)  
- `GEOplatform` — Satellite platform (e.g. MSG_native)  
- `temporalresolution` — Temporal resolution (minutes)  
- `spatialresolution` — Pixel resolution (degrees)  
- `channel` — Infrared channel used (e.g. IR10.8)  
- `nadir` — Nadir longitude  
- `path_ir` — Path to infrared image files  
- `file_navigation` — Navigation file  

⚠️ These paths must be adapted to your local system.

---

### 🌩️ TOOCAN Algorithm Configuration (`fileparam_TOOCAN.dat`)

This file defines algorithm parameters and output settings.

#### 🔎 General

- `version` — Version number of the TOOCAN algorithm implementation  
  *(Used for traceability and output metadata consistency.)*

---

#### 📁 Output

- `pathout_TOOCAN` — Directory where TOOCAN outputs are written  
  *(Tracking files and intermediate products.)*

---

### 🌡️ Brightness Temperature Thresholds

- `minBT_threshold` — Minimum brightness temperature threshold (K)  
  Defines the coldest threshold used to initiate segmentation.  
  Lower values restrict detection to the coldest convective cores.

- `maxBT_threshold` — Maximum brightness temperature threshold (K)  
  Upper bound for the spreading process.  
  Identifies the warm boundary limit of detected systems.

- `stepBT_threshold` — Increment (K) between successive brightness temperature levels during segmentation.

- `deltaBT_Spread` — Maximum brightness temperature difference allowed between adjacent pixels for the expansion to proceed.

---

### 🌱 Convective Core Detection

- `minAreaSeed` — Minimum area (pixels) required for a seed to be considered a valid convective cluster.  
  Prevents noise or small cold artifacts from being classified as DCS.

- `minLifetime` — Minimum lifetime (number of time steps) required for a system to be retained.  
  Filters short-lived transient systems.

- `firstlabel` — Starting label index used for cluster identification.  
  Primarily used for internal tracking consistency.

---

### 📦 Image Volume Processing

These parameters control temporal tracking and memory handling.

- `VolumeImage` — Number of time steps processed in a single volume block.  
  - Higher values → better temporal continuity, higher memory usage  
  - Lower values → lower memory footprint, possible edge effects  

- `overlap_window_size` — Number of time steps overlapping between successive volume blocks.  
  Ensures temporal continuity at volume boundaries.

- `nbMaxCluster` — Maximum number of clusters handled within a single volume.  
  Acts as a safety limit to prevent memory overflow.

---

⚠️ **Scientific note**

These parameters directly control the detection and tracking behaviour of TOOCAN for Deep Convective Systems.

Modifying these values may alter the scientific output and should be carefully validated.
---

## ▶️ Running TOOCAN
PYTHONPATH=src python scripts/run_toocan.py \
config/fileparam_TOOCAN.dat \
config/fileparam_GEO.dat \
yearBegin monthBegin dayBegin hourBegin minBegin \
yearEnd monthEnd dayEnd hourEnd minEnd \
lonmin lonmax latmin latmax

## 📂 Creating FileTracking Output
PYTHONPATH=src python scripts/run_write_filetracking.py \
config/fileparam_TOOCAN.dat \
config/fileparam_GEO.dat \
yearBegin monthBegin dayBegin hourBegin minBegin \
yearEnd monthEnd dayEnd hourEnd minEnd \
lonmin lonmax latmin latmax

---

## 🌍 Contributing

Contributions are welcome via GitLab Merge Requests.

Workflow

---

## 🍴 Fork the repository

---

## 🌿 Create a feature branch

---

## 📤 Submit a Merge Request to main


All scientific modifications should clearly document:

algorithmic changes

validation strategy

impact on results

---

## ⚠️ Scientific Disclaimer

This software is provided for research purposes.

Users are responsible for verifying scientific validity for their specific application and dataset.

---

## 📥 Get the Code
git clone https://gitlab.in2p3.fr/thomas.fiolleau1/toocan-python.git

---

## 📜 License

Copyright (C) 2026 Thomas Fiolleau

This project is licensed under the GNU General Public License v3.0 or later (GPL-3.0-or-later).
See the LICENSE file for details.

---