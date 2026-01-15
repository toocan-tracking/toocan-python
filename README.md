# TOOCAN - Tracking Organized Deep Convection

TOOCAN (Tracking Organized deep COnvection ANalysis) is a Python implementation of the algorithm used to detect, segment, and track deep convective systems (DCS) using infrared brightness temperature data from geostationary satellites.

This repository contains the open-source, modular version of TOOCAN.
---

## 🌐 Website and References

- **Website**: [https://toocan.ipsl.fr](https://toocan.ipsl.fr)

- **Main Reference**:  
  Fiolleau, T. & Roca, R. (2024). A database of deep convective systems derived from the intercalibrated meteorological geostationary satellite fleet and the TOOCAN algorithm (2012–2020). Earth Syst. Sci. Data, 16, 4021–4050.  
  [https://doi.org/10.5194/essd-16-4021-2024](https://doi.org/10.5194/essd-16-4021-2024)

  Fiolleau, T., R. Roca, S. Cloché, D. Bouniol, P. Raberanto, 2020: Homogenization of geostationary infrared imager channels for cold cloud studies using Megha-Tropiques/ScaRaB. IEEE Trans. Geosci. Remote Sens., vol 58, no. 9, pp. 6609-6622. doi: 10.1109/TGRS.2020.2978171

  Fiolleau, T. and R. Roca, 2013: An Algorithm for the Detection and Tracking of Tropical Mesoscale Convective Systems Using Infrared Images From Geostationary Satellite, IEEE Trans. Geosci. Remote Sens., vol. 51, no. 7, pp. 4302–4315. doi: 10.1109/TGRS.2012.2227762


## 📦 Structure
pyTOOCAN/
├── config/   # Parameter files
├── notebook/ # Jupyter notebooks
├── scripts/  # CLI runner
├── src/
│ └── toocan/
│ ├── main.py
│ ├── io/
│ ├── detection_spreading/
│ ├── preprocessing/
│ └── utils/
├── tests/
└── README.md

## ▶️ How to Run
```bash
PYTHONPATH=src python scripts/run_toocan.py config/fileparam_TOOCAN.dat config/fileparam_GEO.dat