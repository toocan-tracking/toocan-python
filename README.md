🌩️ TOOCAN — Tracking Of Organized Convective Algorithm using 3-dimensional segmentatioN 
 
TOOCAN (Tracking Of Organized Convective Algorithm using 3-dimensional segmentatioN) is a Python implementation of the algorithm used to detect, segment, and track Deep Convective Systems (DCS) using infrared brightness temperature data from geostationary satellites.

This repository contains the open-source, modular implementation of TOOCAN.

🌐 Official Website

🔗 https://toocan.ipsl.fr

📚 Scientific References

If you use TOOCAN in scientific work, please cite:

📄 Fiolleau, T. and R. Roca, 2013: An Algorithm for the Detection and Tracking of Tropical Mesoscale Convective Systems Using Infrared Images From Geostationary Satellite, IEEE Trans. Geosci. Remote Sens., vol. 51, no. 7, pp. 4302–4315.
🔗 https://doi.org/10.1109/TGRS.2012.2227762

📄 Fiolleau et al. (2020)
Fiolleau, T., R. Roca, S. Cloché, D. Bouniol, P. Raberanto, 2020: Homogenization of geostationary infrared imager channels for cold cloud studies using Megha-Tropiques/ScaRaB. IEEE Trans. Geosci. Remote Sens., vol 58, no. 9, pp. 6609-6622. 
🔗 https://doi.org/10.1109/TGRS.2020.2978171

📄 Fiolleau, T. and R. Roca, 2024: A database of deep convective systems derived from the intercalibrated meteorological geostationary satellite fleet and the TOOCAN algorithm (2012–2020), Earth Syst. Sci. Data, 16, 4021–4050.
🔗 https://doi.org/10.5194/essd-16-4021-2024


📦 Repository Structure
toocan-python/
├── config/        # Parameter files
├── notebooks/     # Jupyter notebooks
├── scripts/       # Command-line runners
├── src/
│   └── toocan/
│       ├── main.py
│       ├── io/
│       ├── detection_spreading/
│       ├── preprocessing/
│       └── utils/
├── tests/
├── pyproject.toml
└── README.md

⚙️ Installation
1️⃣ Clone the repository
git clone https://gitlab.in2p3.fr/thomas.fiolleau1/toocan-python.git
cd toocan-python

2️⃣ Install in editable mode
pip install -e .

3️⃣ Compile the C extension (required)
cd src/toocan/detection_spreading
gcc -O3 -fPIC -shared -o label.so label.c


⚠️ Make sure gcc is installed on your system.

▶️ Running TOOCAN
PYTHONPATH=src python scripts/run_toocan.py \
config/fileparam_TOOCAN.dat \
config/fileparam_GEO.dat \
yearBegin monthBegin dayBegin hourBegin minBegin \
yearEnd monthEnd dayEnd hourEnd minEnd \
lonmin lonmax latmin latmax

📂 Creating FileTracking Output
PYTHONPATH=src python scripts/run_write_filetracking.py \
config/fileparam_TOOCAN.dat \
config/fileparam_GEO.dat \
yearBegin monthBegin dayBegin hourBegin minBegin \
yearEnd monthEnd dayEnd hourEnd minEnd \
lonmin lonmax latmin latmax

🌍 Contributing

Contributions are welcome via GitLab Merge Requests.

Workflow

🍴 Fork the repository

🌿 Create a feature branch

📤 Submit a Merge Request to main

All scientific modifications should clearly document:

algorithmic changes

validation strategy

impact on results

⚠️ Scientific Disclaimer

This software is provided for research purposes.

Users are responsible for verifying scientific validity for their specific application and dataset.

📥 Get the Code
git clone https://gitlab.in2p3.fr/thomas.fiolleau1/toocan-python.git

📜 License

Copyright (C) 2026 Thomas Fiolleau

This project is licensed under the GNU General Public License v3.0 or later (GPL-3.0-or-later).
See the LICENSE file for details.

