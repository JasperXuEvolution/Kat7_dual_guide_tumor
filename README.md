# Kat7 dual guide
project-root/
├── 01_data_collection/
│   ├── auxiliary_code/            # useful but not necessary code
|   ├── main_code/       # necessary code
│   └── data/            # Raw data files, logs, initial data dumps
├── 02_data_cleaning_and_QC/
|   |–– 01-NGS_read_QC.ipynb # main scripts
|   |–– A1-guide_arrray_reconstruction.ipynb # helper scripts
│   ├── data/            # data
│   └── fig/            # figs
├── 03_bootstrapping/
│   ├── 01-Bootstrapping-Kat7_DualGuide-KTHC-Method1.sh # pipeline 
│   |── data/          # data used and result
│   │── main_code/         # Python analysis code
├── config/                 # Global configuration files (e.g., parameters, paths)
├── logs/                   # Consolidated logs from various steps
├── README.md               # Overview and instructions for the project
└── LICENSE                 # Licensing information
