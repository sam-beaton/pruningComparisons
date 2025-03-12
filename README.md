# pruningComparisons

## Overview
This repository contains scripts for fNIRS data pruning  using different methods and parameters, plus analysis of the effects on data quality and retention.

- **MATLAB (`pruningComparisonsMatlab/`)**: Processes raw data and generates output tables.
- **R (`pruningComparisonsMatlab/`)**: Runs multilevel models (MLMs) on processed data. Also includes a parameter choice script.

## Repository Structure
```
pruningComparisons/
│── pruningComparisonsMatlab/
│   ├── pruneRunAll.m
│   ├── +pruneTools/ # modular helpfer functions
|   ├── collationScripts/ # collates results form pruning into tables which are called using R scripts
│   ├── stats/  # Folder for statistical outcomes
│── pruningComparisonsR/
│   ├── compareCVandSCI.R #compares CV, SCI Only and full QT-NIRS pruning
│   ├── prepBrightMotionMLM.R # uses MLM to invesitgate predictorsd associated with motion occurrence
│   ├── prepBrightSciPspMLM.R # uses MLM to investigate effects on average SCi and PSP measures
│   ├── prepBrightQtParamsMLM  # uses MLM to invesitgate effects of changing SCI- and PSP Thresholds in QT-NIRS
│   ├── parameterChoice.R # script to aid QT-NIRS parameter selection
│── README.md

```

## Requirements
- MATLAB (Tested on version R2022b)
- R (Tested on version 4.4.1)
- Required Matlab toolboxes: qt-nirs, Homer2
- Required R packages: `lme4`, `tidyr`, `foreach`, `doParallel`, `dplyr`, `effectsize`, `MuMIn`, `ggplot2`, `viridis`, `arules`, `data.table`, `effects`, `lmerTest`, `scales`, `Matrix`,  `gridExtra`, `MASS`, `car`,    (install using `install.packages("package_name")`)


## License
[MIT License](LICENSE)