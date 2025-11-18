# MAPClean – Microstructurally Adaptive Pixel-Level Cleaning

**A modular MATLAB pipeline for automated cleaning of EBSD datasets**

MAPClean is a modular MATLAB pipeline, built on the MTEX open-source toolbox, for automated cleaning of Electron Backscatter Diffraction (EBSD) datasets. It applies **Mean Angular Deviation (MAD) filtering**, phase and orientation **wild spike removal (WSR)**, and iterative hole filling to produce high-quality, microstructure-consistent EBSD data. The pipeline includes visualisation tools and checkpoint support for long-running datasets.

---
## Features
 - Modular, stage-wise EBSD cleaning
 - **MAD filtering** for noisy pixel removal
 - **Phase WSR** and **Orientation WSR** to correct spikes
 - Iterative, multi-radius **hole filling** for unindexed pixels
 - Automated generation of **phase maps** and **IPF maps**
 - Checkpoint-based workflow for reproducibility and resuming interrupted runs
 - Fully customizable parameters for different datasets

---
## Requirements

- **MATLAB** (R2016b or newer required) *  
- **MTEX Toolbox** ([https://mtex-toolbox.github.io/](https://mtex-toolbox.github.io/))  
- **Image Processing Toolbox** *  
- **Statistics and Machine Learning Toolbox** *  

<sub>* Proprietary software and toolboxes; a valid license is required to run the code.</sub>

 ---

 ## Installation
 1. Clone this repository:
 ```bash
git clone https://github.com/rahulrox9/RahulS_MAPClean
```

2. Add the repository to your MATLAB path:
```matlab
addpath(genpath('path_to_MAPClean'));
 ```

3. Ensure MTEX and MATLAB toolboxes are installed and added to your MATLAB path.

---

## Usage

1. Place raw EBSD `.ctf` files in the `DataFiles` directory.
2. Open `MAPClean.m` and set **stage control flags**:
```matlab
runMAD      = true;
runPhaseWSR = true;
runOriWSR   = true;
runHoleFill = true;
runSaveFile = true;
```

3. Adjust parameters in the `params` structure.
   
4. Run the pipeline:
```matlab
MAPClean
 ```
5. Check outputs in the `exports` directory.
---

## Workflow Overview

 1. **Initialisation** – load EBSD data, set parameters, assign phase colours.
 2. **MAD filtering** – remove high Mean Angular Deviation pixels.
 3. **Data quality assessment** – determine strict vs. relaxed mode.
 4. **Phase WSR** – remove phase wild spikes using neighbourhood analysis.
 5. **Orientation WSR** – remove orientation wild spikes, including twin handling for Anorthite.
 6. **Hole filling** – iterative, radius-based filling of unindexed pixels.
 7. **Export** – cleaned EBSD file.

 ---
## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `thresholdFrac` | 0.75 | Minimum fraction of dominant cluster required to compute mean orientation during WSR |
| `exportRes` | 300 dpi | Resolution for exported figures (phase maps and IPF maps) |
| `madThreshold` | 0.9 rad | Maximum MAD allowed for a pixel; higher MAD pixels are set to notIndexed |
| `radius_phase` | 2 | Neighbourhood radius for phase Wild Spike Removal (WSR) |
| `radius_ori` | 2 | Neighbourhood radius for orientation Wild Spike Removal (WSR) |
| `misTol_ori` | 5° | Maximum misorientation tolerated when comparing neighboring orientations |
| `minFrac_ori` | 0.25 | Minimum fraction of similar neighbours required for orientation WSR |
| `numMaxPasses` | 50 | Maximum iterations for hole filling procedure |
| `radius_fill` | [6 5 4 3 2 1] | Sequence of neighborhood radii used for multi-pass hole filling |
| `min_neighbours` | 3 | Minimum number of valid neighbours required to fill a hole |
| `min_dom_frac` | 0.50 | Minimum fraction of dominant phase among neighbours to fill a hole |
| `phaseFrac` | Adaptive per radius | Two values per radius: minimum neighbour fraction and minimum dominant fraction used to assign a phase during hole filling |

---
## Outputs
- **Cleaned EBSD files** (`*_clean.ctf`)
- **Phase maps** and **IPF maps** as PNGs
- **Checkpoint MAT files** for each stage: MAD, WSR, hole filling
- **Phase statistics** printed to the console
---


## Contributing
- Fork the repository
- Create a new branch for your changes
- Submit a pull request with detailed description

---

## License

This code is licensed under **GPLv3** (see [LICENSE](LICENSE)).
> Note: This project depends on MTEX (GPLv3) and proprietary MATLAB toolboxes. Users must have a valid MATLAB license to run the code.

