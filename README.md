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

%% ## Workflow Overview
%%
%% 1. **Initialisation** – load EBSD data, set parameters, assign phase colours.
%% 2. **MAD filtering** – remove high Mean Angular Deviation pixels.
%% 3. **Data quality assessment** – determine strict vs. relaxed mode.
%% 4. **Phase WSR** – remove phase wild spikes using neighborhood analysis.
%% 5. **Orientation WSR** – remove orientation wild spikes, including twin handling for Anorthite.
%% 6. **Hole filling** – iterative, radius-based filling of unindexed pixels.
%% 7. **Export** – cleaned EBSD files and visualisations.
%%
%% ---
%%
%% ## Parameters
%%
%% | Parameter | Default | Description |
%% |-----------|---------|-------------|
%% | `madThreshold` | 0.9 rad | Maximum MAD allowed for a pixel |
%% | `radius_phase` | 2 | Neighborhood radius for phase WSR |
%% | `radius_ori` | 2 | Neighborhood radius for orientation WSR |
%% | `misTol_ori` | 5° | Orientation misorientation tolerance |
%% | `numMaxPasses` | 50 | Maximum iterations for hole filling |
%% | `radius_fill` | [6 5 4 3 2 1] | Multi-radius hole filling sequence |
%% | `phaseFrac` | Adaptive | Minimum fraction of neighbors to assign a phase |
%%
%% > Full list of adjustable parameters is in the `params` struct inside the main script.
%%
%% ---
%%
%% ## Outputs
%%
%% - **Cleaned EBSD files** (`*_clean.ctf`)
%% - **Phase maps** and **IPF maps** as PNGs
%% - **Checkpoint MAT files** for each stage: MAD, WSR, hole filling
%% - **Phase statistics** printed to the console
%%
%% ---
%%
%% ## Example
%%
%% Before and after cleaning phase maps and IPF maps can be generated automatically in the `exports` directory. Example visualizations illustrate how MAPClean removes noise and fills unindexed regions while preserving microstructure.
%%
%% ---
%%
%% ## Contributing
%%
%% - Fork the repository
%% - Create a new branch for your changes
%% - Submit a pull request with detailed description
%%
%% ---
%%
%% ## License
%%
%% MIT License – see [LICENSE](LICENSE)
