# MAPClean – Microstructurally Adaptive Pixel-Level Cleaning

**A modular MATLAB pipeline for automated cleaning of Electron Backscatter Diffraction datasets**

MAPClean is a modular MATLAB pipeline, built on the open-source MTEX toolbox, for automated cleaning of Electron Backscatter Diffraction (EBSD) datasets. It employs an intelligent **"Strict" vs. "Relaxed"** protocol selection based on data quality to apply the most appropriate cleaning strategy. The pipeline features **Mean Angular Deviation (MAD) filtering**, **phase and orientation wild spike removal**, and adaptive **Standard (unprotected) Hole Filling** (using Breadth-First Search or Multi-Pass Filling) and **Protected Pixel Filling** to recover valid grains lost during filtering. This combination produces high-quality, microstructure-consistent EBSD data.

## Key Features
* **Modular Stage Control:** Enable or disable specific cleaning steps (MAD, WSR, Hole Filling, etc.) via simple flags.
* **Automated Protocol Selection:** Automatically selects **Strict** mode (for well-indexed data >60%) or **Relaxed** mode (for sparse data <60%) to balance preservation with restoration.
* **MAD Filtering:** Removes noisy pixels based on high Mean Angular Deviation thresholds.
* **Phase Wild Spike Removal (WSR):** Corrects misindexed pixels using local phase comparison (Conservative or Aggressive based on protocol).
* **Orientation WSR:** Corrects orientation spikes with specific twin-boundary handling (currently optimised for Anorthite).
- **Adaptive Hole Filling:**
  - **BFS (Strict):** Breadth-First Search cluster discovery for high-quality datasets.
  - **MPF (Relaxed):** Multi-Pass Filling for recovering data in sparse maps.
* **Protected Pixel Filling:** A dedicated stage to recover "Protected" pixels (real grains lost to earlier filters) using ring-based orientation checks.
* **Visualisation:** Automated generation of **Phase maps** and **Inverse Pole Figure (IPF) maps** at every stage.
* **Checkpointing:** Automatic state saving allows runs to be resumed from the last successful stage if interrupted.

## Requirements
* **MATLAB** (Tested on 2024b)
* **MTEX Toolbox** (Tested on 6.0.0) – [Download MTEX](https://mtex-toolbox.github.io/)
* **Image Processing Toolbox**
* **Statistics and Machine Learning Toolbox**

> **Note:** Proprietary MATLAB toolboxes require a valid licence to run.

## Installation & Setup
Since MAPClean is part of a 3-stage ecosystem, we recommend setting up a single **Project Workspace**.

1. Create a new folder (e.g., `My_EBSD_Project/`).
2. Clone/Download this repository and place `MAPClean.m` in that folder.
3. (Recommended) Download **GRaMC** and **GRaFT** and place them in the same folder.
4. Add the repository to your MATLAB path:
```matlab
addpath(genpath('path_to_MAPClean'));
```
5. Ensure MTEX is installed and initialised (`startup_mtex`).

## Usage
1. Place your raw EBSD `.ctf` files in the `DataFiles` directory.
2. Open `MAPClean.m`. The script is configured to process all `*.ctf` files in the folder.
3. Set the **Stage Control Flags** to determine which steps to run:
```matlab
runStart    = true;    % Initial plots
runMAD      = true;    % MAD Filter
runCrop     = true;    % Sample Mask/Cropping
runPhaseWSR = true;    % Phase Wild Spike Removal
runOriWSR   = true;    % Orientation Wild Spike Removal
runHoleFill = true;    % Hole Filling (BFS/Strict or MPF/Relaxed)
runProFill  = true;    % Protected Pixel Filling
runSaveFile = true;    % Export final data
```
4. Adjust the **Global Parameters** if necessary (see defaults below).
5. Run the pipeline:
```matlab
MAPClean
```
6. Outputs (plots, logs, and cleaned files) will appear in the `exports` directory.

## Workflow Details
1. **Initialisation:** Loads `.ctf` data, assigns phase colours, and initialises parameters.
2. **MAD Filtering:** Pixels with a Mean Angular Deviation > `madThreshold` are set to `notIndexed`.
3. **Cropping:** Generates a sample mask to exclude the mounting background from calculations.
4. **Data Quality Assessment:**
* **Strict Mode:** Activated if indexed fraction > 60%. Uses conservative WSR and BFS Hole Filling.
* **Relaxed Mode:** Activated if indexed fraction < 60%. Uses aggressive WSR and Multi-Pass Filling (MPF).
5. **Wild Spike Removal (Phase & Orientation):**
* Removes single-pixel spikes based on neighbour consensus.
* Includes specific logic to avoid removing valid twin boundaries.
6. **Hole Filling (Unprotected):** Iterative filling using a descending radius sequence (e.g., 7 down to 1).
7. **Protected Pixel Filling:** Specifically targets pixels removed by the MAD filter that are likely part of real grains, using a "Ring-based" consistency check to restore them.
8. **Export:** Saves the final `_clean.ctf` file and parameter logs.

## Parameters
The parameters are stored in a global `params` structure. Defaults used in the code are:

| Parameter | Default | Description |
| :--- | :----- | :--- |
| `madThreshold` | 0.9 rad | Threshold for MAD filtering. |
| `radius_phase` | 3 | Kernel radius for Phase WSR. |
| `min_dom_frac` | 0.5 | Min dominant phase fraction for Relaxed phase flipping. |
| `radius_ori` | 2 | Kernel radius for Orientation WSR. |
| `misTol_ori` | 5° | Misorientation tolerance for neighbour comparison. |
| `radius_fill` | `[6 5 4 3 2 1]` | Descending radii sequence for iterative hole filling. |
| `phaseFrac` | `Map` | A container Map linking radius to `[Ni, Frac]`. Default is `[0.4, 0.75]` (40% indexed neighbours, 75% dominant phase required). |
| `thresholdFracRing` | 2/3 | Min dominant cluster fraction required in a ring (Protected Fill). |
| `coverageFrac` | 0.5 | Min indexed neighbour coverage required (Protected Fill). |

## Directory Structure
To ensure data flows automatically between stages, your folder should look like this:
```text
My_EBSD_Project/
├── DataFiles/          # [Input] Place raw .ctf files here
├── checkpoints/        # [Shared] Auto-generated intermediate saves
├── exports/            # [Shared] Auto-generated outputs
│   ├── MAPClean/       # Stage 1 Outputs
│   ├── GrainClean/     # Stage 2 Outputs (GRaMC)
│   └── Textures/       # Stage 3 Outputs (GRaFT)
├── MAPClean.m          # (Stage 1) This Script
├── GRaMC.m             # (Stage 2) Grain Reconstruction
└── GRaFT.m             # (Stage 3) Texture Analysis
```

## The MAPClean Ecosystem
MAPClean is the first step in a modular three-stage pipeline for EBSD analysis.
1. **MAPClean:** (This Repository) Pixel-level noise removal and data restoration.
2. **GRaMC:** (Grain Reconstruction and Multi-stage Cleaning) Takes MAPClean output to reconstruct grains, merge twins, and remove inclusions.
3. **GRaFT:** (Grain Reporting and Fabric Texture) Takes GRaMC output to perform batch statistical analysis, texture quantification (J/M indices), and shape analysis.

## Contributing
1. Fork the repository.
2. Create a new branch for your feature.
3. Submit a pull request with a detailed description of changes.

## Licence
This code is licensed under **GPL version 3** (see [LICENSE](LICENSE)).
> **Note:** This project depends on MTEX (GPLv3). Users must have valid licences for any proprietary MATLAB toolboxes used.
