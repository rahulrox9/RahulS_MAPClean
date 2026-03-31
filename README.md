 # MAPClean – Microstructurally Adaptive Pixel-Level Cleaning
 
 **A modular MATLAB pipeline for automated cleaning of Electron Backscatter Diffraction datasets**
 
 MAPClean is a modular MATLAB pipeline, built on the open-source MTEX toolbox, for automated cleaning of Electron Backscatter Diffraction (EBSD) datasets. It employs an intelligent **Strict vs. Relaxed** protocol selection based on data quality to apply the most appropriate cleaning strategy. The pipeline features **Mean Angular Deviation (MAD) filtering**, **sample-mask generation** for physically grounded restriction of the mapped specimen region, **phase and orientation wild spike removal (WSR)**, and adaptive **hole filling** using parallelised Breadth-First Search (BFS, Strict mode) or Multi-Pass Filling (MPF, Relaxed mode). Parallel execution via MATLAB's Parallel Computing Toolbox is used in the strict BFS hole-filling stage and in the orientation wild-spike detection stage.
 
 ## Key Features
 
 - **Modular Stage Control:** Enable or disable individual pipeline stages via simple boolean flags.
 - **Automated Protocol Selection:** Automatically selects Strict mode (well-indexed data) or Relaxed mode (sparse data) based on the fraction of notIndexed pixels within the sample mask.
 - **MAD Filtering:** Removes unreliable pixels where the Mean Angular Deviation exceeds a user-defined threshold.
 - **Sample Mask Generation:** Physically constrained sample support is estimated directly from indexed pixels. A `sampleMask` is generated to restrict subsequent reassignment and filling to the mapped specimen region.
 - **Phase WSR:** Corrects misindexed pixels by comparing each pixel's phase to its local neighbourhood. Conservative thresholds are used in Strict mode; more aggressive flipping/removal is used in Relaxed mode.
 - **Orientation WSR:** Detects and corrects orientation wild spikes using local misorientation clustering. Includes twin-boundary protection for Anorthite.
 - **BFS Hole Filling (Strict):** Parallelised cluster-based filling. Disconnected hole clusters are identified, screened for sufficient neighbourhood support, and then filled in parallel. Each worker updates its local patch in sequence so that newly filled pixels can support later fills within the same cluster.
 - **MPF Hole Filling (Relaxed):** Multi-pass iterative filling for sparse maps with large connected hole regions. Vectorised Ni and fracDom pre-filters identify valid candidates each pass. Surviving candidates are then checked serially using local dominant-phase and orientation-coherence logic, with in-pass updates after every accepted fill.
 - **Orientation Clustering:** Both BFS and MPF use hierarchical single-linkage clustering of neighbour quaternions to assign a mean orientation. MPF additionally uses a ring-based fallback when the full neighbourhood does not yield a sufficiently strong dominant orientation cluster.
 - **Visualisation:** Phase maps and Inverse Pole Figure (IPF) maps are exported automatically at every major pipeline stage.
 - **Checkpointing:** Each stage saves a `.mat` checkpoint. Subsequent runs load from the last saved checkpoint, allowing interrupted runs to resume without reprocessing.
 
 ## Requirements
 
 - **MATLAB** (Tested on 2024b)
 - **MTEX Toolbox** (Tested on 6.0.0) — https://mtex-toolbox.github.io/
 - **Parallel Computing Toolbox** — required for strict BFS hole filling and `parfor`-based orientation spike detection
 - **Image Processing Toolbox** — required for `fspecial`, `imdilate`, `bwdist`
 - **Statistics and Machine Learning Toolbox** — required for `linkage`, `cluster`
 
 > **Note:** Proprietary MATLAB toolboxes require a valid licence to run.

## Installation & Setup

Since MAPClean is part of a three-stage ecosystem, we recommend setting up a single project workspace.

1. Create a new folder (e.g., `My_EBSD_Project/`).
2. Place `MAPClean.m` in that folder.
3. (Recommended) Place **GRaMC** and **GRaFT** in the same folder.
4. Add the folder to your MATLAB path:
```matlab
addpath(genpath('path_to_MAPClean'));
```
5. Ensure MTEX is installed and initialised (`startup_mtex`).
 
 ## Usage
 
 1. Place your raw EBSD `.ctf` files in the `DataFiles/` directory.
 2. Open `MAPClean.m` and set the Stage Control Flags.
 3. Set the sample list.
 4. Adjust global parameters if required.
 5. Run the script.
 6. Outputs are written to diagnostics folders and optionally exported as cleaned `.ctf`.
 
 ## Pipeline Stages
 
 1. MAD Filter (`doMADFilter`)
 2. Sample Mask (`buildSampleMask`)
 3. Phase WSR (`doPhaseWSR_strict` / `doPhaseWSR_relaxed`)
 4. Orientation WSR (`doOrientationWSR`)
 5. Hole Filling (`doHoleFillingBFS` / `doHoleFillingMPF`)
 
 ## Parallelisation
 
 BFS (Strict):
 - Hole clusters are detected serially.
 - Clusters are processed in parallel across workers.
 - Each worker operates on an isolated patch and updates locally.
 - A convergence loop ensures missed pixels are retried.
 
 MPF (Relaxed):
 - Uses vectorised screening (`conv2`) to identify candidates.
 - Filling is applied sequentially with immediate updates.
 - Parallelism is limited to screening and orientation detection stages.
 
 Orientation WSR:
 - Wild spike detection is parallelised using `parfor`.
 - Correction is applied sequentially.
 
 ## Outputs
 
 - Checkpoints (`.mat`)
 - Phase maps and IPF maps (`.png`)
 - Log file
 - Optional cleaned `.ctf`

## Directory Structure

```text
My_EBSD_Project/
├── DataFiles/          # [Input] Raw .ctf files
├── checkpoints/        # [Auto] Stage checkpoints (.mat)
├── exports/
│   ├── MAPClean/       # [Auto] MAPClean plots and logs
│   ├── GrainClean/     # Stage 2 outputs (GRaMC)
│   └── Textures/       # Stage 3 outputs (GRaFT)
├── MAPClean.m          # (Stage 1) This script
├── GRaMC.m             # (Stage 2) Grain Reconstruction
└── GRaFT.m             # (Stage 3) Texture Analysis
```

## The MAPClean Ecosystem

MAPClean is the first step in a modular three-stage pipeline for EBSD analysis.

1. **MAPClean** — Pixel-level noise removal and data restoration (this script).
2. **GRaMC** — Grain Reconstruction and Multi-stage Cleaning: reconstructs grains, merges twins, removes inclusions.
3. **GRaFT** — Grain-Resolved Fabric and Texture analysis: batch statistical analysis, texture quantification (J/M indices), shape analysis.

## Contributing

1. Fork the repository.
2. Create a new branch for your feature.
3. Submit a pull request with a detailed description of changes.
 
 ## Licence

This code is licensed under **GPL version 3** (see [LICENSE](LICENSE)).

> **Note:** This project depends on MTEX (GPLv3). Users must have valid licences for any proprietary MATLAB toolboxes used.
