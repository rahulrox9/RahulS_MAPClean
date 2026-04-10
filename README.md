 # MAPClean – Microstructurally Adaptive Pixel-Level Cleaning
 
 **A modular MATLAB pipeline for automated cleaning of Electron Backscatter Diffraction datasets**
 
MAPClean is a modular MATLAB pipeline, built on the open-source MTEX toolbox, for automated cleaning of Electron Backscatter Diffraction (EBSD) datasets. The workflow operates at the pixel level, combining phase-aware and orientation-aware logic to remove noise, correct misindexed pixels, and reconstruct missing data. It employs an intelligent **Strict vs. Relaxed** protocol selection based on data quality to apply the most appropriate cleaning strategy. The pipeline features **Mean Angular Deviation (MAD) filtering**, **sample-mask generation** for physically grounded restriction of the mapped specimen region, **phase and orientation wild spike removal (WSR)**, adaptive **hole filling** using parallelised Breadth-First Search (BFS, Strict mode) or Multi-Pass Filling (MPF, Relaxed mode), and a **protected hole filling** for conservative recovery of uncertain grain interior regions. Parallel execution via MATLAB's Parallel Computing Toolbox is used in the strict BFS hole-filling stage and in the orientation wild-spike detection stage.
 
 ## Key Features
 
 - **Modular Stage Control:** Enable or disable individual pipeline stages via simple boolean flags.
 - **Automated Protocol Selection:** Automatically selects Strict mode (well-indexed data) or Relaxed mode (sparse data) based on the fraction of notIndexed pixels within the sample mask.
 - **MAD Filtering:** Removes unreliable pixels where the Mean Angular Deviation exceeds a user-defined threshold.
 - **Sample Mask Generation:** Physically constrained sample support is estimated directly from indexed pixels. A `sampleMask` is generated to restrict subsequent reassignment and filling to the mapped specimen region.
 - **Phase WSR:** Corrects misindexed pixels by comparing each pixel's phase to its local neighbourhood. Conservative thresholds are used in Strict mode; more aggressive flipping/removal is used in Relaxed mode.
 - **Orientation WSR:** Detects and corrects orientation wild spikes using local misorientation clustering. Includes twin-boundary protection for Anorthite.
 - **BFS Hole Filling (Strict):** Parallelised cluster-based filling. Disconnected hole clusters are identified, screened for sufficient neighbourhood support, and then filled in parallel. Each worker updates its local patch in sequence so that newly filled pixels can support later fills within the same cluster.
 - **MPF Hole Filling (Relaxed):** Multi-pass iterative filling for sparse maps with large connected hole regions. Vectorised Ni and fracDom pre-filters identify valid candidates each pass. Surviving candidates are then checked serially using local dominant-phase and orientation-coherence logic, with in-pass updates after every accepted fill.
 - **Protected Hole Filling:** A final conservative filling stage applied to pixels previously excluded from standard hole filling (e.g. MAD-removed or low-confidence WSR pixels). Uses strict neighbourhood coverage, dominant-phase enforcement, and ring-based orientation validation to recover only physically reliable pixels.
 - **Orientation Clustering:** Both BFS and MPF use hierarchical single-linkage clustering of neighbour quaternions to assign a mean orientation. MPF additionally uses a ring-based fallback when the full neighbourhood does not yield a sufficiently strong dominant orientation cluster.
 - **Visualisation:** Phase maps and Inverse Pole Figure (IPF) maps are exported automatically at every major pipeline stage.
 - **Checkpointing:** Each stage saves a `.mat` checkpoint. Subsequent runs load from the last saved checkpoint, allowing interrupted runs to resume without reprocessing.

## Parameters

| Parameter | Default | Description |
|:---|:---|:---|
| `exportRes` | 300 | Figure export resolution (dpi) |
| `madThreshold` | 1.0 | MAD filter threshold (radians) |
| `qcfrac` | 0.6 | notIndexed fraction threshold for Strict/Relaxed selection |
| `radius_phase` | 3 | Kernel radius for Phase WSR (pixels) |
| `min_dom_frac` | 0.5 | Minimum dominant phase fraction for Relaxed phase flipping |
| `misTol_ori` | 5° | Misorientation tolerance for orientation clustering |
| `thresholdFrac` | 0.75 | Minimum dominant cluster fraction for orientation assignment |
| `minLead` | 2 | Minimum absolute lead count for Relaxed orientation clustering |
| `scaleLead` | 0.1 | Scaled lead requirement as fraction of neighbourhood size (Relaxed) |
| `minFrac_ori` | 0.25 | Minimum fraction of similar neighbours before a pixel is flagged as a wild spike |
| `radius_ori` | 2 | Kernel radius for Orientation WSR (pixels) |
| `radius_fill_strict` | `[6 5 4 3 2 1]` | Descending radius sequence for BFS hole filling |
| `radius_fill_relaxed` | `[7 6 5 4 3 2 1]` | Descending radius sequence for MPF hole filling |
| `phaseFrac_strict` | radius-specific map | BFS support thresholds stored as `[Ni threshold, fracDom threshold]` |
| `phaseFrac_relaxed` | radius-specific map | MPF support thresholds stored as `[Ni threshold, fracDom threshold]` |
| `coverageFrac`   | 0.60 | Minimum neighbourhood coverage required for protected filling |
| `domFracProt`    | 0.90 | Minimum dominant phase fraction for protected filling |
| `threshFracRing` | 0.75 | Orientation coherence threshold for ring-based validation |

Note: `[Ni threshold, fracDom threshold]`: `[min indexed support / neighbourhood size, min dominant-phase support / min indexed support]` 
### Strict phase-fraction thresholds
- params.phaseFrac_strict(6) = [62/136 45/62];
- params.phaseFrac_strict(5) = [45/100 34/45];
- params.phaseFrac_strict(4) = [30/68 23/30];
- params.phaseFrac_strict(3) = [22/44 17/22];
- params.phaseFrac_strict(2) = [13/20 11/13];
- params.phaseFrac_strict(1) = [6/8 6/6];

### Relaxed phase-fraction thresholds
- params.phaseFrac_relaxed(7) = [72/184 56/72];
- params.phaseFrac_relaxed(6) = [50/136 36/50];
- params.phaseFrac_relaxed(5) = [34/100 26/34];
- params.phaseFrac_relaxed(4) = [24/68 18/24];
- params.phaseFrac_relaxed(3) = [16/44 12/16];
- params.phaseFrac_relaxed(2) = [10/20 8/10];
- params.phaseFrac_relaxed(1) = [4/8 4/4];
 
 ## Requirements
 
 - **MATLAB** (Tested on 2024b)
 - **MTEX Toolbox** (Tested on 6.0.0) – https://mtex-toolbox.github.io/
 - **Parallel Computing Toolbox** – required for strict BFS hole filling and `parfor`-based orientation spike detection
 - **Image Processing Toolbox** – required for `fspecial`, `imdilate`, `bwdist`
 - **Statistics and Machine Learning Toolbox** – required for `linkage`, `cluster`
 
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
 6. Proetected Hole Filling (`doProtectedFilling`)
 
 ## Parallelisation
 
 BFS (Strict):
 - Hole clusters are detected serially.
- Clusters are processed in parallel across workers using a conflict-aware scheduler to prevent overlapping updates.
 - Each worker operates on an isolated patch and updates locally.
 - A convergence loop ensures missed pixels are retried.
 
 MPF (Relaxed):
 - Uses vectorised screening (`conv2`) to identify candidates.
- Filling is applied sequentially with immediate in-pass updates, allowing newly filled pixels to support subsequent decisions within the same pass.
 - Parallelism is limited to screening and orientation detection stages.
 
 Orientation WSR:
 - Wild-spike detection is parallelised using `parfor`.
 
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

1. **MAPClean** – Pixel-level noise removal and data restoration (this script).
2. **GRaMC** – Grain Reconstruction and Multi-stage Cleaning: reconstructs grains, merges twins, removes inclusions.
3. **GRaFT** – Grain-Resolved Fabric and Texture analysis: batch statistical analysis, texture quantification (J/M indices), shape analysis.

## Contributing

1. Fork the repository.
2. Create a new branch for your feature.
3. Submit a pull request with a detailed description of changes.
 
 ## Licence

This code is licensed under **GPL version 3** (see [LICENSE](LICENSE)).

> **Note:** This project depends on MTEX (GPLv3). Users must have valid licences for any proprietary MATLAB toolboxes used.
