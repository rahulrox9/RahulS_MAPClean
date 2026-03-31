# MAPClean – Microstructurally Adaptive Pixel-Level Cleaning

**A modular MATLAB pipeline for automated cleaning of Electron Backscatter Diffraction datasets**

MAPClean is a modular MATLAB pipeline, built on the open-source MTEX toolbox, for automated cleaning of Electron Backscatter Diffraction (EBSD) datasets. It employs an intelligent **Strict vs. Relaxed** protocol selection based on data quality to apply the most appropriate cleaning strategy. The pipeline features **Mean Angular Deviation (MAD) filtering**, **Band Contrast (BC) gating** for physically grounded sample masking, **phase and orientation wild spike removal (WSR)**, and adaptive **hole filling** using parallelised Breadth-First Search (BFS, Strict mode) or Multi-Pass Filling (MPF, Relaxed mode). Parallel execution via MATLAB's Parallel Computing Toolbox is used throughout the hole filling stages to reduce runtime on large maps.

## Key Features

- **Modular Stage Control:** Enable or disable individual pipeline stages via simple boolean flags.
- **Automated Protocol Selection:** Automatically selects Strict mode (well-indexed data) or Relaxed mode (sparse data) based on the fraction of notIndexed pixels within the BC-gated sample mask.
- **MAD Filtering:** Removes unreliable pixels where the Mean Angular Deviation exceeds a user-defined threshold.
- **BC Gating:** Physically motivated sample masking based on Band Contrast. A crossover detection algorithm identifies the BC range where indexed and notIndexed pixel distributions separate. A minimum range check triggers automatic fallback to 1st/99th percentile of indexed BC if the crossover range is too narrow. This produces a `sampleMask` (valid fill candidates) and a `neverFillMask` (permanently excluded pixels), replacing the legacy row/column sweep.
- **Phase WSR:** Corrects misindexed pixels by comparing each pixel's phase to its local neighbourhood. Conservative thresholds in Strict mode; more aggressive flipping in Relaxed mode.
- **Orientation WSR:** Detects and corrects orientation wild spikes using local misorientation clustering. Includes twin boundary protection for Anorthite (Albite, Pericline, Carlsbad, Manebach, and Baveno twin laws).
- **BFS Hole Filling (Strict):** Parallelised cluster-based filling. Disconnected hole clusters are discovered serially via 8-connectivity BFS, then filled in parallel (`parfor` over clusters). Intra-cluster fills propagate sequentially in BFS queue order. A convergence loop repeats discovery and fill until no new pixels are filled, catching pixels whose neighbourhood was improved by adjacent clusters.
- **MPF Hole Filling (Relaxed):** Multi-pass iterative filling for sparse maps with large connected hole regions. Vectorised Ni and fracDom pre-filters identify valid candidates each pass. Remaining candidates are filled in parallel (`parfor` over pixels within each pass). Passes repeat until convergence or `maxPasses_mpf` is reached.
- **Orientation Clustering:** Both BFS and MPF use hierarchical single-linkage clustering of neighbour quaternions to assign a mean orientation. BFS additionally uses a ring-based fallback: when the full kernel produces a weak dominant cluster, progressively closer rings are checked for stronger spatial coherence.
- **Visualisation:** Phase maps and Inverse Pole Figure (IPF) maps are exported automatically at every pipeline stage.
- **Checkpointing:** Each stage saves a `.mat` checkpoint. Subsequent runs load from the last saved checkpoint, allowing interrupted runs to resume without reprocessing.

## Requirements

- **MATLAB** (Tested on 2024b)
- **MTEX Toolbox** (Tested on 6.0.0) — [Download MTEX](https://mtex-toolbox.github.io/)
- **Parallel Computing Toolbox** — required for `parfor` in BFS and MPF hole filling
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
2. Open `MAPClean.m` and set the **Stage Control Flags**:
```matlab
runStart    = false;   % Initial raw plots
runMAD      = true;    % MAD filter
runBC       = true;    % BC gating (sample mask)
runPhaseWSR = true;    % Phase Wild Spike Removal
runOriWSR   = true;    % Orientation Wild Spike Removal
runFill     = true;    % Hole Filling (BFS or MPF)
```
3. Adjust **Global Parameters** if needed (see table below).
4. Run the pipeline:
```matlab
MAPClean
```
5. Outputs (plots, logs, cleaned `.mat` checkpoints) will appear in the `exports/MAPClean/` directory.

> **Tip:** To rerun only hole filling from a previous WSR checkpoint, set `runMAD`, `runBC`, `runPhaseWSR`, and `runOriWSR` to `false`. The pipeline will load from saved checkpoints automatically.

## Pipeline Stages

| # | Stage | Function | Mode |
|---|-------|----------|------|
| 1 | MAD Filter | `doMADFilter` | Both |
| 2 | BC Gating | `doBCGating` | Both |
| 3 | Data Quality Assessment | (inline) | Both |
| 4 | Phase WSR | `doPhaseWSR_strict` / `doPhaseWSR_Relaxed` | Strict / Relaxed |
| 5 | Orientation WSR | `doOrientationWSR` | Both |
| 6 | Hole Filling | `doHoleFillingBFS` / `doHoleFillingMPF` | Strict / Relaxed |

## Parameters

| Parameter | Default | Description |
|:---|:---|:---|
| `exportRes` | 300 | Figure export resolution (dpi) |
| `madThreshold` | 0.9 | MAD filter threshold (radians) |
| `qcfrac` | 0.6 | notIndexed fraction threshold for Strict/Relaxed selection |
| `bcMinRange` | 30 | Minimum acceptable BC crossover range; triggers percentile fallback if narrower |
| `radius_phase` | 3 | Kernel radius for Phase WSR (pixels) |
| `min_dom_frac` | 0.5 | Minimum dominant phase fraction for Relaxed phase flipping |
| `misTol_ori` | 5° | Misorientation tolerance for orientation clustering |
| `thresholdFrac` | 0.75 | Minimum dominant cluster fraction for WSR orientation assignment |
| `minLead` | 2 | Minimum absolute lead count for Relaxed orientation clustering |
| `scaleLead` | 0.1 | Scaled lead requirement as fraction of neighbourhood size (Relaxed) |
| `minFrac_ori` | 0.25 | Minimum fraction of similar neighbours before a pixel is flagged as a wild spike |
| `radius_ori` | 2 | Kernel radius for Orientation WSR (pixels) |
| `radius_fill_strict` | `[6 5 4 3 2 1]` | Descending radius sequence for BFS hole filling |
| `radius_fill_relaxed` | `[7 6 5 4 3 2 1]` | Descending radius sequence for MPF hole filling |
| `phaseFrac_strict` | `[0.4, 0.75]` | BFS: `[Ni threshold, fracDom threshold]` |
| `phaseFrac_relaxed` | `[0.35, 0.75]` | MPF: `[Ni threshold, fracDom threshold]` |
| `misTol_fill` | 5° | Misorientation tolerance for hole filling orientation clustering |
| `thresholdFrac_fill` | 0.75 | Minimum dominant cluster fraction for hole filling orientation assignment |
| `maxPasses_mpf` | 500 | Maximum passes per radius in MPF |

## Parallelisation

MAPClean uses MATLAB's Parallel Computing Toolbox. A local parallel pool of 10 workers is started automatically at runtime. The following assumptions underpin parallel execution in the hole filling stages:

**BFS (Strict):**
- *Inter-cluster independence:* Disconnected clusters are filled in parallel. Fills in one cluster do not update the neighbourhood seen by another within the same convergence iteration. This is valid because disconnected clusters share no pixels; any information lost at cluster boundaries is recovered in subsequent convergence iterations.
- *Intra-cluster sequential propagation:* Within each cluster, pixels are filled sequentially in BFS queue order. Each worker maintains a local copy of the padded maps, updated live as fills are made, so interior pixels benefit from boundary fills within the same cluster.
- *Convergence loop:* After each parallel fill round, a new cluster discovery pass runs. Pixels that failed due to insufficient neighbourhood context are retried if adjacent clusters have since filled. The loop terminates when no new pixels are filled.
- *Snapshot consistency:* Each worker receives a read-only snapshot of the global maps taken before the parallel fill. Workers write only to their own cluster's pixels; no two clusters share pixels so there are no write conflicts.

**MPF (Relaxed):**
- Vectorised `conv2`-based pre-filters compute Ni (indexed neighbour fraction) and fracDom (dominant phase fraction) for all hole pixels simultaneously at the start of each pass, identifying valid fill candidates without per-pixel looping.
- Valid candidates within each pass are filled in parallel (`parfor` over pixels). Fills within a pass do not update the shared maps mid-pass; the sequential update assumption is accepted as a negligible approximation given that hole pixels do not contribute valid orientation information to each other's neighbourhood.

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
