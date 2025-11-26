# MAPClean – Microstructurally Adaptive Pixel-Level Cleaning

**A modular MATLAB pipeline for automated cleaning of Electron Backscatter Diffraction datasets**

MAPClean is a modular MATLAB pipeline, built on the open-source MTEX toolbox, for automated cleaning of Electron Backscatter Diffraction (EBSD) datasets. It applies **Mean Angular Deviation (MAD) filtering**, **phase and orientation wild spike removal**, and iterative **Hole Filling** using Breadth-First Search (BFS) based cluster discovery or Multi-Pass Filling (MPF) to produce high-quality, microstructure-consistent EBSD data. The pipeline includes visualisation tools and checkpoint support for long-running datasets.


## Features
- Modular, stage-wise cleaning of Electron Backscatter Diffraction datasets  
- **Mean Angular Deviation filtering:** Removes noisy pixels based on high mean angular deviation.  
- **Phase Wild Spike Removal:** Corrects misindexed pixels by local phase comparison.
- **Orientation Wild Spike Removal:** Corrects orientation spikes, including dedicated logic for twin handling.  
- **Iterative Hole Filling:** Multi-radius filling of unindexed pixels using:
  -- **Breadth-First Search** (BFS) based cluster discovery for strict cleaning.
  -- **Multi-Pass Filling** (MPF) for relaxed cleaning. 
- Automated generation of **phase maps** and **inverse pole figure (IPF) maps**.  
- Checkpoint-based workflow for reproducibility and resuming interrupted runs.  
- Fully customisable parameters for different datasets.  

## Requirements
- **MATLAB** (version 2016b or newer)  
- **MTEX Toolbox** ([https://mtex-toolbox.github.io/](https://mtex-toolbox.github.io/))  
- **Image Processing Toolbox**  
- **Statistics and Machine Learning Toolbox**  
> Note: Proprietary MATLAB toolboxes require a valid licence to run the code.  

## Installation

1. Clone this repository:
```bash
git clone https://github.com/rahulrox9/RahulS_MAPClean
```
2. Add the repository to your MATLAB path:
 ```matlab
 addpath(genpath('path_to_MAPClean'));
 ```
3. Ensure the MTEX and required proprietary MATLAB toolboxes are installed and properly added to your MATLAB path.

## Usage
1. Place raw EBSD `.ctf` files in the `DataFiles` directory.  
2. Open `MAPClean.m` and set stage control flags:
```matlab
runStart    = true;    % Initial plots
runMAD      = true;    % Mean Angular Deviation Filter
runCrop     = true;    % Sample Mask/Cropping
runPhaseWSR = true;    % Phase Wild Spike Removal
runOriWSR   = true;    % Orientation Wild Spike Removal
runHoleFill = true;    % Standard Hole Filling (BFS/MPF)
runProFill  = true;    % Protected Pixel Filling (Protected Holes)
runSaveFile = true;    % Export final EBSD and parameters
```
3. Adjust parameters in the `params` structure.  
4. Run the pipeline:
```matlab
MAPClean
```
5. Check outputs in the `exports` directory.

## Workflow Overview
1. **Initialisation** – load Electron Backscatter Diffraction data and set parameters.  
2. **Mean Angular Deviation filtering** – remove pixels with high deviation from neighbours.  
3. **Data quality assessment** – determine strict or relaxed cleaning mode based on dataset quality.  
4. **Phase wild spike removal** – remove misindexed pixels by comparing to neighbouring pixels.  
5. **Orientation wild spike removal** – remove pixels with unusual orientations, including twin handling.  
6. **Filling unindexed pixels** – Breadth-First Search based, iterative filling:  
- Clusters of connected unindexed pixels are discovered using 8-connectivity  
- Only pixels that are not protected are included in cluster discovery  
- Each cluster is processed iteratively based on the dominant phase fraction and the minimum number of valid neighbours  
- Cluster information is logged only if pixels are successfully filled (default: clusters larger than ten pixels)  
7. **Export** – save cleaned EBSD file, phase maps, and inverse pole figure maps  

## Checkpoints and Resumable Workflow
- At the end of each cleaning stage, the pipeline automatically saves a checkpoint as a `.mat` file in the `exports` directory.
- These checkpoint files store the current state of the EBSD data, phase maps, and orientation data.
- If a stage is skipped (for example, by setting its control flag to `false`), the pipeline can automatically load the corresponding checkpoint to resume processing from the last completed stage.
- This ensures that long-running datasets do not need to be reprocessed from the beginning if interrupted.
- Checkpoints are named according to the stage: `*_ebsd_mad.mat`, `*_ebsd_phase.mat`, `*_ebsd_ori.mat`,`*_ebsd_fill.mat`,`*_ebsd_pro.mat`.

## Parameters
| Parameter | Default | Description |
|-----------|---------|-------------|
| `exportRes` | 300 | (dots per inch) Resolution for exported figures (phase maps and inverse pole figure maps) |
| `madThreshold` | 0.9 | (radians) Maximum deviation allowed for a pixel; higher deviation pixels are set to unindexed |
| `radius_phase` | 2 | Neighbourhood radius for phase wild spike removal |
| `min_dom_frac` | 0.50 | Minimum fraction of dominant phase among neighbours required to fill a pixel |
| `radius_ori` | 2 | Neighbourhood radius for orientation wild spike removal |
| `misTol_ori` | 5° | Maximum misorientation tolerated when comparing neighbouring orientations |
| `thresholdFrac` | 0.75 | Minimum fraction of dominant cluster required to compute mean orientation during orientation wild spike removal |
| `minLead` | 2 | Minimum lead count for WSR (Relaxed/Aggressive) |
| `scaleLead` | 0.1 | Scaling factor for required lead (Relaxed/Aggressive) |
| `minFrac_ori` | 0.25 | Minimum fraction of similar neighbours required for orientation wild spike removal |
| `radius_fill` | [6 5 4 3 2 1] | Sequence of neighbourhood radii used for multi-pass BFS filling of unindexed pixels |
| `phaseFrac` | set individually for each radius | Adaptive phase fraction based on the neighbourhood radius; two-element vector `[a b]`, where `a` is the minimum fraction of indexed neighbours required for hole filling, and `b` is the minimum fraction of the dominant phase among neighbours |
| `thresholdFracRing` | 2/3 | Minimum dominant cluster fraction in a ring |
| `coverageFrac` | 2/3 | Minimum indexed neighbour coverage in the kernel |

## Outputs
- **Cleaned EBSD files** (`*_clean.ctf`)  
- **Visualisations:** Phase maps and Inverse Pole Figure maps in PNG format.
- **Checkpoints:** fMAT files saved in the `checkpoints` directory for each stage. 
- **Logging:** Cluster-level statistics are printed using MATLAB diary and exported to text files in the exports directory for review. 

## Contributing
- Fork the repository  
- Create a new branch for your changes  
- Submit a pull request with a detailed description  

## Licence
This code is licensed under **GPL version 3** (see [LICENSE](LICENSE)).  
> Note: This project depends on MTEX (GPLv3) and proprietary MATLAB toolboxes. Users must have a valid MATLAB licence to run the code.
