# MAPClean – Microstructurally Adaptive Pixel-Level Cleaning
> Work in progress ...

**A modular MATLAB pipeline for automated cleaning of Electron Backscatter Diffraction datasets**

MAPClean is a modular MATLAB pipeline, built on the open-source MTEX toolbox, for automated cleaning of Electron Backscatter Diffraction datasets. It applies **Mean Angular Deviation filtering**, **phase and orientation wild spike removal**, and iterative hole filling using **Breadth-First Search based cluster discovery** to produce high-quality, microstructure-consistent EBSD data. The pipeline includes visualisation tools and checkpoint support for long-running datasets.


## Features
- Modular, stage-wise cleaning of Electron Backscatter Diffraction datasets  
- **Mean Angular Deviation filtering** for removal of noisy pixels  
- **Phase wild spike removal** to correct misindexed pixels  
- **Orientation wild spike removal** to correct orientation spikes, including handling of twins  
- Iterative, multi-radius **filling of unindexed pixels** using **Breadth-First Search (BFS)** cluster discovery  
- Protected pixels are excluded from cluster discovery and filling  
- Only clusters that have pixels filled above a minimum threshold are logged to prevent excessive console output  
- Automated generation of **phase maps** and **inverse pole figure maps**  
- Checkpoint-based workflow for reproducibility and resuming interrupted runs  
- Fully customisable parameters for different datasets  

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
3. Ensure MTEX and MATLAB toolboxes are installed and added to your MATLAB path.

## Usage
1. Place raw EBSD `.ctf` files in the `DataFiles` directory.  
2. Open `MAPClean.m` and set stage control flags:
```matlab
runMeanAngularDeviationFilter      = true;
runPhaseWildSpikeRemoval           = true;
runOrientationWildSpikeRemoval     = true;
runHoleFilling                     = true;
runSaveFile                        = true;
```
3. Adjust parameters in the `params` structure.  
4. Run the pipeline:
```matlab
MAPClean
```
5. Check outputs in the `exports` directory.

## Workflow Overview
1. **Initialisation** – load Electron Backscatter Diffraction data, set parameters, and assign phase colours.  
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
- At the end of each cleaning stage (Mean Angular Deviation filtering, phase wild spike removal, orientation wild spike removal, and hole filling), the pipeline automatically saves a checkpoint as a `.mat` file in the `exports` directory.
- These checkpoint files store the current state of the EBSD data, phase maps, and orientation data.
- If a stage is skipped (for example, by setting its control flag to `false`), the pipeline can automatically load the corresponding checkpoint to resume processing from the last completed stage.
- This ensures that long-running datasets do not need to be reprocessed from the beginning if interrupted.
- Checkpoints are named according to the stage: `*_ebsd_mad.mat`, `*_ebsd_phase.mat`, `*_ebsd_ori.mat`,`*_ebsd_fill.mat`.

## Parameters
| Parameter | Default | Description |
|-----------|---------|-------------|
| `thresholdFrac` | 0.75 | Minimum fraction of dominant cluster required to compute mean orientation during orientation wild spike removal |
| `exportRes` | 300 dots per inch | Resolution for exported figures (phase maps and inverse pole figure maps) |
| `madThreshold` | 0.9 radians | Maximum deviation allowed for a pixel; higher deviation pixels are set to unindexed |
| `radius_phase` | 2 | Neighbourhood radius for phase wild spike removal |
| `radius_ori` | 2 | Neighbourhood radius for orientation wild spike removal |
| `misTol_ori` | 5° | Maximum misorientation tolerated when comparing neighbouring orientations |
| `minFrac_ori` | 0.25 | Minimum fraction of similar neighbours required for orientation wild spike removal |
| `radius_fill` | [6 5 4 3 2 1] | Sequence of neighbourhood radii used for multi-pass BFS filling of unindexed pixels |
| `min_neighbours` | 3 | Minimum number of valid neighbouring pixels required to fill a pixel in a cluster |
| `min_dom_frac` | 0.50 | Minimum fraction of dominant phase among neighbours required to fill a pixel |
| `phaseFrac` | set individually for each radius | Adaptive phase fraction based on the neighbourhood radius; two-element vector `[a b]`, where `a` is the minimum fraction of indexed neighbours required for hole filling, and `b` is the minimum fraction of the dominant phase among neighbours |


## Outputs
- **Cleaned EBSD files** (`*_clean.ctf`)  
- **Phase maps** and **inverse pole figure maps** in PNG format  
- **Checkpoint MAT files** for each stage: Mean Angular Deviation filter, wild spike removal, and hole filling  
- **Cluster-level statistics** printed in the console for clusters meeting `logClusterSizeThreshold`  

## Logging
- Cluster-level statistics are printed in the following format:  
 ```
Cluster <cluster identifier>: filled <number of pixels filled>/<total number of pixels in cluster>
 ```  
- Only clusters where at least 10 connected holes exist are displayed.  
- For large datasets, the output is captured using MATLAB `diary` and exported to text files for later review.

## Contributing
- Fork the repository  
- Create a new branch for your changes  
- Submit a pull request with a detailed description  

## Licence
This code is licensed under **GPL version 3** (see [LICENSE](LICENSE)).  
> Note: This project depends on MTEX (GPLv3) and proprietary MATLAB toolboxes. Users must have a valid MATLAB licence to run the code.
