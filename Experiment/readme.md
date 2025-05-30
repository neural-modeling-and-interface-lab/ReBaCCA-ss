## Experiment Overview

This repository contains data and analysis scripts for a Delayed-Nonmatch to Sample (DNMS) experiment investigating neural representations across different behavioral conditions and recording sessions. The experiment involved recording spike trains from neurons while animals performed the DNMS task with four trial types:

- **Left Nonmatch (LN)**
- **Left Sample (LS)**
- **Right Nonmatch (RN)**
- **Right Sample (RS)**

The primary goal is to compare neural activity patterns across these conditions and sessions using **ReBaCCA-ss**, followed by visualization with **Multidimensional Scaling (MDS)**. The analysis aims to reveal similarities and differences in neural representations across different sessions or trial types.

## Data Description

The preprocessed data is stored in the `Preprocessed data/` directory. Each file corresponds to a specific recording session and trial type, named as `<session_id>_<trial_type>.mat`, where:

- `<session_id>`: Identifier for the animal_session (e.g., `'1053_6'`, `'1029_1'`).
- `<trial_type>`: One of `'_LEFT_nonmatch'`, `'_LEFT_sample'`, `'_RIGHT_nonmatch'`, `'_RIGHT_sample'`.

Each file contains spike train data for successful trials, organized as a matrix with rows as neurons and columns as time points (in milliseconds).

## How to Run the Analysis

1. **Running the Analysis**:
   - Run `visualizeSpikeRaster.m` for spike raster plots.
   - Run `mainTrialSpecificForMDSAcrossSession.m` to perform ReBaCCA analysis and save results (takes time, the results are already stored).
   - Run `mdsVisualizationAllTogether.m` to generate MDS visualizations.

2. **Parameters**:
   - Key parameters: `alpha = 0.5`, `tol = 1e-2`, `maxIter = 64`, `percentVar = 1-1e-2`. Adjust in scripts as needed.

## Scripts and Their Purposes
- **`visualizeSpikeRaster.m`**  
  Generates spike raster plots for selected trials across datasets and trial types.

- **`mainTrialSpecificForMDSAcrossSession.m`**  
  Main analysis script. Loads spike train data, selects trials, and performs ReBaCCA analysis to compute similarity measures across sessions and trial types. Saves results in `Results/`.

- **`mdsVisualizationAllTogether.m`**  
  Loads ReBaCCA results and visualizes them using MDS. Generates plots showing relationships between conditions and sessions, saved as PDFs in `Fig/`.

- **`gaussian_kernel_smoothing.m`**  
  Function that applies Gaussian kernel smoothing to spike trains, used in ReBaCCA to account for temporal correlations.

- **`plotSpikeRaster.m`**  
  Function to plot spike rasters, visualizing neuron firing times.

- **`export_pdf_figure.m`**  
  Utility function to save MATLAB figures as PDFs with specific formatting.

- **`continuumRegression.m`**  
  Implements continuum regression, a core component of ReBaCCA, based on Xie et al. (2020).
  
- **`continuumRegressionMulti.m`**  
  Extends continuum regression to multi-dimensional responses, used within ReBaCCA.
  
- **`ReBaCCAss.m`**  
  Performs ReBaCCA analysis on pairs of spike train datasets, including permutation tests for significance.


