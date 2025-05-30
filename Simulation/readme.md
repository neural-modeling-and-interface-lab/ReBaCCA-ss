## Simulation Overview

This simulation evaluates the performance of **ReBaCCA-ss** in comparing neural activity patterns (Figure 6). It generates simulated spike train data with known underlying structures (Gaussian firing rate bumps with 2 scenarios) to test ReBaCCA's ability to recover these structures and assess its robustness through permutation tests.

## How to Run the Simulation
1. **Data Generation**:
   - Run `generateSimulatedData.m` or `generateSimulatedData_scenario2.m` to create simulated data for multiple repeats.

2. **Run Algorithm**:
   - Run `main_loop_for_kernel.m` or `main_scenario2_loop_for_kernel.m` to perform the algorithm to select optimal kernel.
   - Run `main_loop_for_alpha.m` or `main_scenario2_loop_for_alpha.m` to perform the algorithm under different alpha with the optimal kernel.

3. **Analysis**:
   - Run `plotResults.m` or `plotResults_scenario_2.m` to perform the analysis of performance versus kernel widths.
   - Run `plotResultsForAlpha.m` or `plotResultsForAlpha_scenario_2.m` to perform the analysis of performance versus alphas.

3. **Theoretical Correlation from CCA**:
   - Run `checkTheoreticalCorrelation.m` or `checkTheoreticalCorrelationForMatrix.m` to perform the simulation to visualize theoretical correlation from CCA

5. **Parameters**:
   - Key parameters include `kernel_width_pool`, `alphaPool`, `totalRepeat`, etc.
   - Adjust in scripts as needed.

## Scripts and Their Functions

Each script has the version for 2 scenarios.

- **`checkTheoreticalCorrelation.m`**:
  - Validates theoretical CCA correlation calculations for smoothed spike trains.
  
- **`checkTheoreticalCorrelationForMatrix.m`**:
  - Validates theoretical CCA correlation calculations for smoothed spike train matrices.

- **`generateSimulatedData.m`**:
  - Generates spike trains and firing rates for the simulation.
  - Saves data for multiple repeats.

- **`main_loop_for_kernel.m`**:
  - Main analysis script for varying kernel bandwidths.

- **`main_loop_for_alpha.m`**:
  - Similar to `main_loop_for_kernel.m` but varies the `alpha` parameter in ReBaCCA-ss.

- **`plotResults.m`**:
  - Visualizes results, including firing rates, spike rasters, PCA/CCA components, and performance metrics.

- **`export_pdf_figure.m`**:
  - Utility function to save figures as PDFs with customizable styles and sizes.

- **`plotResultsForAlpha.m`**:
  - Visualizes how ReBaCCA results change with varying `alpha` values for a specific simulation scenario.
  
- **`continuumRegression.m`**:
  - Implements continuum regression for a single response variable, a core component of ReBaCCA.

- **`continuumRegressionMulti.m`**:
  - Extends continuum regression to multi-dimensional responses, enhancing ReBaCCA's capabilities.

- **`plotSpikeRaster.m`**:
  - Function to plot spike rasters, visualizing neuron firing times with color-coded neuron groups.

- **`gaussian_kernel_smoothing.m`**:
  - Applies Gaussian kernel smoothing to spike trains, preparing data for analysis.


