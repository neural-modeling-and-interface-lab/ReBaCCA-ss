# ReBaCCA-ss: Relevance-Balanced Continuum Correlation Analysis with Smoothing and Surrogating

**Authors:** Xiang Zhang, Chenlin Xu, Zhouxiao Lu, Haonan Wang, Dong Song  
**Corresponding Authors:** xzhang21@usc.edu, dsong@usc.edu

## Overview

**ReBaCCA-ss** is a method for quantifying similarity between population spiking activities.  
It addresses challenges in neural data analysis by:

- **Balancing alignment and variance explained** using continuum canonical correlation
- **Correcting for noise/artifacts** with surrogate spike trains
- **Automatically selecting optimal kernel bandwidth**

The method is validated on both simulated data and hippocampal recordings from rats.  
For full algorithmic details, see our [paper](https://arxiv.org/abs/2505.13748) (add link here).

---

## Repository Structure

- `Simulation/`  
  MATLAB scripts for generating and analyzing synthetic spike data.

- `Experiment/`  
  MATLAB scripts for analyzing real neural spike recordings.

Each folder contains:
- Main scripts to run simulations/analysis
- Supporting functions for ReBaCCA-ss and data processing

## Usage

### Simulation

To reproduce simulation results:

1. Open MATLAB.
2. Navigate to the `Simulation` folder.
3. Run the main script (update the script name as needed):

    ```matlab
    % Example:
    main_simulation.m
    ```

This will generate synthetic data, compute CCA and ReBaCCA-ss similarity, and display figures.

### Experiment

To analyze real neural data:

1. Place your data files in the `Experiment/data` folder (or update paths in the script).
2. Run the main script:

    ```matlab
    % Example:
    main_experiment.m
    ```

- The script will preprocess the data, apply the ReBaCCA-ss analysis pipeline, and save or plot similarity metrics.
- Parameters such as kernel width range, number of surrogates, and the α value can be set at the top of the main script.

---

## Methodology Summary

1. **Smoothing:** Spike trains are converted to continuous rates using Gaussian convolution.
2. **Dimensionality Reduction:** Latent neural dynamics are extracted via PCA.
3. **Continuum Canonical Correlation:** Alignment balances variance and correlation (tuned by parameter α).
4. **Surrogate Correction:** Similarity from permuted spike trains is subtracted to remove artifacts.
5. **Bandwidth Optimization:** Kernel width maximizing informative similarity is selected.

For full details, see the [paper](#).

---

## Example Output

- Similarity scores (0 to 1) between spike patterns
- Plots showing bandwidth optimization, latent dynamics alignment, and comparison to standard CCA
- Example figures are generated in both simulation and experimental analyses

---

## Citing

If you use this code or the ReBaCCA-ss method, please cite:

@misc{zhang2025rebaccassrelevancebalancedcontinuumcorrelation,
      title={ReBaCCA-ss: Relevance-Balanced Continuum Correlation Analysis with Smoothing and Surrogating for Quantifying Similarity Between Population Spiking Activities}, 
      author={Xiang Zhang and Chenlin Xu and Zhouxiao Lu and Haonan Wang and Dong Song},
      year={2025},
      eprint={2505.13748},
      archivePrefix={arXiv},
      primaryClass={q-bio.NC},
      url={https://arxiv.org/abs/2505.13748}, 
}

---

## License

This project is licensed under License CC BY-NC-SA 4.0.

---

## Contact

For questions, suggestions, or collaborations, please contact:  
- Xiang Zhang: xzhang21@usc.edu  
- Dong Song: dsong@usc.edu

---

