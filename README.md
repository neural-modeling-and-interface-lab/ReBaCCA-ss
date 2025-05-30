# ReBaCCA-ss: Relevance-Balanced Continuum Correlation Analysis with Smoothing and Surrogating

**Authors:** Xiang Zhang, Chenlin Xu, Zhouxiao Lu, Haonan Wang, Dong Song  
**Corresponding Authors:** Xiang Zhang, Dong Song

## Overview

**ReBaCCA-ss** is a method for quantifying similarity between population spiking activities.  
It addresses challenges in neural data analysis by:

- **Balancing alignment and variance explained** using continuum canonical correlation
- **Correcting for noise/artifacts** with surrogate spike trains
- **Automatically selecting optimal kernel bandwidth**

The method is validated on both simulated data and hippocampal recordings from rats.  
For full algorithmic details, see our [paper](https://arxiv.org/abs/2505.13748).

---

## Repository Structure

- `Simulation/`  
  MATLAB scripts for generating and analyzing synthetic spike data.

- `Experiment/`  
  MATLAB scripts for analyzing real neural spike recordings.

Each folder contains:
- Readme file for code instruction
- Main scripts to run simulations/analysis
- Supporting functions for ReBaCCA-ss and data processing

## Citing

If you use this code or the ReBaCCA-ss method, please cite:

@misc{zhang2025rebaccassrelevancebalancedcontinuumcorrelation,
      title={ReBaCCA-ss: Relevance-Balanced Continuum Correlation Analysis with Smoothing and Surrogating for Quantifying Similarity Between Population Spiking Activities}, 
      author={Xiang Zhang and Chenlin Xu and Zhouxiao Lu and Haonan Wang and Dong Song},
      year={2025},
      eprint={2505.13748},
      archivePrefix={arXiv},
      primaryClass={q-bio.NC},
      url={ https://arxiv.org/abs/2505.13748 }, 
}

---

## License

This project is licensed under CC BY-NC-SA 4.0.

---

## Contact

For questions, suggestions, or collaborations, please contact:  
- Xiang Zhang: xzhang21@usc.edu  
- Dong Song: dsong@usc.edu

---

