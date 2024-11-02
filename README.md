# Multi-task benchmarking of spatially resolved gene expression simulation models

This repository contains all the scripts required to reproduce the multi-task benchmarking of spatially resolved gene expression simulation models. The repository is organized into four main folders, each corresponding to a distinct stage in the benchmarking pipeline:

<img width="762" alt="overview" src="https://github.com/user-attachments/assets/37094612-d161-429a-ac39-e98b3afae32f">

- 01_run_processing: Scripts for data processing and preparation.
- 02_run_simulation: Scripts for running the simulations by simAdaptor or Non-simAdaptor on the processed data.

<img width="762" alt="Screenshot 2024-11-02 at 17 02 56" src="https://github.com/user-attachments/assets/996da9a5-d0c5-431c-9641-b14e66ca299e">


- 03_run_evaluation: Scripts for evaluating simulation results across various metrics.
- 04_run_msFig: Scripts to generate figures and summaries for the manuscript.

# Processed Datasets
The processed datasets required for reproduction are available on Figshare and can be accessed via this DOI link: https://doi.org/10.6084/m9.figshare.26054188.v2. Please download and store them in the appropriate directories as required by the scripts.

# Reproduction Steps
To reproduce the results:

1. Start with the 01_run_processing scripts to prepare and process the data required for simulations.
2. Run the 02_run_simulation scripts to simulate spatially resolved gene expression data.
3. Run the 03_run_evaluation scripts to evaluate the simulated data using multiple benchmarking metrics.
4. Finally, execute the 04_run_msFig scripts to generate figures and summaries for the manuscript.

# Reference

Liang, X., Cao, Y., Yang, J. Y. H. (2024). Multi-task benchmarking of spatially resolved gene expression simulation models. bioRxiv.
[Preprint]. https://doi.org/10.1101/2024.05.29.596418

Liang, Xiaoqi (2024). SpatialSimBench dataset. figshare. Dataset. https://doi.org/10.6084/m9.figshare.26054188.v2
