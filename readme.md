# Spatial-DC: a robust deep learning-based method for deconvolution of spatial proteomics

<p align="center">
  <img width="80%" src=workflow.jpg>
</p>

Spatial-DC (**S**patial-**D**igital **C**ytometry) is a robust deep learning-based method for deconvolution of spatial proteomics, involving infers cell-type composition in each spot and reconstructs spatially and cell-type resolved proteomic profiles. We achieve this by utilizing transfer learning and self-supervised learning with graph convolutional networks (GCN), which enables the incorporation of target spatial proteomics with reference single-cell or single-cell-type proteomics data.

## Introduction to File Directory
This repository contains various subdirectories and source code files related to downstream analysis. Here's an overview of each subdirectory:

### benchmark_methods
This subdirectory provides the used code for running the benchmark methods on spatial proteomics data. 

### Spatial_DC_V1
The `Spatial_DC_V1` subdirectory contains the source code for running Spatial-DC, which is a specific deep learning-based method for spatial proteomics deconvolution.

### downstream_applications
The `downstream_applications` subdirectory provides the source code for performing downstream analysis and visualization of spatial proteomics data. This includes scripts for data interpretation, visualization techniques, or additional analysis steps that build upon the results obtained from benchmark methods or Spatial-DC.

### tools
The `tools` subdirectory contains code for other common analysis tasks related to spatial proteomics data. This includes utilities for calculating evaluation metrics, transforming data formats (e.g., converting h5seurat to h5ad), or other analysis-related functions that are not specific to a single method or application.


## Tutorial
The details of how to run Spatial-DC are provided [Tutorials]

## Disclaimer
This tool is for research purpose and not approved for clinical use.
This is not an official Tencent product.

## Questions
If you have any suggestions/ideas for Spatial-DC or have issues trying to use it, you can contact us at liy2020@mail.sustech.edu.cn.

[![python >3.8.12](https://img.shields.io/badge/python-3.8.12-brightgreen)](https://www.python.org/) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14386585.svg)](https://doi.org/10.5281/zenodo.14386585)


## Citation
wait for update
