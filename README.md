# Notochord Simulation and Figure Generation

## Overview

This repository contains simulation code, data processing scripts, and plotting routines used to generate all figures in the manuscript. The workflow is organized into a simple and reproducible pipeline, separating data generation, data storage, and figure reproduction.

---

## Workflow

### 1. Generate Simulation Data

Run scripts in:
Data generation/

Main entry file:
Main_data_generation.m

This step produces .mat files containing simulation outputs (Geometry, Pressure, Flux, etc.).

### 2.	Prepare data for plotting
Move or copy all generated .mat files into:
Data_for_plot/
This folder serves as the centralized data source for all figure scripts.

### 3.	Generate figures (main text and supplementary)
Run the following scripts in the root directory: 
o	PART1_BaseResults.m 
o	PART2_Propagation_and_TimeScales.m 
o	PART3_AsymmetricResponse.m 
o	PART4_VolumeChange_DifScenarios.m 
o	PART5_DecayProfiles.m 

### 4.	Growth simulation (Figure 2)
Located in:
Fig2_growth/
Main file:
Figure2_growth_main.m

## Folder Structure
### 1. Data generation
Contains core simulation and model definitions
Data generation/
    Main_data_generation.m
    get_notochord_params_base.m
    Notochord_main_pkg.m
    all_cell_volume.m
    tension.m
    volume_radius_relation.m
### 2. Data storage
Data_for_plot/
    (All generated .mat files used for plotting)
### 3. Growth simulation (Figure 2)
Fig2_growth/
    Figure2_growth_main.m
    Data_processing.m
    VideoWriting_Growth.m
    plot_config_only.m
    plot_notochord_video.m
    Notochord_main_pkg_long.m
    all_cell_volume.m
    tension.m
    viridis.m
    volume_radius_relation.m
    Volume_all_dpf_exp.mat
    Volume_all_dpf_exp_proc_smooth_*.mat
    Growth pattern.mp4
### 4. Root directory (figure generation + utilities)
Root/
    PART1_BaseResults.m
    PART2_Propagation_and_TimeScales.m
    PART3_AsymmetricResponse.m
    PART4_VolumeChange_DifScenarios.m
    PART5_DecayProfiles.m

    Plot1_dynamical_curve_representation_Intra.m
    Plot2_2D_kymograph_Intra_and_Interstitial.m
    Plot_Decay_Profile.m
    plot_config_only.m
    plot_notochord_video.m
    quantify_prop.m

    Fig_Adjustment.m
    VideoWriting.m
    viridis.m

    all_cell_volume.m
    tension.m
    volume_radius_relation.m

## Figure–Code Mapping
### Figure 2 (Growth dynamics)
•	Fig2_growth/Figure2_growth_main.m 
### Figure 3 and Figure S3 (Propagation dynamics and local response)
•	PART2_Propagation_and_TimeScales.m 
o	Fig. 3H–L 
o	Fig. S3 
### Figure 4 and Figure S4 (Curvature dependence)
•	PART3_AsymmetricResponse.m 
o	Fig. 4C–H 
o	Fig. S4 
### Figure 5A–B (Response asymmetry)
•	PART3_AsymmetricResponse.m 
### Figure 5D (Decay length scale)
•	PART5_DecayProfiles.m 
### Figure 5E (Volume loss vs pressure differential)
•	PART4_VolumeChange_DifScenarios.m 
### Figure S5 (Parameter dependence: decay, permeability, initial conditions)
•	PART5_DecayProfiles.m 
### Base simulation results (used across multiple figures)
•	PART1_BaseResults.m 
o	Provides baseline data and comparison panels

## Key Utility Functions
•	Notochord_main_pkg.m
Core simulation solver 
•	Plot1_dynamical_curve_representation_Intra.m
Line plots of intracellular dynamics 
•	Plot2_2D_kymograph_Intra_and_Interstitial.m
Kymographs of spatial-temporal dynamics 
•	Plot_Decay_Profile.m
Spatial decay analysis 
•	quantify_prop.m
Quantification of propagation dynamics 
•	Fig_Adjustment.m
Standardized figure formatting 
## Notes
•	All figure scripts assume data is available in Data_for_plot/. 
•	Simulation outputs follow a consistent structure: 
o	Geometry → volume, radius, tension 
o	Pressure → hydraulic and osmotic pressure 
o	Flux → water and solute transport 
•	Plotting scripts are modular and can be executed independently once data is prepared. 
## Recommended Usage
To fully reproduce all results:
1.	Run Data generation/Main_data_generation.m 
2.	Move outputs to Data_for_plot/ 
3.	Run PART1–PART5 
4.	Run Fig2_growth/Figure2_growth_main.m
