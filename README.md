# Fitted-State-Space Tool (fitss2)

This is a data acquisition-based methodology to construct a state-space model of black-box Voltage Source Converters (VSCs) for application in small-signal stability analysis in modern power-systems using the Component Connection Method (CCM). This approach determines the fitted state-space (FSS) model through the frequency-domain (FD) admittance matrix of the VSC, identified via the input-output responses of voltage and current at the point of connection (POC). By applying a multiple single-input single-output (SISO) vector fitting array, a stable FSS representation is determined with high accuracy. Order reduction methodologies such as Singular Value Decomposition (SVD) and Balanced Realization (BR) are then employed to generate the minimal representation of the VSC, eliminating redundant states and reducing computational burden. 

# Announcement:

**Full access to the tool will be provided once the paper is accepted for the ISGT 2025 conference. For now, this repository only contains the IEEE 9-bus system parameters used in the case study submitted to the ISGT 2025 conference.**

# Network parameters:

- See the file: "IEEE9_GFL_GFM.xlsx"

# How to use:

- Copy the function to the common project folder.
- Use the function "fitss2" to determine the fitted state-space model.

# Tutorial:
  
  Check the file "example_fitss.m" and follow the structure of the inputs and outputs (the functions and the script will be updated).

# How to cite:

Garcia-Reyes, Luis A.; Prieto-Araujo, Eduardo; A. Lacerda, Vinícius; Arévalo-Soler, Josep; Gomis-Bellmunt, Oriol; Martin-Almenta, Macarena; Nuño-Martínez, Edgar; Renedo, Javier, "Data-Driven State-Space Modeling for Small-Signal Stability Analysis of Black-Box Power Converters", 2025 IEEE PES Innovative Smart Grid Technologies Europe (ISGT EUROPE), submmited for review. 
