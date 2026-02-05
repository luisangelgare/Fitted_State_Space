# Fitted State-Space Tool (fitss)

This is a data acquisition-based methodology to construct a state-space model of black-box Voltage Source Converters (VSCs) for application in small-signal stability analysis in modern power-systems using the Component Connection Method (CCM). This approach determines the fitted state-space (FSS) model through the frequency-domain (FD) admittance matrix of the VSC, identified via the input-output responses of voltage and current at the point of connection (POC). By applying a multiple single-input single-output (SISO) vector fitting array, a stable FSS representation is determined with high accuracy. Order reduction methodologies such as Singular Value Decomposition (SVD) and Balanced Realization (BR) are then employed to generate the minimal representation of the VSC, eliminating redundant states and reducing computational burden. 

The **fitss** has been developed as part of the **MSCA-ADOreD** project, funded by the European Union’s Horizon Europe Research and Innovation Programme under the **Marie Skłodowska-Curie grant agreement No. 101073554**.

---

# Validation and Citation Guidelines 
The **fitss** has been validated and documented in the following scientific publication. Users must **cite this work** in any academic, scientific, or industrial document where the tool is used, following the terms of the **GPL-3.0 license** under which the fitss is distributed. 

<div style="background-color:#000000; color:white; padding:18px; border-radius:10px; font-family:monospace; font-size:14px; line-height:1.4;">

<b>BibTeX Citation (SIaD Tool)</b>

@misc{CITCEA-fitted-state-space,<br>
&nbsp;&nbsp;&nbsp;&nbsp;title={Data-Driven State-Space Modeling for Small-Signal Stability Analysis of Black-Box Power Converters},<br>
&nbsp;&nbsp;&nbsp;&nbsp;author={Garcia-Reyes, Luis A. and Prieto-Araujo, Eduardo and Lacerda, Vinicius A. and Arévalo-Soler, Josep and Gomis-Bellmunt, Oriol and Martin-Almenta, Macarena and Nuño-Martínez, Edgar and Renedo, Javier},<br>
&nbsp;&nbsp;&nbsp;&nbsp;year={2025},<br>
&nbsp;&nbsp;&nbsp;&nbsp;booktitle={2025 IEEE PES Innovative Smart Grid Technologies Conference Europe (ISGT Europe)},<br>
&nbsp;&nbsp;&nbsp;&nbsp;pages={1-6},<br>
}
</div>

---

# Network parameters for the example in the paper:

- See the file: "IEEE9_GFL_GFM.xlsx"
- See the file: "New_England.xlsx"

# How to use (example):

- Copy the function files to the common project folder.
- Use the function "fitss2" to determine the fitted state-space model.

# Tutorial:
  
  Check the file "example_fitss.m" and follow the structure of the inputs and outputs (the functions and the script will be updated).


