# Fitted-State-Space Tool

This is a methodology for constructing a state-space representation of a black-box-type Voltage Source Converter (VSC) using its input-output frequency response characterization. The proposed approach operates within the synchronous reference frame (dq0), where the VSC’s admittance (for impedance improvements will be come later) is identified in the frequency domain. By applying a method based on the Pole-Collapsing Column Fitting (PCCF), a linear and stable state-space model is obtained. To ensure model accuracy while reducing computational complexity, order-reduction techniques such as Singular Value Decomposition (SVD) and balanced realization are utilized.

# How to use:

- Copy all the functions to the common project folder.
- Use the function "FittedStateSpace2" to determine the state-space model.

# Tutorial:
  
  Run the file "example_fitss.m" and follow the structure of the inputs and outputs.

# How to cite:

Garcia-Reyes L. A., Prieto-Araujo E., A. Lacerda V., "Reduced-Order State-Space Modeling of VSCs based on Frequency-Domain Identification", 2025 IEEE PES Innovative Smart Grid Technologies Europe (ISGT EUROPE),
