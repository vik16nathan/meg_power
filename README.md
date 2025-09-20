## Lower Bound Estimates for Electrophysiological Power Dissipation in Human Gray Matter 

This repository contains all code used to run analyses for "Lower Bound Estimates for Electrophysiological Power Dissipation in Human Gray Matter" (Nathan et al., 2025, in review). 

The code is grouped into the following sections:

- _process_brainstorm_data_: used to preprocess support areas and cortical thicknesses used in primary dipole power calculations (fig. 1)
- _calculate_parcellate_primary_p_: used to calculate primary dipole P and parcellate into Schaeffer 600 regions for future analysis (fig. 1)
- _phantom_: physical phantom analyses (fig. 2)
- _cell_res_: dendrite-derived resistance calculations in each parcel of the brain (fig. 2D)
- _secondary_p_: code to simulate secondary currents at every point in an FEM mesh before subsequently calculating power (fig. 3)
- _metabolism_correlations_: code to calculate S600-parcellated power's correlation with glucose/oxygen/biologically-defined metabolism (fig. 4/5)
- _supp_: supplementary figures (pie charts of secondary power consumption by tissue)

