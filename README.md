# RoSE-O_POLCAM

This repository contains MATLAB scripts and utilities for the analysis of simulated and experimental POLCAM single-molecule orientation-localisation microscopy (SMOLM) data using RoSE-O, providing accurate and precise measurements of the location and orientation of single molecules.

## Associated Paper and Data

Paper: [POLCAM: Instant molecular orientation microscopy for the life sciences](paper link here)

bioRxiv preprint: [POLCAM: Instant molecular orientation microscopy for the life sciences](https://doi.org/10.1101/2023.02.07.527479)

<!-- Data repository: [POLCAM-SMOLM data](link to data here)  -->

## Code Structure

- `utils\`: Utility functions for the RoSE-O analysis, including background estimation, basis image computation, and the core RoSE-O analysis algorithm.
- `RoSEO_PolCam_simData_demo.m`: Simulate POLCAM data using predefined parameters, analyze the data using RoSE-O, and evaluate the detection performance and estimation accuracy and precision.
- `RoSEO_PolCam_fibrils_demo.m`: Analyze experimental POLCAM-SMOLM images of amyloid fibrils using RoSE-O. The raw images required for this analysis are available in the data repository.
- Sample data required to run this code is available from the [Lew Lab OSF repository](https://osf.io/utjmr/):
  - `fibrils_NileRed_0.tif`: https://osf.io/download/65ededc74dfd8c0b0b27fcf2/
  - `fibrils_NileRed_n45.tif`: https://osf.io/download/65eded5a8992ec0b57788b41/
  - `fibrils_NileRed_p45.tif`: https://osf.io/download/65eded8b4dfd8c0b1227f52b/
  - `fibrils_NileRed_90.tif`: https://osf.io/download/65eded87e5e51c0b6ebc5ef6/

RRID:SCR_025340
