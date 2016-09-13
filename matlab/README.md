dti-qc
------

Usage:

    analyze_dti_phantom(dwi, fa, bval, output, ynNyq)

Arguments:

+ `dwi`: 4D diffusion weighted image
+ `fa`: FA map from `DTIfit`
+ `bval`: b value files form `dcm2nii`
+ `output`: folder to store results
+ `ynNyq`: (`'y'`/`'n'`'). `'y'` for ASSET/accelerated data

Details:

+ add this folder to your `MATLABPATH`.
+ Run `analyze_dti_phantom(dwi.nii.gz, fa.nii.gz, .bval, output/, 'y/n')`.
+ Gather outputs (`.csv` and `.jpg`) from `output/`.

phantom analysis code by Sofia Chavez, 2015.

code packaged and maintained by Joseph Viviano, 2016.

