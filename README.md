# MASIMRS User Guide

Perform 2D PRESS MRSI Processing for Philips Enhanced DICOM and Siemens DICOM data with LCModel

## Contents

* [Overview](#overview)
* [Authors and Reference](#authors-and-reference)
* [Containerization of Source Code](#containerization-of-source-code)
* [Command](#command)
* [Arguments and Options](#arguments-and-options)
* [Outputs (Voxel-wise)](#outputs-voxel-wise)
* [Outputs (Regional)](#outputs-regional)

## Overview

MASIMRS performs an analysis of 2D PRESS MRSI from Philips and Siemens DICOM data with LCModel. It converts the complex time domain data stored in the DICOMs to NIFTI files with the signal (or fourier spectrum) stored in the fourth dimension. It then processes the signal in each voxel with LCModel to compute metabolite peaks and saves each metabolite map as its own NIFTI file. If given a segmentation and a label target, MASIMRS will also perform a regional analysis by computing a weighted average regional signal based on the amount of overlap from the MRSI voxels with consideration of lipid intereference. It will then convert this signal to NIFTI format and run it through LCModel to obtain the peaks and ratios.

## Authors and References

[Cai LY](mailto:leon.y.cai@vanderbilt.edu), Tanase C, Anderson AW, Patel NJ, Lee CA, Jones RS, LeStourgeon LM, Mahon A, Taki I, Juvera J, Pruthi S, Gwal K, Ozturk A, Kang H, Rewers A, Rewers MJ, Alonso GT, Glaser N, Ghetti S, Jaser SS, Landman BA, Jordan LC. Exploratory Multisite MR Spectroscopic Imaging Shows White Matter Neuroaxonal Loss Associated with Complications of Type 1 Diabetes in Children. American Journal of Neuroradiology. 2023 Jul 1;44(7):820-7.

Cai LY, Del Tufo SN, Barquero L, D'Archangel M, Sachs L, Cutting LE, Glaser N, Ghetti S, Jaser SS, Anderson AW, Jordan LC, Landman BA. Spatiospectral image processing workflow considerations for advanced MR spectroscopy of the brain. Medical Imaging: Image Processing. 2023 Aug 1. In submission.

[Medical-image Analysis and Statistical Interpretation (MASI) Lab](https://my.vanderbilt.edu/masi), Vanderbilt University, Nashville, TN, USA

## Containerization of Source Code

MASIMRS is designed to run inside a Singularity container, `masimrs.sif`. The following commands will build the most updated container from the source code hosted on GitHub and install the necessary dependencies.

    git clone https://github.com/MASILab/masimrs.git
    cd masimrs
    git checkout v1.0.0
    sudo singularity build /path/to/masimrs.sif singularity

To build the container, we use Singularity 3.8 Community Edition with root permissions.

## Command

To run the MASIMRS software inside the Singularity container, use the following steps.

* First, place your MRS DICOM file, `my.dcm` into your local inputs folder (i.e., `/host/path/to/inputs/my.dcm`).
* Then, create the output folder on your local machine (i.e., `/host/path/to/outputs`).
* Add any additional inputs into your local inputs folder based on the [below options](#arguments-and-options). 

Of note, with this singularity command, the local input and output folder will be bound into `/INPUTS` and `/OUTPUTS` *inside* the container, respectively. All arguments passed into the singularity command should thus use paths relative to the *inside* of the container. For instance, if you want to provide `/home/Downloads/file.nii.gz` as an argument for an option, you must first move it into `/host/path/to/inputs/file.nii.gz` and then pass `/INPUTS/file.nii.gz` to the argument.

Binding the `/tmp` directory is required when running the container with `--contain` which is preferred.

    singularity run 
    --cleanenv 
    --contain
    --bind /host/path/to/inputs:/INPUTS
    --bind /host/path/to/outputs:/OUTPUTS
    --bind /tmp:/tmp
    /host/path/to/masimrs.sif
    /INPUTS/<my.dcm>
    /OUTPUTS/<my_prefix>
    <scanner_type>
    [options]

## Arguments and Options

**/INPUTS/<my.dcm>**

Path to the Philips or Siemens 2D CSI data in enhanced DICOM format relative to the *inside* of the container.

**/OUTPUTS/<my_prefix>**

Path with prefix of the outputs relative to the *inside* of the container.

**<scanner_type>**

A string indicating the scanner type, either `philips` or `siemens`.

**-x** or **--novox**

Do NOT perform a voxel-wise analysis of metabolite ratios.

**-s /seg.nii.gz 1+2,3** or **--seg /seg.nii.gz 1+2,3**

This allows MASIMRS to perform an optional regional analysis and takes two arguments: `/seg.nii.gz` and `1+2,3`. The first, `/seg.nii.gz`, is a path relative to the *inside* of the container to a tissue segmentation NIFTI file co-registered to `/INPUTS/<my.dcm>`. The second is a plus (`+`) or comma (`,`) separated list of labels in `/seg.nii.gz` corresponding to the regions intended to be analyzed. For `1+2,3`, regions 1 and 2 will be considered as one region separately from region 3.

This option can be used multiple times to analyze multiple regions from multiple segmentation files in one command.

**-m /input/brain_mask.nii.gz** or **--mask /input/brain_mask.nii.gz**
                        
This allows MASIMRS to perform lipid exclusion of voxels included in `/input/brain_mask.nii.gz` before regional analysis. The argument should be provided as the path to an image co-registered to `/INPUTS/<my.dcm>` relative to the *inside* of the container. 

Default = do NOT remove lipid interference

**-b /custom/set.basis** or **--basis /custom/set.basis**

This allows the user to provide a custom basis function for LCModel analysis.

Default = use a Provencher basis set with matching TE

**-r REF** or **--ref REF**

A string indicating which LCModel metabolite to use in the denominator when computing ratios. See LCModel documentation for a list.

Default = `Cr+PCr` (creatine and phosphocreatine)

**-g /target.nii.gz interp** or **--grid /target.nii.gz interp**

This allows MASIMRS to regrid MRS signal and spectral maps (computed as the sum of the absolute value across the time or spectral dimension, respectively) to match `/target.nii.gz` which should be co-registered to `/INPUTS/<my.dcm>` and whose path should be provided relative to the *inside* of the container using interpolation method `interp` for visualization. `interp` can be `nearest`, `linear`, or `cubic`.

Default = do NOT regrid.

**-n N** or **--nthreads N**

A positive integer indicating the number of threads to use when running portions of the pipeline that can be multithreaded (i.e. voxel-wise analysis).

Default = `1` (do NOT multithread)

**-v** or **--verbose**

Output status to console as program runs.

**-h** or **--help**

## Outputs (Voxel-wise)

**Raw data formatted as NIFTI**

These outputs are repeated, swapping `met` for `h2o` in the file names, for the water reference for Philips data.

* **`/OUTPUTS/my_prefix`_met_signal_abs.nii.gz**: A 4D NIFTI file with the magnitude of the raw spectroscopy data in the time domain.

* **`/OUTPUTS/my_prefix`_met_signal_ang.nii.gz**: A 4D NIFTI file with the phase of the raw spectroscopy data in the time domain.

* **`/OUTPUTS/my_prefix`_met_signal_real.nii.gz**: A 4D NIFTI file with the real component of the raw spectroscopy data in the time domain.

* **`/OUTPUTS/my_prefix`_met_signal_imag.nii.gz**: A 4D NIFTI file with the imaginary component of the raw spectroscopy data in the time domain.

* **`/OUTPUTS/my_prefix`_met_signal_abs_sum.nii.gz**: `/OUTPUTS/my_prefix_met_signal_abs.nii.gz` summed in the 4th dimension.

* **`/OUTPUTS/my_prefix`_met_signal_abs_sum_grid.nii.gz**: `/OUTPUTS/my_prefix_met_signal_abs_sum.nii.gz` regridded to match the input to `--grid`.

* **`/OUTPUTS/my_prefix`_met.ppm**: A text file listing the PPM of the raw data.

* **`/OUTPUTS/my_prefix`_met_spectra_abs.nii.gz**: A 4D NIFTI file with the magnitude of the raw spectroscopy data in the fourier domain.

* **`/OUTPUTS/my_prefix`_met_spectra_ang.nii.gz**: A 4D NIFTI file with the phase of the raw spectroscopy data in the fourier domain.

* **`/OUTPUTS/my_prefix`_met_spectra_real.nii.gz**: A 4D NIFTI file with the real component of the raw spectroscopy data in the fourier domain.

* **`/OUTPUTS/my_prefix`_met_spectra_imag.nii.gz**: A 4D NIFTI file with the imaginary component of the raw spectroscopy data in the fourier domain.

* **`/OUTPUTS/my_prefix`_met_spectra_abs_sum.nii.gz**: `/OUTPUTS/my_prefix_met_spectra_abs.nii.gz` summed in the 4th dimension.

* **`/OUTPUTS/my_prefix`_met_spectra_abs_sum_grid.nii.gz**: `/OUTPUTS/my_prefix_met_spectra_abs_sum.nii.gz` regridded to match the input to `--grid`.

**LCModel inputs**

* **`/OUTPUTS/my_prefix`.ctrl**: A summary of the LCModel control parameters used in the voxel-wise analysis.

* **`/OUTPUTS/my_prefix`_sl1.raw**: The RAW file containing metabolite spectra data for LCModel.

* **`/OUTPUTS/my_prefix`_sl1.h2o**: The RAW file containing water reference spectra data for LCModel (Philips only).

**LCModel spectra outputs**

* **`/OUTPUTS/my_prefix`_lcm.ppm**: A text file listing the PPM of the data processed by LCModel

* **`/OUTPUTS/my_prefix`_lcm_spectra.nii.gz**: A 4D NIFTI file with the real spectra data converted by LCModel for each voxel. The fourth dimension contains the spectra sampled at the PPM in `/OUTPUTS/my_prefix_lcm.ppm`.

* **`/OUTPUTS/my_prefix`_lcm_spectra_abs_sum.nii.gz**: The absolute value of `/OUTPUTS/my_prefix_lcm_spectra.nii.gz` summed in the 4th dimension.

* **`/OUTPUTS/my_prefix`_lcm_fit.nii.gz**: A 4D NIFTI file with the spectra data fit by LCModel for each voxel. The fourth dimension contains the spectra sampled at the PPM in `/OUTPUTS/my_prefix_lcm.ppm`.

* **`/OUTPUTS/my_prefix`_lcm_baseline.nii.gz**: A 4D NIFTI file with the spectra baseline fit by LCModel for each voxel. The fourth dimension contains the baseline sampled at the PPM in `/OUTPUTS/my_prefix_lcm.ppm`.

**LCModel metabolite outputs**

* **`MET`** is variable and denotes the labels of metabolites fit with LCModel

* **`/OUTPUTS/my_prefix`_`MET`.nii.gz**: 3D NIFTI images generated during voxel-wise analysis, one for each metabolite peak fit by LCModel, indicating the peak height computed by LCModel.

* **`/OUTPUTS/my_prefix`_`MET`_SD.nii.gz**: 3D NIFTI images generated during voxel-wise analysis, one for each metabolite peak fit by LCModel, indicating the %SD metric computed by LCModel.

* **`/OUTPUTS/my_prefix`.pdf**: A PDF document showing the voxel-wise LCModel outputs summarizing spectra, fit, and ratios.

## Outputs (Regional)

* **`SEG`** is variable and denotes the file name without extension passed into `--seg`
* **`LBL`** is variable and denotes the string of labels passed into `--seg`
* **`REF`** is variable and denotes the string passed into `--ref`

**Raw data formatted as text files**

These outputs are repeated once for each designated region.

* **`/OUTPUTS/my_prefix`\_`SEG`\_`LBL`.signal**: A text file containing the complex time domain signal of the designated region.

* **`/OUTPUTS/my_prefix`\_`SEG`\_`LBL`.spectra**: A text file containing the complex fourier domain signal of the designated region.

* **`/OUTPUTS/my_prefix`\_`SEG`\_`LBL`.ppm**: A text file containing the PPM of the spectra.

**Lipid exclusion**

* **`/OUTPUTS/my_prefix`\_`outliers`.nii.gz**: A 3D NIFTI file masking the MRS voxels excluded per `--mask`.

**Intermediates during weighted average pooling**

These outputs are repeated once for each designated region.

* **`/OUTPUTS/my_prefix`\_`SEG`\_`LBL`_mask.nii.gz**: A 3D NIFTI file masking the designated region. 

* **`/OUTPUTS/my_prefix`\_`SEG`\_`LBL`_weights.nii.gz**: A 3D NIFTI file indicating the contributions of each voxel to the weighted average.

**LCModel inputs**

These outputs are repeated once for each designated region.

* **`/OUTPUTS/my_prefix`\_`SEG`\_`LBL`.ctrl**: A summary of the LCModel control parameters used in the analysis for the designated region.

* **`/OUTPUTS/my_prefix`\_`SEG`\_`LBL`_sl1.raw**: The RAW file containing metabolite spectra data for LCModel for the designated region.

* **`/OUTPUTS/my_prefix`\_`SEG`\_`LBL`_sl1.h2o**: The RAW file containing water reference spectra data for LCModel for the designated region (Philips only).

**LCModel outputs**

These outputs are repeated once for each designated region.

* **`/OUTPUTS/my_prefix`\_`SEG`\_`LBL`.lcm_ppm**: A text file listing the PPM of the data processed by LCModel for the designated region

* **`/OUTPUTS/my_prefix`\_`SEG`\_`LBL`.lcm_spectra**: A text file containing the spectra data converted by LCModel corresponding to the PPM in `/OUTPUTS/my_prefix_SEG_LBL.lcm_ppm` for the designated region.

* **`/OUTPUTS/my_prefix`\_`SEG`\_`LBL`.lcm_fit**: A text file containing the spectra data fit by LCModel corresponding to the PPM in `/OUTPUTS/my_prefix_SEG_LBL.lcm_ppm` for the designated region.

* **`/OUTPUTS/my_prefix`\_`SEG`\_`LBL`.lcm_baseline**: A text file containing the spectra baseline fit by LCModel corresponding to the PPM in `/OUTPUTS/my_prefix_SEG_LBL.lcm_ppm` for the designated region.

* **`/OUTPUTS/my_prefix`\_`SEG`\_`LBL`.pdf**: LCModel print-out summarizing the spectra, fit, and ratios for the designated region.

**Metabolite output**

* **`/OUTPUTS/my_prefix`\_`REF`\_ratios_by_label.csv**: A CSV file indicating the metabolite peaks, ratios, and %SD for all regions designated.
