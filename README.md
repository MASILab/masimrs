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

MASIMRS performs an analysis of 2D PRESS MRSI from Philips and Siemens DICOM data with LCModel. It converts the complex time domain data stored in the DICOMs to NIFTI files with the signal (or fourier spectrum) stored in the fourth dimension. It then processes the signal in each voxel with LCModel to compute metabolite peaks and saves each metabolite map as its own NIFTI file. If given a segmentation and a label target, MASIMRS will also perform a regional analysis by computing a weighted average regional signal based on the amount of overlap from the MRSI voxels. It will then convert this signal to NIFTI format and run it through LCModel to obtain the peaks and ratios.

## Authors and Reference

TBD

[Medical-image Analysis and Statistical Interpretation (MASI) Lab](https://my.vanderbilt.edu/masi), Vanderbilt University, Nashville, TN, USA

## Containerization of Source Code

MASIMRS is designed to run inside a Singularity container, `masimrs.sif`. The following commands will build the most updated container from the source code hosted on GitHub and install the necessary dependencies.

    git clone https://github.com/MASILab/masimrs.git
    cd masimrs
    git checkout v1.0.0
    sudo singularity build /path/to/masimrs.sif singularity

To build the container, we use Singularity 3.8 Community Edition with root permissions.

## Command

To run the MASIMRS software inside the Singularity container, use the following command. Please see [below](#arguments-and-options) for more elaboration of the arguments and options.

Of note, the paths input into the script should be designed relative to the *inside* of the container and files should be bound into the container: an `/INPUTS` and `/OUTPUTS` folder are provided for you. For instance, if you have an input DICOM MRSI file at `/home/Desktop/my.dcm`, it should be bound into the `/INPUTS` folder in the container with the Singularity option `--bind /home/Desktop/my.dcm:/INPUTS/my.dcm` and given as an argument to the container as `/INPUTS/my.dcm`. Similarly, if you want to save outputs with the file prefix `/home/Downloads/my_prefix`, bind `/home/Downloads/` to the `/OUTPUTS` directory in the container and provide the output prefix argument as `/OUTPUTS/my_prefix`. Binding the `/tmp` directory is required when running the container with `--contain` which is preferred.

    singularity run 
    --cleanenv 
    --contain
    --bind /host/path/to/my.dcm:/INPUTS/my.dcm
    --bind /host/path/to/output/directory/:/OUTPUTS
    --bind /tmp:/tmp
    /host/path/to/masimrs.sif
    /INPUTS/my.dcm
    /OUTPUTS/my_prefix
    scanner_type
    [options]

## Arguments and Options

**/INPUTS/my.dcm**

Path to the Philips or Siemens 2D CSI data in enhancd DICOM format relative to the *inside* of the container.

**/OUTPUTS/my_prefix**

Path with prefix of the outputs relative to the *inside* of the container.

**scanner_type**

A string indicating the scanner type, either `philips` or `siemens`.

**-x** or **--novox**

Do NOT perform a voxel-wise analysis of metabolite ratios.

**-s /seg.nii.gz 1+2,3** or **--seg /seg.nii.gz 1+2,3**

This allows MASIMRS to perform an optional regional analysis and takes two arguments: `/seg.nii.gz` and `1+2,3`. The first, `/seg.nii.gz`, is a path relative to the *inside* of the container to a tissue segmentation NIFTI file co-registered to `/INPUTS/my.dcm`. The second is a plus (`+`) or comma (`,`) separated list of labels in `/seg.nii.gz` corresponding to the regions intended to be analyzed. For `1+2,3`, regions 1 and 2 will be considered as one region separately from region 3.

This option can be used multiple times to analyze multiple regions from multiple segmentation files in one command.

**-r REF** or **--ref REF**

A string indicating which LCModel metabolite to use in the denominator when computing ratios. See LCModel documentation for a list.

Default = `Cr+PCr` (creatine and phosphocreatine)

**-n N** or **--nthreads N**

A positive integer indicating the number of threads to use when running portions of the pipeline that can be multithreaded (i.e. voxel-wise analysis).

Default = `1` (do NOT multithread)

**-v** or **--verbose**

Output status to console as program runs.

**-h** or **--help**

## Outputs (Voxel-wise)

* **`/OUTPUTS/my_prefix`** is variable and denotes the prefix output argument
* **`MET`** is variable and denotes the labels of metabolites fit with LCModel

**Raw data formatted as NIFTI**

These outputs are repeated, swapping `met` for `h2o` in the file names, for the water reference for Philips data.

* **`/OUTPUTS/my_prefix`_met_signal_abs.nii.gz**: A 4D NIFTI file with the magnitude of the raw spectroscopy data in the time domain.

* **`/OUTPUTS/my_prefix`_met_signal_ang.nii.gz**: A 4D NIFTI file with the phase of the raw spectroscopy data in the time domain.

* **`/OUTPUTS/my_prefix`_met_signal_real.nii.gz**: A 4D NIFTI file with the real component of the raw spectroscopy data in the time domain.

* **`/OUTPUTS/my_prefix`_met_signal_imag.nii.gz**: A 4D NIFTI file with the imaginary component of the raw spectroscopy data in the time domain.

* **`/OUTPUTS/my_prefix`_met_signal_abs_sum.nii.gz**: `/OUTPUTS/my_prefix_met_signal_abs.nii.gz` summed in the 4th dimension.

* **`/OUTPUTS/my_prefix`_met.ppm**: A text file listing the PPM of the raw data.

* **`/OUTPUTS/my_prefix`_met_spectra_abs.nii.gz**: A 4D NIFTI file with the magnitude of the raw spectroscopy data in the fourier domain.

* **`/OUTPUTS/my_prefix`_met_spectra_ang.nii.gz**: A 4D NIFTI file with the phase of the raw spectroscopy data in the fourier domain.

* **`/OUTPUTS/my_prefix`_met_spectra_real.nii.gz**: A 4D NIFTI file with the real component of the raw spectroscopy data in the fourier domain.

* **`/OUTPUTS/my_prefix`_met_spectra_imag.nii.gz**: A 4D NIFTI file with the imaginary component of the raw spectroscopy data in the fourier domain.

* **`/OUTPUTS/my_prefix`_met_spectra_abs_sum.nii.gz**: `/OUTPUTS/my_prefix_met_spectra_abs.nii.gz` summed in the 4th dimension.

**LCModel inputs**

* **`/OUTPUTS/my_prefix`.ctrl**: A summary of the LCModel control parameters used in the voxel-wise analysis.

* **`/OUTPUTS/my_prefix`_sl1.raw**: The RAW file containing metabolite spectra data for LCModel.

* **`/OUTPUTS/my_prefix`_sl1.h2o**: The RAW file containing water reference spectra data for LCModel (Philips only).

**LCModel spectra outputs**

* **`/OUTPUTS/my_prefix`_lcm.ppm**: A text file listing the PPM of the data processed by LCModel

* **`/OUTPUTS/my_prefix`_lcm_spectra.nii.gz**: A 4D NIFTI file with the spectra data converted by LCModel for each voxel. The fourth dimension contains the spectra sampled at the PPM in `/OUTPUTS/my_prefix_lcm.ppm`.

* **`/OUTPUTS/my_prefix`_lcm_fit.nii.gz**: A 4D NIFTI file with the spectra data fit by LCModel for each voxel. The fourth dimension contains the spectra sampled at the PPM in `/OUTPUTS/my_prefix_lcm.ppm`.

* **`/OUTPUTS/my_prefix`_lcm_baseline.nii.gz**: A 4D NIFTI file with the spectra baseline fit by LCModel for each voxel. The fourth dimension contains the baseline sampled at the PPM in `/OUTPUTS/my_prefix_lcm.ppm`.

**LCModel metabolite outputs**

* **`/OUTPUTS/my_prefix`_`MET`.nii.gz**: 3D NIFTI images generated during voxel-wise analysis, one for each metabolite peak fit by LCModel, indicating the peak height computed by LCModel.

* **`/OUTPUTS/my_prefix`_`MET`_SD.nii.gz**: 3D NIFTI images generated during voxel-wise analysis, one for each metabolite peak fit by LCModel, indicating the %SD metric computed by LCModel.

* **`/OUTPUTS/my_prefix`.pdf**: A PDF document showing the voxel-wise LCModel outputs summarizing spectra, fit, and ratios.

## Outputs (Regional)

* **`seg`** is variable and denotes the file name without extension passed into `--seg`
* **`lbl`** is variable and denotes the string of labels passed into `--seg`
* **`REF`** is variable and denotes the string passed into `--ref`

**Raw data formatted as text files**

These outputs are repeated once for each designated region.

* **`/OUTPUTS/my_prefix`\_`seg`\_`lbl`.signal**: A text file containing the complex time domain signal of the designated region.

* **`/OUTPUTS/my_prefix`\_`seg`\_`lbl`.spectra**: A text file containing the complex fourier domain signal of the designated region.

* **`/OUTPUTS/my_prefix`\_`seg`\_`lbl`.ppm**: A text file containing the PPM of the spectra.

**Intermediates during weighted average pooling**

These outputs are repeated once for each designated region.

* **`/OUTPUTS/my_prefix`\_`seg`\_`lbl`_mask.nii.gz**: A 3D NIFTI file masking the designated region. 

* **`/OUTPUTS/my_prefix`\_`seg`\_`lbl`_weights.nii.gz**: A 3D NIFTI file indicating the contributions of each voxel to the weighted average.

**LCModel inputs**

These outputs are repeated once for each designated region.

* **`/OUTPUTS/my_prefix`\_`seg`\_`lbl`.ctrl**: A summary of the LCModel control parameters used in the analysis for the designated region.

* **`/OUTPUTS/my_prefix`\_`seg`\_`lbl`_sl1.raw**: The RAW file containing metabolite spectra data for LCModel for the designated region.

* **`/OUTPUTS/my_prefix`\_`seg`\_`lbl`_sl1.h2o**: The RAW file containing water reference spectra data for LCModel for the designated region (Philips only).

**LCModel outputs**

These outputs are repeated once for each designated region.

* **`/OUTPUTS/my_prefix`\_`seg`\_`lbl`.lcm_ppm**: A text file listing the PPM of the data processed by LCModel for the designated region

* **`/OUTPUTS/my_prefix`\_`seg`\_`lbl`.lcm_spectra**: A text file containing the spectra data converted by LCModel corresponding to the PPM in `/OUTPUTS/my_prefix_seg_lbl.lcm_ppm` for the designated region.

* **`/OUTPUTS/my_prefix`\_`seg`\_`lbl`.lcm_fit**: A text file containing the spectra data fit by LCModel corresponding to the PPM in `/OUTPUTS/my_prefix_seg_lbl.lcm_ppm` for the designated region.

* **`/OUTPUTS/my_prefix`\_`seg`\_`lbl`.lcm_baseline**: A text file containing the spectra baseline fit by LCModel corresponding to the PPM in `/OUTPUTS/my_prefix_seg_lbl.lcm_ppm` for the designated region.

* **`/OUTPUTS/my_prefix`\_`seg`\_`lbl`.pdf**: LCModel print-out summarizing the spectra, fit, and ratios for the designated region.

**Metabolite output**

* **`/OUTPUTS/my_prefix`\_`REF`\_ratios_by_label.csv**: A CSV file indicating the metabolite peaks, ratios, and %SD for all regions designated.
