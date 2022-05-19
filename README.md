# MASIMRS User Guide

Perform 2D PRESS MRSI Processing for Philips Enhanced DICOM and Siemens DICOM data with LCModel

## Contents

* [Authors and Reference](#authors-and-reference)
* [Containerization of Source Code](#containerization-of-source-code)
* [Command](#command)
* [Arguments and Options](#arguments-and-options)
* [Outputs](#outputs)

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

Of note, the paths input into the script should be designed relative to the *inside* of the container and files should be bound into the container: an `/INPUTS` and `/OUTPUTS` folder are provided for you. For instance, if you have an input DICOM file at `/home/Desktop/my.dcm`, it should be bound into the `/INPUTS` folder in the container with the Singularity option `--bind /home/Desktop/my.dcm:/INPUTS/my.dcm` and given as an argument to the container as `/INPUTS/my.dcm`. Similarly, if you want to save outputs with the file prefix `/home/Downloads/my_prefix`, bind `/home/Downloads/` to the `/OUTPUTS` directory in the container and provide the output prefix argument as `/OUTPUTS/my_prefix`. Binding the tmp directory is required when running the container with `--contain`.

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

This allows MASIMRS to perform an optional regional analysis and takes two arguments: `/seg.nii.gz` and `1+2,3`. The first, `/seg.nii.gz`, is a path to a tissue segmentation NIFTI file co-registered to `/INPUTS/my.dcm`. The second is a plus (`+`) or comma (`,`) separated list of labels in `/seg.nii.gz` corresponding to the regions intended to be analyzed. For `1+2,3`, regions 1 and 2 will be considered as one region separately from region 3.

**-r REF** or **--ref REF**

A string indicating which LCModel metabolite to use in the denominator when computing ratios. See LCModel documentation for a list.

Default = `Cr+PCr` (creatine and phosphocreatine)

**-n N** or **--nthreads N**

A positive integer indicating the number of threads to use when running portions of the pipeline that can be multithreaded (i.e. voxel-wise analysis).

Default = `1` (do NOT multithread)

**-v** or **--verbose**

Output status to console as program runs.

**-h** or **--help**

## Outputs

**/OUTPUTS/my_prefix.pdf**

Voxel-wise LCModel outputs summarizing spectra, fit, and ratios.

**/OUTPUTS/my_prefix_MET.nii.gz**

3D NIFTI images generated during voxel-wise analysis, one for each metabolite peak fit by LCModel, indicating the peak height computed by LCModel.

**/OUTPUTS/my_prefix_MET_SD.nii.gz**

3D NIFTI images generated during voxel-wise analysis, one for each metabolite peak fit by LCModel, indicating the %SD metric computed by LCModel.

**/OUTPUTS/my_prefix.ctrl**

A summary of the LCModel control parameters used in the voxel-wise analysis.

**/OUTPUTS/my_prefix_lcm.ppm**, **/OUTPUTS/my_prefix_lcm_spectra.nii.gz**, **/OUTPUTS/my_prefix_lcm_fit.nii.gz**, **/OUTPUTS/my_prefix_lcm_baseline.nii.gz**