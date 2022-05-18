# MASI Library for Managing Philips and Siemens MRS DICOM Processing
# Leon Cai
# MASI Lab
# June 3, 2021

#########
# TO DO #
#########
# - Tramp for absolute concentration interpretation?

##########
# Set Up #
##########

import suspect as sus
from spec2nii.dcm2niiOrientation.orientationFuncs import dcm_to_nifti_orientation
import pydicom as pyd
import numpy as np
import nibabel as nib
import nibabel.nicom.dicomwrappers as nnd
import argparse as ap

import os
import subprocess
from multiprocessing import Pool
from struct import unpack
from glob import glob
from datetime import datetime
from io import StringIO

####################
# Useful Variables #
####################

basis_sets_dir = '/lcmodel/basis-sets-3t'
lcmodel_exec = 'lcmodel'
gs_exec = 'gs'

VERBOSE = False

#####################
# Class Definitions #
#####################

class MRSException(Exception):
    pass

###################################
# Load DICOMs into MRSData format #
###################################

def load_dicom_philips(dcm_file):

    """
    Load 2D PRESS sequences in the ~enhanced~ Philips DICOM Magnetic Resonance Spectroscopy format.

    Parameters
    ----------
    dcm_file : str
        The name of the file to load
    
    Returns
    -------
    MRSData
        The loaded data from the file
    """

    # Extract metadata based on suspect.io.dicom.load_dicom()

    dataset = pyd.dcmread(dcm_file)

    sw = dataset[0x0018, 0x9052].value                                  # Spectral width
    dt = 1.0 / sw                                                       # Time between points
    f0 = dataset[0x0018, 0x9098].value                                  # Transmitter frequency
    ppm0 = dataset[0x0018, 0x9053].value                                # Chemical shift reference
                                           
    tes = [float(frame[0x0018, 0x9114][0][0x0018, 0x9082].value) for frame in dataset.PerFrameFunctionalGroupsSequence] # Echo Time
    if all([te == tes[0] for te in tes]):
        te = tes[0] 
    else:
        raise MRSException('The different frames have different TEs. Aborting.')
    
    trs = [float(frame[0x0018, 0x9112][0][0x0018, 0x0080].value) for frame in dataset.PerFrameFunctionalGroupsSequence] # Repetition Time
    if all([tr == trs[0] for tr in trs]):
        tr = trs[0]
    else:
        raise MRSException('The different frames have different TRs. Aborting.')

    # Compute data shape based on suspect.io.dicom.load_dicom()

    rows = int(dataset[0x0028, 0x0010].value)                           # Rows
    cols = int(dataset[0x0028, 0x0011].value)                           # Columns

    frames = int(dataset[0x0028, 0x0008].value)                         # Frames: Only ever seen this as 2. Unknown how to handle when not 2.
    if frames != 2:
        raise MRSException('This script only supports spectroscopic data with two frames, one for water and one for metabolites. Aborting.')
    
    num_point_rows = int(dataset[0x0028, 0x9001].value)                 # Data point rows (number of MRS sample rows): Only ever seen this as 1. Unknown how to handle when not 1.
    if num_point_rows != 1:
        raise MRSException('This script only supports spectroscopic data with one row. Aborting.')
    
    num_point_cols = int(dataset[0x0028, 0x9002].value)                 # Data point columns (number of MRS samples)

    # Extract byte data to complex vector based on MATLAB script provided by Philips

    byte_data = dataset[0x5600, 0x0020].value                           # Spectroscopy data in bytes
    num_floats = int(len(byte_data)/4)                                  # Convert to complex floats
    raw_data = np.zeros(num_floats)
    for i in range(num_floats):
        raw_data[i] = unpack('<f', byte_data[4*i:4*(i+1)])[0]
    complex_data_1d = raw_data[0::2] - 1j * raw_data[1::2]              # All dimensions compressed, should be divisible by all dimensions in data_shape

    # Note: Philips DICOMs store the real and imaginary parts of the data such that they need to be recombined as a-bi (as opposed to a+bi). LCModel deals with this 
    #       with the BRUKER = T line in the RAW file which tells LCModel to complex conjugate the time data. However, in addition to running LCModel, we save NIFTIs
    #       with this library. So, we'll keep BRUKER = F in the RAW file and conjugate the Philips data ourselves (as above).
    
    # Check Shape

    data_shape = (frames, rows, cols, num_point_rows, num_point_cols)   # Final Shape
    if len(complex_data_1d) / frames / rows / cols / num_point_rows / num_point_cols != 1:
        raise MRSException('Complex data of length {} extracted from DICOM cannot be reshaped into an array of shape {}. Aborting.'.format(len(complex_data_1d), data_shape))

    # Reshape data

    complex_data_2d = np.reshape(complex_data_1d, (frames*rows*cols, num_point_cols))   # (Frames, rows, and columns) by samples
    complex_data_3d = np.stack((complex_data_2d[:rows*cols, :],                         # Frames by (rows and columns) by samples. Frame 0 = metabolites, frame 1 = water reference
                                complex_data_2d[rows*cols:, :]), axis=0)
    complex_data_4d = np.empty((frames, rows, cols, num_point_cols)).astype(complex)    # Frames by rows by columns by samples
    for i in range(num_point_cols):
        complex_data_4d[0, :, :, i] = np.reshape(complex_data_3d[0, :, i], (rows, cols))
        complex_data_4d[1, :, :, i] = np.reshape(complex_data_3d[1, :, i], (rows, cols))
    complex_data_5d = np.expand_dims(complex_data_4d, axis=3)                           # Frames (2) by rows by columns by data point rows (1) by data point columns (num samples)

    # Get affine transformation

    aff = _get_affine(dcm_file, 'philips')

    # Get title metadata

    title = _get_title(dcm_file)

    # Return suspect MRSData object

    return sus.MRSData(complex_data_5d, dt, f0, ppm0=ppm0, te=te, tr=tr, transform=aff, metadata=title)

def load_dicom_siemens(dcm_file):

    """
    Load 2D PRESS sequences in the Siemens DICOM Magnetic Resonance Spectroscopy format.

    Parameters
    ----------
    dcm_file : str
        The name of the file to load
    
    Returns
    -------
    MRSData
        The loaded data from the file
    """

    data = sus.io.load_siemens_dicom(dcm_file)
    if len(data.shape) == 3:
        _print('{} has data with dimension {}. We require four dimensional data (row, col, slice, spec). Adding slice dimension to data.'.format(dcm_file, data.shape))
        data = np.expand_dims(data, axis=2)
        _print('{} now has shape {}.'.format(dcm_file, data.shape))
    data.transform = _get_affine(dcm_file, 'siemens')
    data.metadata = _get_title(dcm_file)
    return data

########################################################
# Convert MRSData to Real and Imaginary Nibabel Images #
########################################################

def nifti_image_philips(data):

    """
    Convert complex Philips 2D PRESS sequences in the suspect MRSData format to the real and imaginary images in Nibabel Nifti1Image format

    Parameters
    ----------
    data : MRSData
        Data in MRSData format
    
    Returns
    -------
    Nibabel.Nifti1Image
        Data in Nibabel format
    """

    re_nii, im_nii = _reim_nii(data[0, :, :, :, :])
    return re_nii, im_nii

def nifti_image_siemens(data):

    """
    Convert complex Siemens 2D PRESS sequences in the suspect MRSData format to the real and imaginary images in Nibabel Nifti1Image format

    Parameters
    ----------
    data : MRSData
        Data in MRSData format
    
    Returns
    -------
    Nibabel.Nifti1Image
        Data in Nibabel format
    """    

    re_nii, im_nii = _reim_nii(data)
    return re_nii, im_nii

############################################
# Save Frequency MRSData into NIFTI format #
############################################

# Note: The Philips spectra seem to come out backwards, even when the time data are 
#       in the same direction (i.e., decaying) as the Siemens data. This doesn't seem
#       to be a problem for LCModel based on QA of LCModel PDFs: it takes the time 
#       domain data directly and seems to sort it out during the basis fitting process. 
#       But, when computing the spectra and looking at them directly w.r.t. the PPM 
#       computed by Suspect (descending), the Philips spectra need to be flipped.
#       
#       Update, this is because Philips stores the data as a-bi, not a+bi and X(-w) = conj(x(t)).
#       We account for this now in the Philips loader function.
#
#       Update 2, turns out we already accounted for this in LCModel with the BRUKER = T line
#       in the RAW file, so updating the Philips loader function double accounted for it and 
#       canceled it. So in order to keep the NIFTI and LCModel data the same, we're going to
#       keep BRUKER = F and do the conjugation ourselves.

def save_nifti_philips(nii_dir, nii_prefix, data):

    mrs_data = data[0, :, :, :, :] # 0 = metabolite, 1 = water; treat "data rows = 1" dim as 1 slice (r, c, s, freq)
    h2o_data = data[1, :, :, :, :]

    _write_data_nii(nii_dir, '{}_met'.format(nii_prefix), mrs_data, data.transform, data.frequency_axis_ppm())
    _write_data_nii(nii_dir, '{}_h2o'.format(nii_prefix), h2o_data, data.transform, data.frequency_axis_ppm())

def save_nifti_siemens(nii_dir, nii_prefix, data):

    _write_data_nii(nii_dir, '{}_met'.format(nii_prefix), data, data.transform, data.frequency_axis_ppm())

####################################
# Save MRSData into LCModel format #
####################################

def save_lcmodel_philips(lcm_dir, lcm_prefix, data):

    """
    Save a Philips MRSData object from load_dicom_philips() for LCModel.

    Parameters
    ----------
    lcm_dir : str
        The directory where the LCModel files should be saved
    lcm_prefix : str
        The root of the file name
    data : MRSData object
        The output of load_dicom_philips()
    
    Returns
    -------
    lcm_files : dict
        Dictionary of LCModel files/strings and list of strings (control files only)
    """

    # LCModel expects only one spectral dimension and the rest as spatial
    # and the water reference and metabolite frames separately

    mrs_data = data[0, :, :, :, :] # Separate frames, treat "spectral rows = 1" as slice dim
    h2o_data = data[1, :, :, :, :]
    basis_file = _get_basis(data.te)
    return _write_lcmodel_files(lcm_dir, lcm_prefix, basis_file, mrs_data, h2o_data=h2o_data, conj=True)

def save_lcmodel_siemens(lcm_dir, lcm_prefix, data):

    """
    Save a Siemens MRSData object from load_dicom_siemens() for LCModel.

    Parameters
    ----------
    lcm_dir : str
        The directory where the LCModel files should be saved
    lcm_prefix : str
        The root of the file name
    data : MRSData object
        The output of load_dicom_siemens()

    Returns
    -------
    lcm_files : dict
        Dictionary of LCModel files/strings and list of strings (control files only)
    """

    # Unlike Philips, Siemens data doesn't come with a water reference

    basis_file = _get_basis(data.te)
    return _write_lcmodel_files(lcm_dir, lcm_prefix, basis_file, data, conj=False)

###############
# Run LCModel #
###############

def run_lcmodel(control_file):

    """
    Run LCModel on a control file as saved by save_lcmodel_*().

    Parameters
    ----------
    control_file : str
        The control file to run, as output by save_lcmodel_*().
    """

    lcmodel_cmd = '{} < {}'.format(lcmodel_exec, control_file)
    _run(lcmodel_cmd)

###################
# Post-processing #
###################

def ps2pdf(ps_files):

    """
    Wrapper for ghostscript to merge PS files from LCModel into a PDF in increments of 500 pages

    Parameters
    ----------
    ps_files : list
        List of strings

    Returns
    -------
    pdf_file : str
        Output PDF file
    """

    # Extract PS prefix, considering edge case with only one PS file

    pdf_file = '{}.pdf'.format(ps_files[0].split('_sl' if len(ps_files) > 1 else '.ps')[0])

    # Break PS files into intervals (GS can handle up to 505 stitches in one doc)

    intervals = list(range(0, len(ps_files)+500, 500))
    intervals[-1] = len(ps_files)

    # For each interval, make a PDF

    interval_pdf_files = []
    for i in range(1, len(intervals)):
        interval_pdf_file = pdf_file.replace('.pdf', '_{}.pdf'.format(i))
        interval_ps_files_str = ' '.join(ps_files[intervals[i-1]:intervals[i]])
        gs_cmd = '{} -q -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -o {} {}'.format(gs_exec, interval_pdf_file, interval_ps_files_str)
        _run(gs_cmd)
        interval_pdf_files.append(interval_pdf_file)

    # Merge interval PDF files into one

    interval_pdf_files_str = ' '.join(interval_pdf_files)
    gs_cmd = '{} -q -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -o {} {}'.format(gs_exec, pdf_file, interval_pdf_files_str)
    _run(gs_cmd)

    # Remove PS Files and interval PDF files

    ps_files_str = ' '.join(ps_files)
    rm_cmd = 'rm {}'.format(ps_files_str)
    _run(rm_cmd)

    rm_cmd = 'rm {}'.format(interval_pdf_files_str)
    _run(rm_cmd)

    # Return PDF file

    return pdf_file

def csv2nii(csv_files, aff):
      
    # Get set up

    dims = _lcm2idxs(csv_files[-1], 'csv') + 1                                                      # Extract spatial dimensions ("shape", index + 1) from last CSV file (r, c, s)
    
    nii_prefix = csv_files[0].split('_sl')[0]
    nii_labels = np.loadtxt(csv_files[0], delimiter=',', skiprows=0, max_rows=1, dtype=str)         # Extract NIFTI (metabolite) labels,
    nii_labels = nii_labels[2:]                                                                     # ...ignoring row and col labels,
    nii_files = []
    nii_isratio = []                                                                                # ...tracking which were ratios,
    for nii_label in nii_labels:                                                                    # ...and reformat them as file names
        nii_isratio.append('/' in nii_label)
        nii_files.append('{}{}.nii.gz'.format(nii_prefix, nii_label.replace(' ', '_')
                                                                   .replace('%', '')
                                                                   .replace('/', '_')))             # *** CAREFUL, THEY'RE NOT ALL /CR or /CR+PCr, even in the same dataset!!! Best to directly compute ratios after the fact.

    img_4d = np.empty((dims[0], dims[1], dims[2], len(nii_labels)))                                 # Build NaN 4D array to hold row, col, slice, met data
    img_4d[:] = np.nan

    # Iterate through each CSV file, and assign values to voxels in 4d image

    for csv_file in csv_files:
        idxs = _lcm2idxs(csv_file, 'csv')
        vals = np.loadtxt(csv_file, delimiter=',', dtype=str)
        if len(vals) != 0:                                                                          # If there is data...
            img_4d[idxs[0], idxs[1], idxs[2], :] = vals[1, 2:].astype(float)                        # ...skip the label row and first two columns (row/col label)

    # Iterate through fourth dimension and save each 3D metabolite and %SD volume as its own NII

    for nii_idx, nii_file in enumerate(nii_files):
        if not nii_isratio[nii_idx]:
            img = img_4d[:, :, :, nii_idx]
            nii = _nii(img, aff)
            nib.save(nii, nii_file)

    # Remove CSV Files

    rm_cmd = 'rm {}'.format(' '.join(csv_files))
    _run(rm_cmd)

    return nii_files

def coo2txt(coo_files, dummy=False):

    assert len(coo_files) == 1, 'COO files can only be converted to text files if only one spectra was run through LCModel. Aborting.'

    # Reformat COO files as Numpy arrays

    ppm, data, fit, bg = _coo2np(coo_files)
    if dummy:
        data[:] = np.nan
        fit[:] = np.nan
        bg[:] = np.nan
        
    # Build and save text file

    txt_prefix = coo_files[-1].split('.coord')[0]

    ppm_file = '{}.lcm_ppm'.format(txt_prefix)
    np.savetxt(ppm_file, ppm)

    data_file = '{}.lcm_spectra'.format(txt_prefix)
    np.savetxt(data_file, data)

    fit_file = '{}.lcm_fit'.format(txt_prefix)
    np.savetxt(fit_file, fit)

    bg_file = '{}.lcm_baseline'.format(txt_prefix)
    np.savetxt(bg_file, bg)

    # Clean Up

    rm_cmd = 'rm {}'.format(' '.join(coo_files))
    _run(rm_cmd)    

    return ppm_file, data_file, fit_file, bg_file

def coo2nii(coo_files, aff):

    assert len(coo_files) > 1, 'COO files can only be converted to NIFTI files if more than one spectra was run through LCModel. Aborting.'

    # Reformat COO files as Numpy arrays

    ppm, data_img, fit_img, bg_img = _coo2np(coo_files)

    # Compute total "energy"

    sum_img = np.sum(data_img, axis=3)

    # Build and Save NIIs

    nii_prefix = coo_files[-1].split('_sl')[0]

    ppm_file = '{}_lcm.ppm'.format(nii_prefix)
    np.savetxt(ppm_file, ppm)

    data_nii = _nii(data_img, aff)
    data_nii_file = '{}_lcm_spectra.nii.gz'.format(nii_prefix)
    nib.save(data_nii, data_nii_file)

    sum_nii = _nii(sum_img, aff)
    sum_nii_file = '{}_lcm_spectra_sum.nii.gz'.format(nii_prefix)
    nib.save(sum_nii, sum_nii_file)

    fit_nii = _nii(fit_img, aff)
    fit_nii_file = '{}_lcm_fit.nii.gz'.format(nii_prefix)
    nib.save(fit_nii, fit_nii_file)

    bg_nii = _nii(bg_img, aff)
    bg_nii_file = '{}_lcm_baseline.nii.gz'.format(nii_prefix)
    nib.save(bg_nii, bg_nii_file)

    # Clean Up

    rm_cmd = 'rm {}'.format(' '.join(coo_files))
    _run(rm_cmd)    

    return ppm_file, data_nii_file, fit_nii_file, bg_nii_file

def ctrl2txt(ctrl_files):

    # Go through all control files and save unique fields into dictionary

    ctrl_dict = {}
    for ctrl_file in ctrl_files:
        ctrl_file_idx_str = ctrl_file.split('.control')[0].split('_')[-1]
        with open(ctrl_file, 'r') as ctrl_fobj:
            ctrl_lines = ctrl_fobj.read().splitlines()
            ctrl_lines = [ctrl_line.replace('\'', '').split(' = ') for ctrl_line in ctrl_lines]
            for ctrl_line in ctrl_lines:
                if len(ctrl_line) == 2:
                    ctrl_key = ctrl_line[0].replace(' ', '')
                    if ctrl_key == 'FILCSV':
                        ctrl_val = ctrl_line[1].replace(ctrl_file_idx_str, 'ROW-COL')
                    else:
                        ctrl_val = ctrl_line[1]
                    if not ctrl_key in ctrl_dict:
                        ctrl_dict[ctrl_key] = [ctrl_val]
                    elif not ctrl_val in ctrl_dict[ctrl_key]:
                        ctrl_dict[ctrl_key].append(ctrl_val)

    # Convert dictionary into text and save control summary

    ctrl_lines = []
    for ctrl_key in ctrl_dict:
        ctrl_line = '{} = {}'.format(ctrl_key, ', '.join(ctrl_dict[ctrl_key]))
        ctrl_lines.append(ctrl_line)
    ctrl_str = '\n'.join(ctrl_lines)
    
    txt_prefix = ctrl_files[0].split('_sl')[0]
    txt_file = '{}.ctrl'.format(txt_prefix)
    with open(txt_file, 'w') as txt_fobj:
        txt_fobj.write('NOTE: THIS FILE IS A SUMMARY OF LCMODEL CONTROL PARAMETERS AND IS NOT DIRECTLY COMPATIBLE WITH LCMODEL.\n\n')
        txt_fobj.write(ctrl_str)
    
    # Clean up control files

    rm_cmd = 'rm {}'.format(' '.join(ctrl_files))
    _run(rm_cmd)

    return txt_file

####################
# Helper Functions #
####################

def _get_basis(te):

    basis_file_dict = {
        20: 'gamma_press_te20_3t_v1.basis',
        40: 'gamma_press_te40_3t_v1.basis',
        68: 'gamma_press_te68_3t_v1.basis',
        80: 'gamma_press_te80_3t_v1.basis',
        144: 'gamma_press_te144_3t_v1.basis',
        270: 'gamma_press_te270_3t_v1.basis',
        288: 'gamma_press_te288_3t_v1.basis',
        25: 'press_te25_3t_v3.basis',
        30: 'press_te30_3t_v3.basis',
        35: 'press_te35_3t_v3.basis',
        135: 'press_te135_3t_v3.basis'
    }
    basis_file_key = np.round(te).astype(int)
    if basis_file_key not in basis_file_dict:
        raise MRSException('No basis file found for a TE of {}. Aborting.'.format(basis_file_key))
        
    return os.path.join(basis_sets_dir, basis_file_dict[basis_file_key])

def _get_affine(dcm_file, scanner):

    if scanner == 'philips': # Philips Enhanced DICOM

        # Load
        dcm = pyd.dcmread(dcm_file)
        frames = dcm.PerFrameFunctionalGroupsSequence

        # Orientation + Position
        pt_orientation = np.array(frames[0].PlaneOrientationSequence[0].ImageOrientationPatient).reshape(2, 3)
        pt_position = np.array(frames[0].PlanePositionSequence[0].ImagePositionPatient)

        # Voxel Size
        px_size = np.array(frames[0].PixelMeasuresSequence[0].PixelSpacing)
        slice_thickness = float(frames[0].PixelMeasuresSequence[0].SliceThickness)
        vox_size = np.array([px_size[0], px_size[1], slice_thickness])
        
        # Spatial Dimensions
        spatial_rows = int(dcm[0x0028, 0x0010].value)
        spatial_cols = int(dcm[0x0028, 0x0011].value)
        spatial_slices = 1 # 2D PRESS
        spatial_dims = np.array([spatial_rows, spatial_cols, spatial_slices])

        # Half Shift
        half_shift = False

    elif scanner == 'siemens_svs': # Siemens Single Voxel Spectroscopy, keeping for posterity

        # Load
        dcm = nnd.wrapper_from_file(dcm_file)

        # Orientation + Position
        pt_orientation = np.array(dcm.csa_header['tags']['ImageOrientationPatient']['items']).reshape(2, 3)
        pt_position = dcm.csa_header['tags']['VoiPosition']['items'] # VoiPosition - this does not have the FOV shift that imagePositionPatient has

        # Voxel Size
        vox_size = np.array([dcm.csa_header['tags']['VoiPhaseFoV']['items'][0],
                             dcm.csa_header['tags']['VoiReadoutFoV']['items'][0],
                             dcm.csa_header['tags']['VoiThickness']['items'][0]])
        
        # Spatial Dimensions
        spatial_dims = np.ones((3,)) # Single voxel

        # Half Shift
        half_shift = False

    elif scanner == 'siemens': # Siemens Chemical Shift Imaging

        # Load
        dcm = nnd.wrapper_from_file(dcm_file)

        # Orientation + Position
        pt_orientation = np.array(dcm.csa_header['tags']['ImageOrientationPatient']['items']).reshape(2, 3)
        pt_position = np.array(dcm.csa_header['tags']['ImagePositionPatient']['items'])

        # Voxel Size
        vox_size = np.array([dcm.csa_header['tags']['PixelSpacing']['items'][0],
                             dcm.csa_header['tags']['PixelSpacing']['items'][1],
                             dcm.csa_header['tags']['SliceThickness']['items'][0]])

        # Spatial Dimensions
        spatial_rows = dcm.csa_header['tags']['Rows']['items'][0]
        spatial_cols = dcm.csa_header['tags']['Columns']['items'][0]
        spatial_slices = 1 # dcm.csa_header['tags']['NumberOfFrames']['items'][0] # 2D PRESS
        spatial_dims = np.array([spatial_rows, spatial_cols, spatial_slices])

        # Half Shift
        half_shift = True

    aff = dcm_to_nifti_orientation(pt_orientation,
                                   pt_position,
                                   vox_size,
                                   spatial_dims,
                                   half_shift=half_shift,
                                   verbose=VERBOSE)

    return aff.Q44

def _get_title(dcm_file):

    dataset = pyd.dcmread(dcm_file)

    date = dataset[0x0008, 0x0021].value
    yyyy = date[:4]
    mm = date[4:6]
    dd = date[6:]

    time = dataset[0x0008, 0x0031].value.split('.')[0]
    hr = time[:2]
    min = time[2:4]
    sec = time[4:]
    
    desc = dataset[0x0008, 0x103e].value
    
    subj = dataset[0x0010, 0x0020].value

    return '{}, scanned on {}/{}/{} at {}:{}:{} ({})'.format(subj, mm, dd, yyyy, hr, min, sec, desc)

def _write_lcmodel_files(lcm_dir, lcm_prefix, basis_file, mrs_data, h2o_data=None, conj=False):

    """
    Write an MRSData object to LCModel files.

    Parameters
    ----------
    lcm_dir : str
        The directory where the LCModel files should be saved
    lcm_prefix : str
        The root of the file name
    basis_file : str
        The basis file to run LCModel with
    mrs_data : MRSData object
        Contains spectra
    h2o_data : MRSData object
        Contains water reference
    conj : bool
        Indicates whether data must be complex conjugated by LCModel (i.e., Bruker, Philips)

    Returns
    -------
    lcm_files : dict
        Dictionary of LCModel files with a string for the basis file and list of strings for input (slice-wise) and output (voxel-wise) files
    """

    # We will run LCModel with 3 spatial dimensions (row, col, slice) and 1 spectral dimension

    if len(mrs_data.shape) != 4:
        raise MRSException("We run LCModel with exactly 3 spatial dimensions (row, col, slice) and 1 spectral dimension. Aborting.")
    spatial_dims = mrs_data.shape[0:-1]

    _print('{} has {} rows, {} columns, and {} slices.'.format(mrs_data.metadata, spatial_dims[0], spatial_dims[1], spatial_dims[2]))

    # Run LCModel to a maximium ppm range of 0.2 to 4.0

    ppm = mrs_data.frequency_axis_ppm()
    ppm_min = np.max((0.2, np.min(ppm)))
    ppm_max = np.min((4.0, np.max(ppm)))

    # Define and write LCModel inputs slice-wise and LCModel controls and outputs voxel-wise,
    # following LCModel naming convention ([prefix]_sl[s]_[r]-[c].[ext]).

    lcm_files = {}
    lcm_files['basis'] = basis_file                                                                     # Initialize lists
    lcm_files['raw'] = []
    lcm_files['h2o'] = []
    lcm_files['ctrl'] = []
    lcm_files['ps'] = []
    lcm_files['csv'] = []
    lcm_files['coo'] = []
    lcm_files['tab'] = []

    for slice_idx in range(1, spatial_dims[2] + 1):                                                     # For each slice...

        lcm_slice_prefix = "{}_sl{}".format(lcm_prefix, slice_idx)                                      # - Create common prefix

        lcm_files['raw'].append(os.path.join(lcm_dir, "{}.raw".format(lcm_slice_prefix)))               # - Add RAW file to list,
        _write_lcmodel_raw(lcm_files['raw'][-1], mrs_data, conj=conj)                                   #   and write to disk

        if h2o_data is not None:                                                                        # - If needed,
            lcm_files['h2o'].append(os.path.join(lcm_dir, "{}.h2o".format(lcm_slice_prefix)))           #   add H2O file to list,
            _write_lcmodel_raw(lcm_files['h2o'][-1], h2o_data, conj=conj)                               #   and write to disk

        for row_idx in range(1, spatial_dims[0] + 1):                                                   # Then for each voxel in the slice...

            for col_idx in range(1, spatial_dims[1] + 1):
            
                lcm_voxel_prefix = "{}_{}-{}".format(lcm_slice_prefix, row_idx, col_idx)                # - Create common prefix

                lcm_files['ctrl'].append(os.path.join(lcm_dir, "{}.control".format(lcm_voxel_prefix)))  # - Set control file
                lcm_files['csv'].append(os.path.join(lcm_dir, "{}.csv".format(lcm_voxel_prefix)))       # - Set output files
                lcm_files['ps'].append(os.path.join(lcm_dir, "{}.ps".format(lcm_voxel_prefix)))
                lcm_files['coo'].append(os.path.join(lcm_dir, "{}.coord".format(lcm_voxel_prefix)))

                _write_lcmodel_control(lcm_files, mrs_data, spatial_dims, slice_idx, row_idx, col_idx, h2o_data=h2o_data, ppm_range=(ppm_min, ppm_max))

    # Account for edge case where there's only one trace being analyzed

    if len(lcm_files['ps']) == 1:
        lcm_files['ps'][0] = lcm_files['ps'][0].replace('_sl1_1-1', '')
    if len(lcm_files['coo']) == 1:
        lcm_files['coo'][0] = lcm_files['coo'][0].replace('_sl1_1-1', '')

    return lcm_files

def _write_lcmodel_raw(raw_file, data, conj=False):

    """
    Write an MRSData object to a RAW file for LCModel.
    Writes both the MRS and H2O data
    Adapted from suspect.lcmodel.save_raw(). Changes are indicated with ***

    Parameters
    ----------
    raw_file : str
        The name of the RAW file to save
    data : MRSData object
        Contains metadata to write
    conj : bool
        Indicates whether data must be complex conjugated by LCModel (i.e., Bruker, Philips)
    """

    with open(raw_file, 'w') as fout:
        
        # Write sequence parameters

        fout.write(" $SEQPAR\n")
        fout.write(" ECHOT = {}\n".format(data.te))
        fout.write(" HZPPPM = {}\n".format(data.f0))
        fout.write(" SEQ = 'PRESS'\n")
        fout.write(" $END\n")

        # Write namelist

        fout.write(" $NMID\n")
        # if conj:                              # *** This line is TYPICALLY required to be T for Philips scans per LCModel manual (data must be complex conjugated)
        #     fout.write(" BRUKER = T\n")       # *** HOWEVER, we are accounting for this in the Philips DICOM loader, so we comment it out but leave it for posterity.
        fout.write(" FMTDAT = '(2E15.6)'\n")
        fout.write(" TRAMP = 1.\n")             # *** This line added to explicitly provide default value
        if data.transform is None:              # *** This block rearranged:
            vol = 1                             # *** - If we don't know the volume, use LCModel default of 1 mL.
            _print("WARNING: Saving LCModel data without a transform, using default voxel volume of 1 ml.")
        else:                                   # *** - Otherwise convert the volume from mm^3 to cc
            vol = data.voxel_volume() * 1e-3
        fout.write(" VOLUME = {}\n".format(vol))
        fout.write(" $END\n")

        # Write data

        for point in np.nditer(data, order='C'):
            fout.write("  {0: 4.6e}  {1: 4.6e}\n".format(float(point.real), float(point.imag)))

def _write_lcmodel_control(lcm_files, mrs_data, spatial_dims, slice_idx, row_idx, col_idx, h2o_data=None, ppm_range=(0.2, 4.0)):

    """
    Write an MRSData object and associated parameters to a control file for LCModel. 
    Adapted from suspect.lcmodel.write_all_files(). Changes are indicated with ***

    Parameters
    ----------
    lcm_files : dictionary
        Contains input/output file names
    mrs_data : MRSData object
        Spectra to analyze
    spatial_dims : NumPy array
        Contains the size of the three spatial dimensions
    slice_idx : int
        Indicates the slice to be analyzed
    row_idx : int
        Indicates the row to be analyzed
    col_idx : int
        Indicates the column to be analyzed
    h2o_data : MRSData object
        Water reference
    ppm_range : Tuple
        Contains the min and max ppm to analyze
    """

    with open(lcm_files['ctrl'][-1], 'wt') as fout:
        
        # Start

        fout.write(" $LCMODL\n")

        # Enable LCModel to run with correct OWNER/KEY pair

        fout.write(" OWNER = ''\n")
        fout.write(" KEY = 210387309\n")

        # Title of session

        fout.write(" TITLE = '{}'\n".format(mrs_data.metadata))

        # Input file names

        fout.write(" FILBAS = '{}'\n".format(lcm_files['basis']))           # Basis file input
        fout.write(" FILRAW = '{}'\n".format(lcm_files['raw'][-1]))         # Raw file input
        if h2o_data is not None:
            fout.write(" FILH2O = '{}'\n".format(lcm_files['h2o'][-1]))     # H2O file input

        # Output file names
        # Note: LCModel automatically appends sl[s]_[r]-[c] to PS and COO files, so drop them here

        # PS

        fout.write(" FILPS = '{}'\n".format('{}.ps'.format(lcm_files['ps'][-1].split('_sl')[0])))
        fout.write(" PGNORM = 'US'\n")                                      # US page size
        fout.write(" IPAGE2 = 0\n")                                         # Don't print the second page

        # CSV

        fout.write(" FILCSV = '{}'\n".format(lcm_files['csv'][-1]))
        fout.write(" LCSV = 11\n")                                          # 0 = no csv output, 11 = yes

        # COO

        fout.write(" FILCOO = '{}'\n".format('{}.coord'.format(lcm_files['coo'][-1].split('_sl')[0])))
        fout.write(" LCOORD = 9\n")                                         # 0 = no coordinates output, 9 = yes

        # Instrument and Acquisition

        fout.write(" HZPPPM = {}\n".format(mrs_data.f0))                    # Hz per ppm
        fout.write(" NUNFIL = {}\n".format(mrs_data.np))                    # Number of data points in one scan
        fout.write(" DELTAT = {}\n".format(mrs_data.dt))                    # Time (s) between two consecutive points
        fout.write(" ECHOT = {}\n".format(mrs_data.te))                     # Echo time

        # Data Shape
        # *** Swapped rows for columns here. Rows should be first dimension

        fout.write(" NDROWS = {}\n".format(spatial_dims[0]))                # Number of rows
        fout.write(" IROWST = {}\n".format(row_idx))                        # First row to analyze
        fout.write(" IROWEN = {}\n".format(row_idx))                        # Last row to analyze

        fout.write(" NDCOLS = {}\n".format(spatial_dims[1]))                # Number of columns
        fout.write(" ICOLST = {}\n".format(col_idx))                        # First column to analyze
        fout.write(" ICOLEN = {}\n".format(col_idx))                        # Last column to analyze

        fout.write(" NDSLIC = {}\n".format(spatial_dims[2]))                # Number of slices
        fout.write(" ISLICE = {}\n".format(slice_idx))                      # Slice to analyze

        # Analysis Window

        fout.write(" PPMST = {}\n".format(ppm_range[1]))                    # UPPER limit (left edge of decreasing ppm axis)
        fout.write(" PPMEND = {}\n".format(ppm_range[0]))                   # LOWER limit (right edge of decreasing ppm axis)

        # Finish up

        fout.write(" $END\n")

def _reim_nii(data):

    aff = data.transform
    re_img = np.real(data)
    im_img = np.imag(data)
    re_nii = _nii(re_img, aff) # account for transposed row/cols for MRS data
    im_nii = _nii(im_img, aff)
    return re_nii, im_nii

def _nii(img, aff):

    """
    Wrapper to make nibabel NIFTI with transposed rows/cols for Philips, Siemens, and LCModel data

    Parameters
    ----------
    img : numpy array
        image data
    aff : numpy array
        affine matrix
    """

    transpose_idxs = [i for i in range(len(img.shape))] # 0, 1, 2, ...
    transpose_idxs[0] = 1                               # 1, 1, 2, ...
    transpose_idxs[1] = 0                               # 1, 0, 2, ...

    # Return nibabel image where rows and cols are flipped
    #
    # Note: The functions load_dicom_philips and load_dicom_siemens both output 
    #       data transposed in the 0th and 1st dimensions compared to NIFTI 
    #       standards. Both need to be transposed before saving as NIFTIs. 
    #       This is done to be consistent with the way LCModel indexes. The 
    #       LCModel images ALSO need to be transposed prior to saving as NIFTIs. 
    #       Thus, this one function handles all three cases: Philips DICOM,
    #       Siemens DICOM, and LCModel output.
    #
    #       As such, to convert between FSLEyes indices and LCModel indices, 
    #       account for the different 0/1-based indexing and transposition with
    #       (x, y, z) in FSLEyes => (y+1, x+1, z+1) in LCModel

    return nib.Nifti1Image(np.transpose(img, transpose_idxs), affine=aff) 

def _write_data_nii(nii_dir, nii_prefix, time_img, aff, ppm):

    # Save Complex Time Domain Data

    time_real_nii = _nii(np.real(time_img), aff)
    time_real_nii_file = os.path.join(nii_dir, '{}_signal_real.nii.gz'.format(nii_prefix))
    nib.save(time_real_nii, time_real_nii_file)

    time_imag_nii = _nii(np.imag(time_img), aff)
    time_imag_nii_file = os.path.join(nii_dir, '{}_signal_imag.nii.gz'.format(nii_prefix))
    nib.save(time_imag_nii, time_imag_nii_file)

    time_abs_nii = _nii(np.abs(time_img), aff)
    time_abs_nii_file = os.path.join(nii_dir, '{}_signal_abs.nii.gz'.format(nii_prefix))
    nib.save(time_abs_nii, time_abs_nii_file)

    time_ang_nii = _nii(np.angle(time_img), aff)
    time_ang_nii_file = os.path.join(nii_dir, '{}_signal_ang.nii.gz'.format(nii_prefix))
    nib.save(time_ang_nii, time_ang_nii_file)

    time_abs_sum_nii = _nii(np.sum(np.abs(time_img), axis=3), aff)
    time_abs_sum_nii_file = os.path.join(nii_dir, '{}_signal_abs_sum.nii.gz'.format(nii_prefix))
    nib.save(time_abs_sum_nii, time_abs_sum_nii_file)

    # Save Complex Frequency Domain Data to Match Descending PPM

    spec_img = time_img.spectrum()

    spec_real_nii = _nii(np.real(spec_img), aff)
    spec_real_nii_file = os.path.join(nii_dir, '{}_spectra_real.nii.gz'.format(nii_prefix))
    nib.save(spec_real_nii, spec_real_nii_file)

    spec_imag_nii = _nii(np.imag(spec_img), aff)
    spec_imag_nii_file = os.path.join(nii_dir, '{}_spectra_imag.nii.gz'.format(nii_prefix))
    nib.save(spec_imag_nii, spec_imag_nii_file)

    spec_abs_nii = _nii(np.abs(spec_img), aff)
    spec_abs_nii_file = os.path.join(nii_dir, '{}_spectra_abs.nii.gz'.format(nii_prefix))
    nib.save(spec_abs_nii, spec_abs_nii_file)

    spec_ang_nii = _nii(np.angle(spec_img), aff)
    spec_ang_nii_file = os.path.join(nii_dir, '{}_spectra_ang.nii.gz'.format(nii_prefix))
    nib.save(spec_ang_nii, spec_ang_nii_file)

    spec_abs_sum_nii = _nii(np.sum(np.abs(spec_img), axis=3), aff)
    spec_abs_sum_nii_file = os.path.join(nii_dir, '{}_spectra_abs_sum.nii.gz'.format(nii_prefix))
    nib.save(spec_abs_sum_nii, spec_abs_sum_nii_file)

    # Save PPM (descending)

    assert ppm[0] > ppm[-1], 'PPM must be saved in descending format to match plotting convention. Aborting.'
    ppm_file = os.path.join(nii_dir, '{}.ppm'.format(nii_prefix))
    np.savetxt(ppm_file, ppm)

def _lcm2idxs(lcm_file, ext):

    # Extract voxel location from each LCM (csv or coord) file as row, col, slice, changing from 1 (LCModel) to 0 (python) indexed

    return np.array([int(i)-1 for i in lcm_file.split('.{}'.format(ext))[0].split('_sl')[-1].replace('_', '-').split('-')])[[1, 2, 0]]

def _nii2ref(in_nii, ref_nii, out_dir, interp='cubic'): # "regrid"

    # Make temporary directory
    
    tmp_dir = _mktmpdir(out_dir)

    # Regrid with MRTrix in temporary directory

    nii_file = os.path.join(tmp_dir, 'img.nii.gz')
    nib.save(in_nii, nii_file)
    ref_file = os.path.join(tmp_dir, 'ref.nii.gz')
    nib.save(ref_nii, ref_file)
    nii_regrid_file = os.path.join(tmp_dir, 'img_regrid.nii.gz')
    quiet_str = '-quiet' if not VERBOSE else ''
    grid_cmd = 'mrgrid {} regrid -template {} -strides {} {} -interp {} -fill nan -force {}'.format(nii_file, ref_file, ref_file, nii_regrid_file, interp, quiet_str)
    _run(grid_cmd)

    # Load image into memory before clearing from disk

    out_nii = nib.load(nii_regrid_file)
    out_img = out_nii.get_fdata()
    out_nii = nib.Nifti1Image(out_img, out_nii.affine)

    # Clear temporary directory

    _rmtmpdir(tmp_dir)

    return out_nii

def _mrs2idx(data):

    assert np.isin(len(data.shape), [4, 5]), 'Parameter data is expected to be 4- or 5-dimensional. Aborting.'

    mrs_img = data[:, :, :, 0] if len(data.shape) == 4 else data[0, :, :, :, 0]
    mrs_aff = data.transform
    idx_img = np.reshape(np.linspace(1, mrs_img.size, mrs_img.size), mrs_img.shape)
    idx_nii = _nii(idx_img, mrs_aff)
    return idx_nii

def _idx2mrs(idx_nii, data):

    idx_img = idx_nii.get_fdata()

    assert len(idx_img.shape) == 3 and np.isin(len(data.shape), [4, 5]), 'Parameter idx_nii is expected to be 3-dimensional, and data is expected to be 4- or 5-dimensional. Aborting.'
    if len(data.shape) == 4:    
        assert idx_img.shape[0] == data.shape[1] and idx_img.shape[1] == data.shape[0], 'Parameters idx_nii and data are expected to be relatively transposed in the first two spatial dimensions. Aborting.'
    else:
        assert idx_img.shape[0] == data.shape[2] and idx_img.shape[1] == data.shape[1], 'Parameters idx_nii and data are expected to be relatively transposed in the first two spatial dimensions. Aborting.'
    
    idx_img = np.transpose(idx_img, (1, 0, 2))
    if len(data.shape) == 5:
        idx_img = np.expand_dims(idx_img, axis=0)
    idx_img = np.expand_dims(idx_img, axis=len(idx_img.shape))
    idx_img = np.tile(idx_img, np.divide(data.shape, idx_img.shape).astype(int))

    return idx_img

def _mktmpdir(parent_dir=None):

    tmp_dir = os.path.join(parent_dir, 'tmp_{}'.format(datetime.now().strftime("%m%d%Y_%H%M%S")))
    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)
    return tmp_dir

def _rmtmpdir(tmp_dir):

    rm_cmd = 'rm -r {}'.format(tmp_dir)
    _run(rm_cmd)

def _run(cmd):

    """
    Wrapper to run command in command line

    Parameters
    ----------
    cmd : str
        Command to run
    """

    _print('RUNNING: {}'.format(cmd))
    subprocess.check_call(cmd, shell=True, executable='/bin/bash', stdout=None if VERBOSE else subprocess.DEVNULL)

def _coo2np(coo_files):

    # Get set up, considering one file edge case

    if len(coo_files) == 1:
        dims = (1, 1, 1)
    else:
        dims = _lcm2idxs(coo_files[-1], 'coord') + 1    # Extract spatial dimensions ("shape", index + 1) from last COO file (r, c, s)
    
    ppm_key = 'points on ppm-axis'                      # Define key lines needed to parse out data
    data_key = 'phased data points follow'
    fit_key = 'points of the fit to the data follow'
    bg_key = 'background values follow'
    diag_key = 'in following diagnostic table'

    ppm_start_idx = None                                # Extract indices and PPM, searching all possible COO files (should only be a problem if ALL voxels fail)
    data_start_idx = None
    fit_start_idx = None
    bg_start_idx = None
    diag_start_idx = None
    for coo_file in coo_files:
        with open(coo_file, 'r') as coo_fobj:           
            coo_lines = coo_fobj.read().splitlines()
        for i, coo_line in enumerate(coo_lines):
            if ppm_key in coo_line:
                ppm_start_idx = i+1
            elif data_key in coo_line:
                data_start_idx = i+1
            elif fit_key in coo_line:
                fit_start_idx = i+1
            elif bg_key in coo_line:
                bg_start_idx = i+1
            elif diag_key in coo_line:
                diag_start_idx = i+1
        if ppm_start_idx is not None and data_start_idx is not None and fit_start_idx is not None and bg_start_idx is not None and diag_start_idx is not None: # If all are found, stop searching
            break
        else:
            ppm_start_idx = None
            data_start_idx = None
            fit_start_idx = None
            bg_start_idx = None
            diag_start_idx = None
    
     # Build NaN 4D arrays to hold row, col, slice, spec data if proper data found

    if ppm_start_idx is not None and data_start_idx is not None and fit_start_idx is not None and bg_start_idx is not None and diag_start_idx is not None:
        _parse_coo_pts = lambda lines : [float(s) for s in ' '.join(lines).split()]
        ppm = _parse_coo_pts(coo_lines[ppm_start_idx:data_start_idx-1])
    else:
        ppm = [np.nan]

    data_img = np.empty((dims[0], dims[1], dims[2], len(ppm))) 
    data_img[:] = np.nan

    fit_img = np.empty((dims[0], dims[1], dims[2], len(ppm))) 
    fit_img[:] = np.nan

    bg_img = np.empty((dims[0], dims[1], dims[2], len(ppm))) 
    bg_img[:] = np.nan

    # Loop through all COO files (considering edge cases) and assign values to array

    if ppm_start_idx is not None and data_start_idx is not None and fit_start_idx is not None and bg_start_idx is not None and diag_start_idx is not None:
        for coo_file in coo_files:
            with open(coo_file, 'r') as coo_fobj:
                coo_lines = coo_fobj.read().splitlines()
            if len(coo_lines) > diag_start_idx:
                idxs = _lcm2idxs(coo_file, 'coord') if len(coo_files) > 1 else (0, 0, 0)
                data_img[idxs[0], idxs[1], idxs[2], :] = _parse_coo_pts(coo_lines[data_start_idx:fit_start_idx-1])
                fit_img[idxs[0], idxs[1], idxs[2], :] = _parse_coo_pts(coo_lines[fit_start_idx:bg_start_idx-1])
                bg_img[idxs[0], idxs[1], idxs[2], :] = _parse_coo_pts(coo_lines[bg_start_idx:diag_start_idx-1])

    # Clean outputs

    ppm = np.expand_dims(np.array(ppm), axis=0) # Add dimension to ppm so it prints horizontally
    if len(coo_files) == 1:
        data_img = data_img[0, 0, :] # Keep third dimension so they print horizontally
        fit_img = fit_img[0, 0, :]
        bg_img = bg_img[0, 0, :]

    return ppm, data_img, fit_img, bg_img

def _csv2line(csv_files, file_name, lbl_name, dummy=False):

    # Read CSV file

    with open(csv_files[0], 'r') as csv_fobj:
        csv_raw = csv_fobj.read().splitlines()
    
    # If empty, return header info with empty data

    if np.size(csv_raw) == 0:
        csv_lbls_split = ['File', 'Label(s)']
        csv_line_split = [file_name, lbl_name]
        return csv_lbls_split, csv_line_split

    # Get labels, replacing the first and second column

    csv_lbls = csv_raw[0]
    csv_lbls_split = csv_lbls.split(', ')
    csv_lbls_split[0] = 'File'
    csv_lbls_split[1] = 'Label(s)'
    
    # Get Data

    csv_line = csv_raw[1]
    csv_line_split = csv_line.split(',')
    csv_line_split[0] = file_name
    csv_line_split[1] = lbl_name
    for i in range(2, len(csv_line_split)):
        csv_line_split[i] = float(csv_line_split[i]) if not dummy else np.nan

    # Clear CSV file

    rm_cmd = 'rm {}'.format(' '.join(csv_files))
    _run(rm_cmd)

    return csv_lbls_split, csv_line_split

def _lines2csv(csv_lbls, csv_lines, ratio_ref, out_dir, out_prefix):

    # Find inputs with non-missing data and fill in missing lines

    for csv_lbl in csv_lbls: # What happens if they're all missing??
        if np.size(csv_lbl) > 2:
            csv_lbls = csv_lbl
            break
    for i in range(len(csv_lines)):
        if np.size(csv_lines[i]) == 2:
            csv_nan = [np.nan for n in range(len(csv_lbl)-2)]
            csv_lines[i] = csv_lines[i] + csv_nan

    # Convert lines to dictionary to compute desired ratios

    csv_dict = {}
    for k in range(len(csv_lines)):
        for i, csv_lbl in enumerate(csv_lbls):
            if csv_lbl in csv_dict:
                csv_dict[csv_lbl].append(csv_lines[k][i])
            else:
                csv_dict[csv_lbl] = [csv_lines[k][i]]

    # Compute desired ratios

    csv_keys = list(csv_dict.keys()) # Iterative keys must be kept constant as we are inherently changing keys
    den_key = ratio_ref
    assert den_key in csv_keys, 'Ratio reference was not found in LCModel output. Aborting.'
    for csv_key in csv_keys:
        if '/' in csv_key:
            csv_dict.pop(csv_key) # Remove old ratio
            num_key = csv_key.split('/')[0]
            num_list = csv_dict[num_key]
            den_list = csv_dict[ratio_ref]
            assert len(num_list) == len(den_list), 'A different number of elements in the numerator and denominator concentrations was detected. Aborting.'
            ratio_list = [np.float64(num_list[i]) / den_list[i] for i in range(len(num_list))]
            ratio_key = '{}/{}'.format(num_key, den_key)
            csv_dict[ratio_key] = ratio_list

    # Convert back to CSV

    csv_hdrs = list(csv_dict.keys())
    csv_data = []
    for k in range(len(csv_lines)):
        csv_data.append([])
        for csv_hdr in csv_hdrs:
            csv_data[k].append(str(csv_dict[csv_hdr][k]))
        csv_data[k] = ','.join(csv_data[k])
    csv_hdrs = ','.join(csv_hdrs)
    csv_text = '\n'.join([csv_hdrs] + csv_data)
    
    # Save File

    csv_file = os.path.join(out_dir, '{}_{}_ratios_by_label.csv'.format(out_prefix, ratio_ref))
    with open(csv_file, 'w') as csv_fobj:
        csv_fobj.write(csv_text)
    return csv_file    

def _check_nii_hdrs(nii1, nii2):

    img1 = nii1.get_fdata()
    aff1 = nii1.affine

    img2 = nii2.get_fdata()
    aff2 = nii2.affine

    return np.allclose(aff1, aff2) and np.allclose(img1.shape[:3], img2.shape[:3])

def _write_complex_signal(data, data_file):

    assert np.isin(len(data.shape), [4, 5]), 'MRSData is expected to be 4- or 5-dimensional. Aborting.'
    signal = data
    if len(data.shape) == 5:
        signal = data[0]
    assert np.prod(signal.shape) == signal.shape[-1], 'Only one complex signal can be written to text. Aborting.'
    np.savetxt(data_file, signal.squeeze(), fmt='%.6e%+.6ej')

def _dummy_signal(data):

    return np.cos(np.linspace(-2*np.pi, 2*np.pi, data.shape[-1]))

def _print(msg):

    if VERBOSE:
        print('masimrs.py: {}'.format(msg))
        
###############
# Main Script #
###############

if __name__ == '__main__':

    ##############
    # Get Inputs #
    ##############

    parser = ap.ArgumentParser(description='2D PRESS MRSI Processing for Philips Enhanced DICOM and Siemens DICOM data with LCModel')
    parser.add_argument('dcm_file', metavar='/input/file.dcm', help='path to the input DICOM file')
    parser.add_argument('out_prefix', metavar='/output/prefix', help='path (with prefix) of the output files')
    parser.add_argument('scanner_type', metavar='scanner_type', help='string indicating scanner type ("philips" or "siemens")')
    parser.add_argument('-x', '--novox', action='store_true', help='omit voxel-wise analysis')
    parser.add_argument('-s', '--seg', metavar=('/seg.nii.gz', '1+2,3'), action='append', nargs=2, help='path to segmentation file and label(s) indicating regions to analyze together (+) or separately (,) (can be used multiple times, default = do NOT perform regional analysis)')
    parser.add_argument('-r', '--ref', metavar='REF', default='Cr+PCr', help='string indicating which LCModel metabolite to use when computing ratios (default = Cr+PCr)')
    parser.add_argument('-n', '--nthreads', metavar='N', default=1, help='positive integer indicating number of threads to use for voxel-wise fitting (default = 1)')
    parser.add_argument('-v', '--verbose', action='store_true', help='print progress to console')
    args = parser.parse_args()

    ################
    # Parse Inputs #
    ################

    _print('PARSING INPUTS...')

    dcm_file = args.dcm_file
    if not os.path.exists(dcm_file):
        raise MRSException('Input DICOM {} does not exist. Aborting.'.format(dcm_file))
    
    out_dir = os.path.dirname(args.out_prefix)
    if out_dir == '':
        out_dir = '.'
    if not os.path.exists(out_dir):
        raise MRSException('Output directory {} does not exist. Aborting.'.format(out_dir))
    out_prefix = os.path.basename(args.out_prefix)
    
    scanner = args.scanner_type
    if scanner != 'philips' and scanner != 'siemens':
        raise MRSException('Scanner must be either "philips" or "siemens". Aborting.')
    if scanner == 'philips':
        load_dicom = load_dicom_philips
        nifti_image = nifti_image_philips
        save_nifti = save_nifti_philips
        save_lcmodel = save_lcmodel_philips
    elif scanner == 'siemens':
        load_dicom = load_dicom_siemens
        nifti_image = nifti_image_siemens
        save_nifti = save_nifti_siemens
        save_lcmodel = save_lcmodel_siemens

    num_threads = int(args.nthreads)
    if num_threads < 1:
        raise MRSException('Number of threads must be a positive integer. Aborting.')

    if args.seg is not None:
        seg_arg = args.seg
        seg_files = [s[0] for s in seg_arg]
        seg_lbls = [s[1] for s in seg_arg]
        for i, seg_file in enumerate(seg_files):
            if not os.path.exists(seg_file):
                raise MRSException('Segmentation file {} does not exist. Aborting.'.format(seg_file))
            if not all([s in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '0', '+', ',', '-'] for s in seg_lbls[i]]):
                raise MRSException('Segmentation labels must be integers and separated by either , or +. Aborting.')
            str_lbls = seg_lbls[i].split(',')
            int_lbls = []
            for str_lbl in str_lbls:
                int_lbls.append([int(s) for s in str_lbl.split('+')])
            seg_lbls[i] = int_lbls
    else:
        seg_files = []
        seg_lbls = []

    ratio_ref = args.ref

    no_voxel = args.novox

    VERBOSE = args.verbose

    #############
    # Load Data #
    #############

    _print('*****************************')
    _print('* CONVERTING DICOM TO NIFTI *')
    _print('*****************************')

    _print('LOADING DICOM...')                              # Load DICOM
    data = load_dicom(dcm_file)
    
    _print('SAVING AS NIFTI...')                            # Save data as NIFTI
    save_nifti(out_dir, out_prefix, data)

    ###########################
    # Run Voxel-wise Analysis #
    ###########################

    if not no_voxel:

        _print('*******************************')
        _print('* RUNNING VOXEL-WISE ANALYSIS *')
        _print('*******************************')

        _print('SAVING LCMODEL FILES...')                       # Save LCModel Files
        lcm_files = save_lcmodel(out_dir, out_prefix, data)

        _print('RUNNING LCMODEL...')                            # Run LCModel
        with Pool(num_threads) as p:
            p.map(run_lcmodel, lcm_files['ctrl'])

        _print('REFORMATTING LCMODEL FILES...')                         # Clean up LCModel files
        lcm_files['pdf'] = ps2pdf(lcm_files['ps'])                      # Clean up PS files
        lcm_files['nii'] = csv2nii(lcm_files['csv'], data.transform)    # Clean up CSV files
        coo2nii(lcm_files['coo'], data.transform)                       # Clean up COORD files
        txt_file = ctrl2txt(lcm_files['ctrl'])                          # Clean up CONTROL files

    else:

        _print('********************************')
        _print('* SKIPPING VOXEL-WISE ANALYSIS *')
        _print('********************************')

    ############################
    # Run Region-wise Analysis #
    ############################   

    if len(seg_files) > 0 and len(seg_lbls) > 0:

        _print('********************************')
        _print('* RUNNING REGION-WISE ANALYSIS *')
        _print('********************************')

        # Generate index image for MRS slab

        _print('IDENTIFYING VOXEL INDICES...')

        idx_nii = _mrs2idx(data)
        idx_img = idx_nii.get_fdata()
        idx_aff = idx_nii.affine

        # Prep for one output CSV file

        csv_lbls = []
        csv_lines = []

        # For each segmentation file...
        
        for i, seg_file in enumerate(seg_files):
            
            seg_file_name = os.path.basename(seg_file).replace('.nii.gz', '')

            # Regrid index image to segmentation file

            _print('REGRIDDING VOXEL INDICES TO {}...'.format(seg_file))

            seg_nii = nib.load(seg_file)
            idx_regrid_nii = _nii2ref(idx_nii, seg_nii, out_dir, interp='nearest') # Any time converting between MRS and NII, need to account for transposition
            idx_regrid_img = idx_regrid_nii.get_fdata()

            # Then for each label in the segmentation file

            for seg_lbl in seg_lbls[i]:

                # Prep output names

                seg_lbl_name = '+'.join([str(l) for l in seg_lbl[:20]])
                if len(seg_lbl) > 20:
                    seg_lbl_name = '{}...'.format(seg_lbl_name)

                seg_prefix = '{}_{}_{}'.format(out_prefix, seg_file_name, seg_lbl_name) 
                
                # Identify how much each MRS voxel overlaps with the segmentation label

                _print('IDENTIFYING MRS OVERLAP WITH LABEL(S) {}...'.format(seg_lbl_name))

                seg_img = np.isin(seg_nii.get_fdata(), seg_lbl)
                nib.save(nib.Nifti1Image(seg_img.astype(float), seg_nii.affine), os.path.join(out_dir, '{}_mask.nii.gz'.format(seg_prefix)))
                seg_idxs = idx_regrid_img[seg_img] 
                seg_idxs = seg_idxs[np.logical_not(np.isnan(seg_idxs))].astype(int)
                seg_idxs_unique = np.unique(seg_idxs)
                idx_count_img = np.zeros_like(idx_img)
                for seg_idx in seg_idxs_unique:
                    idx_count_img[idx_img == seg_idx] = np.count_nonzero(seg_idxs == seg_idx) # zero if no overlap found

                # Begin pooling analysis, convert MRS voxel overlap count to percentage of total overlapping voxels

                _print('RUNNING SIGNAL POOLING ANALYSIS...')

                idx_weight_img = idx_count_img / np.sum(idx_count_img) # All NaNs if no overlap found due to 0 divided by 0
                idx_weight_nii = nib.Nifti1Image(idx_weight_img, idx_aff) 
                overlap_found = np.sum(idx_weight_img) > 0
                
                # If overlap found, normalize weights so they sum to 1 for weighted averaging
                # and computed weighted average signal based on overlap with segmentation label

                if overlap_found:

                    _print('OVERLAP FOUND, POOLING TIME DOMAIN SIGNALS WITH WEIGHTED AVERAGE...')
                    data_weights = _idx2mrs(idx_weight_nii, data) # Any time converting between MRS and NII, need to account for transposition
                    data_weighted = data * data_weights
                    data_weighted_avg = np.mean(data_weighted, axis=(-2, -3, -4), keepdims=True)
                    data_weighted_avg.metadata = 'Weighted average signal in label(s) {} from {} in {}'.format(seg_lbl_name, seg_file_name, data_weighted_avg.metadata)

                # If no overlapping idxs are found, generate nonsense data so that LCModel doesn't break
                # and make sure this is obvious in the PDF

                else:

                    _print('WARNING: NO OVERLAP FOUND, GENERATING DUMMY DATA TO ENSURE CODE COMPLETION...')
                    data_weighted = data - data + _dummy_signal(data)
                    data_weighted_avg = np.mean(data_weighted, axis=(-2, -3, -4), keepdims=True)
                    data_weighted_avg.metadata = 'THIS IS A DUMMY SIGNAL! No overlap found in label(s) {} from {} in {}'.format(seg_lbl_name, seg_file_name, data_weighted_avg.metadata)

                # Save Intermediates, considering if no overlap was found

                _print('SAVING INTERMEDIATES FOR QA...')
                nib.save(idx_weight_nii, os.path.join(out_dir, '{}_weights.nii.gz'.format(seg_prefix)))
                data_weighted_avg_file = os.path.join(out_dir, '{}.signal'.format(seg_prefix))
                _write_complex_signal(data_weighted_avg if overlap_found else np.nan*np.zeros_like(data_weighted_avg), data_weighted_avg_file)
                data_weighted_spec_file = os.path.join(out_dir, '{}.spectra'.format(seg_prefix))
                _write_complex_signal(data_weighted_avg.spectrum() if overlap_found else np.nan*np.zeros_like(data_weighted_avg), data_weighted_spec_file)
                data_weighted_ppm_file = os.path.join(out_dir, '{}.ppm'.format(seg_prefix))
                np.savetxt(data_weighted_ppm_file, np.expand_dims(data_weighted_avg.frequency_axis_ppm(), axis=0))

                # Run LCModel for the weighted average, with consideration that there's only one spectra

                _print('RUNNING LCMODEL...')
                lcm_files = save_lcmodel(out_dir, seg_prefix, data_weighted_avg)                # Save LCModel files
                run_lcmodel(lcm_files['ctrl'][0])                                               # Run LCModel

                # Clean LCModel Outputs, with consideration that there's only one spectra or if no overlap was found

                _print('REFORMATTING LCMODEL FILES...')
                lcm_files['pdf'] = ps2pdf(lcm_files['ps'])                                                              # Convert PS to PDF                      
                coo2txt(lcm_files['coo'], dummy=not overlap_found)                                                      # Convert COO to TXT
                txt_file = ctrl2txt(lcm_files['ctrl'])                                                                  # Convert CTRL to TXT
                csv_lbl, csv_line = _csv2line(lcm_files['csv'], seg_file_name, seg_lbl_name, dummy=not overlap_found)   # Save and clean CSV contents
                csv_lbls.append(csv_lbl)
                csv_lines.append(csv_line)

        # Merge regional CSV contents

        _print('MERGING RATIOS FOR ALL REGIONS INTO ONE CSV FILE...')
        csv_file = _lines2csv(csv_lbls, csv_lines, ratio_ref, out_dir, out_prefix)
        
    else:

        _print('*********************************')
        _print('* SKIPPING REGION-WISE ANALYSIS *')
        _print('*********************************')

    #############
    # Finish Up #
    #############

    _print('*********')
    _print('* DONE! *')
    _print('*********')
