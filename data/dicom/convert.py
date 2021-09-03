# Manuel A. Morales, PhD (mmorale5@bidmc.harvard.edu)
# CMR Enthusiast

import os
import glob
import pydicom
import warnings
import numpy as np
import pandas as pd
import nibabel as nib

from . import folder 

def export_niftis(headers_dir, 
                   savedir,
                   SeriesDescriptionContainsTrues=[], 
                   SeriesDescriptionContainsFalse=[],
                   select_subject=None):

    dicoms_header_subjects = glob.glob(os.path.join(headers_dir, '*'))              

    for dicoms_header_subject in dicoms_header_subjects:
        
        # PatientID should have format YYYY_MM_DD_PatientName or YYYY_MM_DD_FNLNS
        PatientID = os.path.basename(dicoms_header_subject).strip('_dicoms_headers.csv')

        if select_subject is not None:
            if PatientID not in select_subject: continue
        
        subject_dicoms_header_complete = pd.read_csv(dicoms_header_subject, index_col=0)
        
        print('='*50)
        print(subject_dicoms_header_complete.PatientName.unique())
        print('='*50)

        # it is better to treat each desired Series description separately 
        for SeriesDescriptionContainsTrue in SeriesDescriptionContainsTrues:
            subject_dicoms_header = subject_dicoms_header_complete.copy()

            
            # only include dicoms matching the `SeriesDescriptionContainsTrue` string. From those, exclude the ones
            # matching strings in the  `SeriesDescriptionContainsFalse`. 
            SeriesDescriptions_unique = subject_dicoms_header.SeriesDescription.unique()
            if (len(SeriesDescriptions_unique)==1)&(SeriesDescriptionContainsTrue in SeriesDescriptions_unique):
                # I had to add this modification because for some reason the mask was failing when there was
                # only a single SeriesDescription.
                SeriesDescription = subject_dicoms_header.SeriesDescription.unique()[0]
                Targets = (subject_dicoms_header.SeriesDescription==SeriesDescription).to_numpy()
            else:
                Targets = subject_dicoms_header.SeriesDescription.str.contains(SeriesDescriptionContainsTrue, regex=True).to_numpy()
            
            if SeriesDescriptionContainsFalse:
                Exclude = ~subject_dicoms_header.SeriesDescription.str.contains('|'.join(SeriesDescriptionContainsFalse)).to_numpy()
                
                subject_dicoms_header = subject_dicoms_header.iloc[Targets&Exclude,:]
            else:
                subject_dicoms_header = subject_dicoms_header.iloc[Targets,:]

            os.makedirs(os.path.join(savedir, PatientID, 'niftis'), exist_ok=True)
            os.makedirs(os.path.join(savedir, PatientID, 'slice_locations'), exist_ok=True)
            os.makedirs(os.path.join(savedir, PatientID, 'nifti2dicom_paths'), exist_ok=True)
            for SeriesDescription in subject_dicoms_header.SeriesDescription.unique():
                index = '_'.join(SeriesDescription.split(' '))
                index = index.replace('/', '_')
                save_name = PatientID + '_' + index
                
                series_dicoms_header = subject_dicoms_header[subject_dicoms_header.SeriesDescription==SeriesDescription]

                for SeriesInstanceUID in sorted(series_dicoms_header['SeriesInstanceUID'].unique()):

                    instance_dicoms_header = series_dicoms_header[series_dicoms_header.SeriesInstanceUID==SeriesInstanceUID]   

                    print('instance shape', instance_dicoms_header.shape)
                    AcquisitionTime = min(instance_dicoms_header.AcquisitionTime)

                    save_npy_dic = os.path.join(savedir, PatientID, 'nifti2dicom_paths', save_name+'_dicom_paths_AcqTime_%d.npy'%(AcquisitionTime))
                    save_npy = os.path.join(savedir, PatientID, 'slice_locations', save_name+'_slice_locations_AcqTime_%d.npy'%(AcquisitionTime))
                    save_nii = os.path.join(savedir, PatientID, 'niftis', save_name+'_AcqTime_%d.nii.gz'%(AcquisitionTime))
                    
                    #print('Saving slices locations to:', save_npy)
                    print('Saving nifti to:', save_nii)

                    # read and save locations (in mm) of image slices in nifti. This is useful to align with other datasets. 
                    SliceLocations = instance_dicoms_header[instance_dicoms_header.SeriesInstanceUID==SeriesInstanceUID].SliceLocation.unique()     
                    np.save(save_npy, SliceLocations)
                    
                    # convert dicoms to nifti and save. 
                    sax_nifti, dicom_4D_paths = read_cine_protocol(series_dicom_header=instance_dicoms_header)
                    # TODO: handle QA errors when creating nifti
                    if sax_nifti is None: continue
                    sax_nifti.to_filename(save_nii)
                    np.save(save_npy_dic, dicom_4D_paths)

                    print('nifti shape:', sax_nifti.shape)
                    print('slices shape:', SliceLocations.shape)

def read_cine_protocol(series_dicom_header):
    """" Read a cine protocol and convert to 4D NIFTI format. This function can fail if basic assumptions about 
         the cine acquisition are violated. 

    Input
    -----
    series_dicom_header: a pandas DataFrame containing the necessary dicom header information to construct 
                        the 4D array. 

    Return
    ------
    sax_4D : concatenated dicom images into a 4D nifti array containing affine information. 

    dicom_4D_paths : dictionary of dicom paths containing sorted by slice and time. This is useful to convert
                     the niftis back to dicom (e.g., to generate segmentations!)
    """
    assert len(series_dicom_header.StudyInstanceUID.unique()) == 1, 'Trying to read dicoms from multiple studies!'

    SliceLocations = series_dicom_header.SliceLocation.unique()
    number_of_slices = len(SliceLocations) 

    phases_per_slice = [len(series_dicom_header[series_dicom_header.SliceLocation==SliceLocation].InstanceNumber) 
                        for SliceLocation in series_dicom_header.SliceLocation.unique()]
    number_of_phases = phases_per_slice[0]

    if len(np.unique(phases_per_slice)) != 1:
        warnings.warn('Number of phases is variable across slice locations! Could be real or error, check!.')
        return None
   
    print('Found cine study with (number_of_slices, number_of_phases)', number_of_slices, number_of_phases)
    pixel_array = pydicom.read_file(series_dicom_header.iloc[0].FileName).pixel_array
   
    sax_4D = np.zeros((pixel_array.shape +(number_of_slices, number_of_phases)), dtype=pixel_array.dtype)
    
    dicom_4D_paths = {}
    for SliceIndex, SliceLocation in enumerate(sorted(SliceLocations)):
        slice_header = series_dicom_header[series_dicom_header.SliceLocation==SliceLocation]
        dicom_4D_paths[SliceIndex] = []
        for InstanceIndex, InstanceNumber in enumerate(sorted(slice_header.InstanceNumber)):
            DicomFileName = slice_header[slice_header.InstanceNumber==InstanceNumber].FileName.item()
            dicom = pydicom.read_file(DicomFileName)
            sax_4D[:,:,SliceIndex,InstanceIndex] += dicom.pixel_array

            dicom_4D_paths[SliceIndex] += [DicomFileName]

    affine = read_affine(series_dicom_header.iloc[series_dicom_header.SliceLocation.argmin()])

    return nib.Nifti1Image(sax_4D, affine=affine), dicom_4D_paths

def _string_to_list_of_floats(x): return list(np.array(x.strip("'[].").split(','), dtype=float))

def read_affine(df):
    """ Read the 4x4 affine matrix from a pandas DataFrame `df` containing dicom header information. 
    """
    SliceThickness          = [df.SliceThickness]
    PixelSpacing            = _string_to_list_of_floats(df.PixelSpacing)
    ImageOrientationPatient = _string_to_list_of_floats(df.ImageOrientationPatient)
    ImagePositionPatient    = _string_to_list_of_floats(df.ImagePositionPatient)

    Zooms                   = np.array(PixelSpacing+SliceThickness, dtype=float)
    ImageOrientationPatient = np.array(ImageOrientationPatient, dtype=float)
    ImagePositionPatient    = np.array(ImagePositionPatient, dtype=float)
    
    ijk2ras = extract_cosines(ImageOrientationPatient)

    ijk2ras = (ijk2ras*np.array([-1,-1,1])).T
    ImagePositionPatient = ImagePositionPatient*np.array([-1,-1,1])

    affine  = np.stack((ijk2ras[:,0]*Zooms[0],
                        ijk2ras[:,1]*Zooms[1],
                        ijk2ras[:,2]*Zooms[2],
                        ImagePositionPatient), axis=1)

    return np.vstack((affine,[[0,0,0,1]])) 

def extract_cosines(ImageOrientationPatient):
    """ Extract the cos(theta) values describing the orientation of the acquisition plane. 
    """
    row_cosine    = np.array(ImageOrientationPatient[:3])
    column_cosine = np.array(ImageOrientationPatient[3:])
    slice_cosine  = np.cross(row_cosine, column_cosine)
    return np.stack((row_cosine, column_cosine, slice_cosine))


