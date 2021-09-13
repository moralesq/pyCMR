import os
import glob
import numpy as np
import pandas as pd
import nibabel as nib

def create_single_slice_dicoms_dataframe(nifti_dir):
    """ Creates a dataframe containing information of single-slice dicoms. This is information is needed 
        to accurately link dicoms that belong to a single acquisition but whose SeriesInstanceUID is different. 
    
    """

    subject_folders  = glob.glob(os.path.join(nifti_dir, '*'))   
    nifti_dir_concat = nifti_dir + '_concatenated_slices'

    group_df = {'SubjectID':[], 'FileID':[], 'NiftiPath':[], 'AcqTime':[], 'Shape':[], 'SliceLocation':[]}
    for subject_folder in subject_folders:

        nifti_paths = glob.glob(os.path.join(subject_folder, 'niftis', '*'))

        for nifti_path in nifti_paths:
            FileID, AcqTime = os.path.basename(nifti_path).strip('.nii.gz').split('_AcqTime_')
            shape  = nib.load(nifti_path).shape

            assert len(shape) == 4
            # code is intendent to treat slices not linked, therefore slice dim must be 1!
            if shape[2] == 1:
                slice_location = np.load(os.path.join(subject_folder, 
                                                    'slice_locations', 
                                                    FileID + '_slice_locations_AcqTime_%s.npy' %(AcqTime)))

                group_df['SubjectID']      += [os.path.basename(subject_folder)]
                group_df['FileID']         += [FileID]
                group_df['NiftiPath']      += [nifti_path]
                group_df['AcqTime']        += [AcqTime]
                group_df['Shape']          += [shape]
                group_df['SliceLocation']  += [slice_location[0]]

    group_df = pd.DataFrame(group_df)
    group_df.to_csv(nifti_dir_concat+'.csv')

    return group_df

def concatenate_dicoms(file_df):
    """ Concatenate multiple 3D niftis (2D slice + time) into a single 4D nifti
    """

    nifti_ref = nib.load(file_df.NiftiPath.iloc[0])

    nx, ny, _, n_frames = nifti_ref.shape

    n_slices = len(file_df)

    nifti_4D = np.zeros((nx, ny, n_slices, n_frames), dtype=nifti_ref.get_fdata().dtype)

    for slice_index in range(n_slices):
        image_slice = nib.load(file_df.NiftiPath.iloc[slice_index]).get_fdata()
        nifti_4D[:,:,slice_index,:] += image_slice.squeeze()

    nifti_4D = nib.Nifti1Image(nifti_4D, affine=nifti_ref.affine)

    return nifti_4D


def link_dicoms_with_different_UIDs(nifti_dir):

    df = create_single_slice_dicoms_dataframe(nifti_dir)
    nifti_dir = '/Users/laptop/Desktop/Maaike/MRI_Cines_2_niftis'
    nifti_dir_concat = nifti_dir + '_concatenated_slices'

    for SubjectID in df.SubjectID.unique():

        # dataframe for dicoms with the same SubjectID
        subject_df = df[df.SubjectID==SubjectID]

        for FileID in subject_df.FileID.unique():
            
            # dicoms with the same FileID and SubjectID
            file_df = subject_df[subject_df.FileID==FileID].sort_values(by=['AcqTime'])
            
            for Shape in file_df.Shape.unique():
                
                # dicoms with the same Shape, FileID, and SubjectID
                shape_df = file_df[file_df.Shape==Shape]

                if len(shape_df) == len(shape_df.SliceLocation.unique()):
                    # all FileIDs have the same shape with slices at different location, 
                    # therefore assume all belong to same acquisition
                    nifti_4D = concatenate_dicoms(shape_df)

                    if nifti_4D.shape[2] >= 5:
                        # save file
                        output_folder = os.path.join(nifti_dir_concat, os.path.basename(SubjectID))

                        os.makedirs(os.path.join(output_folder, 'nifti2dicom_paths'), exist_ok=True)
                        os.makedirs(os.path.join(output_folder, 'niftis'), exist_ok=True)
                        os.makedirs(os.path.join(output_folder, 'slice_locations'), exist_ok=True)

                        
                        nifti_4D.to_filename(os.path.join(output_folder, 
                                                        'niftis', 
                                                        os.path.basename(list(shape_df.NiftiPath)[0])))

                        print(nifti_4D.shape)
                else:
                    # at this point we assume this are likely from the same acquisition. 
                    # We will loop through them in order of acquisition time; if a slice 
                    # is repeated then assume is a different acquisition. 

                    slice_locations = []
                    acquisitions    = []
                    for index, slice_location in zip(shape_df.index, shape_df.SliceLocation):
                        
                        if slice_location not in slice_locations:
                            slice_locations += [slice_location]
                            acquisitions    += [shape_df[shape_df.index==index]]
                        else:

                            nifti_4D = concatenate_dicoms(pd.concat(acquisitions))

                            
                            if nifti_4D.shape[2] >= 5:
                                # save file
                                output_folder = os.path.join(nifti_dir_concat, os.path.basename(SubjectID))

                                os.makedirs(os.path.join(output_folder, 'nifti2dicom_paths'), exist_ok=True)
                                os.makedirs(os.path.join(output_folder, 'niftis'), exist_ok=True)
                                os.makedirs(os.path.join(output_folder, 'slice_locations'), exist_ok=True)

                                
                                nifti_4D.to_filename(os.path.join(output_folder, 
                                                                'niftis', 
                                                                os.path.basename(acquisitions[0].NiftiPath.item())))
    
                                print(nifti_4D.shape)
                            
                            slice_locations = []
                            acquisitions    = []