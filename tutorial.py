import os
from data.dicom import convert, folder, link_UIDs

EXPORT_HEADERS = True
CONVERT_TO_NIFTI = True
LINK_UIDs_DICOMS = True

for dicomfolder in ['MRI_Cines_2', 'MRI_Cines_1/TO-BE-DONE']:
        
    maindir     = '/Users/laptop/Dropbox (MIT)/MARTINOS_REPRODUCIBILITY/Pigs'
    dir         = os.path.join(maindir, dicomfolder)
    headers_dir = dir + '_dicoms_headers'
    nifti_dir   = dir + '_niftis'
        
    if EXPORT_HEADERS:
        folder.export_dicom_headers(dicom_dir=dir, 
                                    savedir=headers_dir, 
                                    subject_dirs_target=None,
                                    max_dataset_size=float("inf"))

    if CONVERT_TO_NIFTI:
        convert._export_niftis(headers_dir=headers_dir, 
                               savedir=nifti_dir,
                               SeriesDescriptionContainsTrues=['cine', 'tf2d'], 
                               SeriesDescriptionContainsFalse=['T1_map', 'long_axis', '4CH', '2CH','LongAx','Long'],
                               select_subject=None)

    if LINK_UIDs_DICOMS:

        link_UIDs.link_dicoms_with_different_UIDs(nifti_dir)

   