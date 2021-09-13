# Manuel A. Morales, PhD (mmorale5@bidmc.harvard.edu)
# CMR Enthusiast

import os
import glob
import pydicom
import warnings
import pandas as pd

def export_dicom_headers(dicom_dir, savedir, subject_dirs_target=None, max_dataset_size=float("inf"), encode_name=False):
  """Read selected header information from dicoms within a `dicom_dir` folder.
  
     The following first-level subfolders are expected:
     dicom_dir/
     |-----subject_dir_1
     |
     |-----subject_dir_2
     |
     |
     |
     |-----subject_dir_N

     All dicoms within each `subject_dir_n` subfolder will be searched, independent of folder structure. 
     If only a few subject_dir_n folders are needed, set subject_dirs_target=['subject_dir_2', 'subject_dir_5'].
     To test this script, use `max_dataset_size=m` to load only a small amount of `m` files for each subject. 

  Exports
  -------
  dicom headers:  If encode_name=True, exported individually for each subject with structure:

                  savedir/
                  |-----YYYY_MM_DD_FNLNS_dicoms_headers.csv
                  |
                  
                  where YYYY = year, MM = month, DD = day refer to the acquisition date; FN=first name, LN=last mame, and S=sex.

                  If encode_name=False, the entire name without sex will be used, i.e., YYYY_MM_DD_SubjectName_dicoms_headers.csv

  """
  os.makedirs(savedir, exist_ok=True)

  subject_dirs = sorted(glob.glob(os.path.join(dicom_dir, '*')))

  for subject_dir in subject_dirs:
    
    subject_folder = os.path.basename(subject_dir)
    if subject_dirs_target is not None:
      if not any([subject_folder==subject_dir_n for subject_dir_n in subject_dirs_target]): continue
    print('Reading:', subject_folder)

   # try:
    dicoms_headers_folder = make_dataset(dir=subject_dir, max_dataset_size=max_dataset_size)

    # ideally there should be a single subject within each folder, but this is not always the case. One could 
    # reject the folder altogether all try to remove the incorrectly named subject (e.g., `Retro Recon`)
    for PatientName in dicoms_headers_folder.PatientName.unique():
      dicoms_headers = dicoms_headers_folder[dicoms_headers_folder.PatientName==PatientName]

      assert len(dicoms_headers.PatientName.unique()) == 1
      assert len(dicoms_headers.PatientSex.unique()) == 1
      assert len(dicoms_headers.StudyDate.unique()) == 1

      if encode_name:
        PatientSex = dicoms_headers.PatientSex.unique().item()
        PatientLastName, PatientFirstName = dicoms_headers.PatientName.unique().item().split(' ')[:2]
        PatientEncode = PatientFirstName[:2]+PatientLastName[:2]+PatientSex
      else:
        PatientEncode = dicoms_headers.PatientName.unique().item()

      StudyDate = dicoms_headers.StudyDate.unique().item()
      year, month, day = str(int(StudyDate[:4])), str(int(StudyDate[4:6])), str(int(StudyDate[6:]))

      # PatientID should have format YYYY_MM_DD_FNLNS
      PatientID = '_'.join([year, month, day, PatientEncode]).upper()

      dicoms_headers['PatientID'] = PatientID
      dicoms_headers.to_csv(os.path.join(savedir, PatientID + '_' + 'dicoms_headers.csv'))
      print('Exported dicoms_headers.csv of shape', dicoms_headers.shape)

      #except:
      #  warnings.warn('============== Could not read dicoms in folder %s =============== '%(subject_folder))

      
def make_dataset(dir, max_dataset_size=float("inf"), 
                 ext_exclude=['.json','.nii.gz', '.bvec', '.bval', '.DS_Store']):
  """ Create a pandas DataFrame containing header information from all dicoms withint a `dir`. 
      Header information includes the dicom `filename` and `PatientName` to facility querying and loading. 

  Input
  -----
  dir : path of folder to search. 

  max_dataset_size: int, the maximum amount of files to search. This is useful when testing the code.  

  ext_exclude : list of strings indicating file extensions. This is useful to skip non-dicom files present in 
                the folder. The default includes common extensions found along dicom files. 

  Returns
  -------
  headers : pandas DataFrame with dicom header information sorted by dicoms filenames. Sorting is important
            to identify slice locations in cases where `SeriesInstanceUID` is the same for all dicom files.  
  """
  assert os.path.isdir(dir), '%s is not a valid directory' % dir

  headers = pd.DataFrame(_init_header())
  for root, _, fnames in os.walk(dir):
      for fname in fnames:
          if fname.endswith(tuple(ext_exclude)): continue

          is_valid, header = is_valid_dicom(os.path.join(root, fname))

          if is_valid:
              if isinstance(header, str):
                print('Dicom is not valid:', header, os.path.join(root, fname))
              else:
                headers = pd.concat((headers, header))
          else:
              print('Could not open/read file:', os.path.join(root, fname))
          
          if len(headers) == max_dataset_size: 
              return headers
              
  return headers.sort_values(by='FileName')

def is_valid_dicom(filename):
  """ Verifies:
       (1) a dicom can be loaded using pydicom.read_file(filename)
       (2) the dicom contains primary data (i.e., it rejects dicoms processsed by Siemens in-line) 
  
  Returns
  -------
  is_dicom : bool, True or False

  dicom    : if valid dicom returns dicom header, otherwise returns reason for not valid. 
  """
  try:
    dicom = pydicom.read_file(filename)     
    # ADDITIONAL TECHNICAL CHECKS
    # remove siemens segmentation outputs. 
    if 'InlineVF' in dicom.SeriesDescription: return True, 'InlineVF'
    # check images come primary  data (i.e., MR or CT scanners)
    if 'Secondary' in dicom.file_meta[0x0002, 0x0002].repval: return True, 'Secondary'

    
  except:
      return False, None

  return True, read_header(dicom, filename)

def read_header(dicom, filename):
  """ Reads selected dicom keys as specified in _init_header() into a pandas DataFrame. 
      The `filename` of the dicom is included in the DataFrame for downstream quality control.
  """

  header = _init_header()

  header['FileName']    += [filename]
  header['PatientName'] += [str(dicom[0x0010, 0x0010].value).replace('^', ' ')]
  for key in header.keys():
    if key not in ['FileName', 'PatientName']:
      try:
        header[key] += [dicom[key].value]
      except:
        header[key] += [None]

  return pd.DataFrame(header)

def _init_header():
  """ Returns empty dictionary of selected dicom keys (e.g., PatientSex, StudyDate, SeriesDescription). 
      The `PatientName` and dicom `FileName` are also included in the header dictionary. 
  """
  header = {'FileName':[],
            'PatientName':[], 
            'PatientSex':[], 
            'PatientAge':[], 
            'StudyDate':[], 
            'SeriesTime':[], 
            'AcquisitionTime':[], 
            'SeriesInstanceUID':[], 
            'StudyInstanceUID':[],
            'SOPInstanceUID':[],
            'SeriesDescription':[],
            'ProtocolName':[], 
            'SeriesTime':[], 
            'TriggerTime':[], 
            'InstanceNumber':[], 
            'ImageOrientationPatient':[],  
            'ImagePositionPatient':[],
            'SliceLocation':[], 
            'PixelSpacing':[], 
            'SliceThickness':[], 
            'SpacingBetweenSlices':[]}

  return header