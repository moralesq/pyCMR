# Manuel A. Morales (mmorale5@bidmc.harvard.edu)
# Postdoctoral Research Fellow
# Beth Israel Deaconess Medical Center (BIDMC)
# Nezafat-Med Cardiovascular

# Created: 2021_07_28

import os
import glob
import pydicom
import pandas as pd

def make_dataset(dir, max_dataset_size=float("inf")):
  """ Create a pandas DataFrame containing with header information from all dicoms withint a `dir`. 
      Header information includes the dicom `filename` and `PatientName` to facility querying and loading. 

  Returns
  -------
  headers : pandas DataFrame with dicom header information sorted by dicoms filenames. Sorting is important
            to identify slice locations in cases where `SeriesInstanceUID` is the same for all dicom files.  
  """
  assert os.path.isdir(dir), '%s is not a valid directory' % dir

  headers = pd.DataFrame(_init_header())
  for root, _, fnames in os.walk(dir):
      for fname in fnames:

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
            'SliceThickness':[]}

  return header