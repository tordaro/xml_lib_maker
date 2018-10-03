from xml_libray_maker import *
import numpy as np
import pandas as pd
from os import mkdir, listdir
import xml.etree.ElementTree as et
from scipy.constants import pi, g

material_df = [pd.read_excel(lib_path, sheet_name=sheet_name, dtype=float)]

@pytest.mark.parametrize("material_df", material_df)
def test_col_names(material_df):
    '''Checks that the header is as expected.'''
    header = ['Diameter [mm]',
              'Massetetthet [kg/m]',
              'MBL [tonn]',
              'MBL [kN]',
              'E-modul [Pa]',
              'Materialkoeffisient'] # Expected header
    for col_name in header:
        assert col_name in material_df.columns, \
        '{} not in dataset.'.format(col_name)

@pytest.mark.parametrize("material_df", material_df)
def test_col_types(material_df):
    '''Checks that all columns are float64.'''
    col_types = material_df.dtypes
    for i, col_type in enumerate(col_types):
        assert col_type == np.float64, \
        '{} is not float64.'.format(material_df.columns[i])

@pytest.mark.parametrize("material_df", material_df)
def test_increasing_observations(material_df):
    '''Checks that all, except last two columns, are 
    increasing in value with each observation.'''
    rows, cols = material_df.shape
    for j in range(cols-2):
        for i in range(rows-1):
            assert material_df.iloc[i, j] < material_df.iloc[i+1, j], \
            '{} is not increasing with each observation.'.format(material_df.columns[j])

@pytest.mark.parametrize("material_df", material_df)
def test_interdepence(material_df):
    '''Checks that MBL [tonn] and MBL [kN]
    corresponds within tolerance of 1000 N.'''
    tolerance = 1.0 # kN
    low_bound = (material_df['MBL [kN]'] - tolerance
                 <= material_df['MBL [tonn]'] * g)
    up_bound = (material_df['MBL [kN]'] + tolerance\
                >= material_df['MBL [tonn]'] * g)
    bounds = ~(low_bound & up_bound)
    assert low_bound.all() and up_bound.all(), \
    'Discrepancy between MBL [tonn] and MBL [kN] in \n{}'.format(material_df[bounds])