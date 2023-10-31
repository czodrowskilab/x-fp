#import the essentials

import pytest
import os
from xgboost import XGBClassifier

# define a fixture that reads the DataFrame
@pytest.fixture(scope="module")
def example_data():
    try:
        df = pd.read_csv('./example_data/df_test_herg.tsv', sep='\t')
    except Exception as e:
        pytest.fail(f"Failed to load DataFrame: {e}")
    return df

# define a fixture that loads the XGBoost classifier model
@pytest.fixture(scope="module")
def model():
    model_fp2 = XGBClassifier()
    fp2_model_path = './example_data/herg_model-xgb.json'
    try:
        model_fp2.load_model(fp2_model_path)
    except Exception as e:
        pytest.fail(f"Failed to load model: {e}")
    return model_fp2

def test_imports():
    '''
    Test if all the required modules can be imported.
    '''
    try:
        from collections import defaultdict
        from multiprocessing import Pool, cpu_count
        import numpy as np
        import pandas as pd
        import seaborn as sns
        from matplotlib import pyplot as plt
        from rdkit.Chem import Draw, BondType, AllChem as Chem, rdFingerprintGenerator
        from tqdm import tqdm
    except ImportError:
        assert False, "One or more required modules could not be imported."
    else:
        assert True, "All required modules were imported successfully."

#importing to test the functions 
import pandas as pd
from rdkit.Chem import Draw, BondType, AllChem as Chem, rdFingerprintGenerator
import numpy as np

def test_example_data_folder():
    '''
    Test if the 'example_data' folder exists and contains the expected files.
    '''
    # specify the path to the 'example_data' folder
    example_data_folder = './example_data'

    # check if the folder exists
    assert os.path.isdir(example_data_folder), f"'{example_data_folder}' does not exist."

    # specify the expected files
    expected_files = ['df_test_herg.tsv', 'herg_model-xgb.json']

    # check if each expected file is in the folder
    for filename in expected_files:
        file_path = os.path.join(example_data_folder, filename)
        assert os.path.isfile(file_path), f"'{file_path}' does not exist."

def test_example_data():
    '''
    Test if the example data file can be loaded as a DataFrame and has the expected shape and columns.
    '''
    #check if the file is readable
    df = pd.read_csv('./example_data/df_test_herg.tsv', sep='\t')

    # check if the file has the expected number of rows and columns
    assert df.shape == (1216, 7), f"Expected (1216, 7) but got {df.shape}."

    # check if the file has the expected columns
    expected_columns = ['molecule_chembl_id', 'IC50', 'units', 'pIC50', 'InChIKey', 'Smiles', 'herg_class']
    for column in expected_columns:
        assert column in df.columns, f"'{column}' is not in the file."

def test_imports():
    '''
    Test if all the required modules can be imported.
    '''
    # list of modules to test
    modules = [
        'xfp',
        'xfp.core_xfp',
        'xfp.report_maker',
        'xfp.xfp_tree_explainer',
        'xfp.xfp_xgb',
    ]

    # try to import each module
    for module in modules:
        try:
            __import__(module)
        except ImportError:
            pytest.fail(f"Failed to import '{module}'")

#import X-FP modules now
from xfp.core_xfp import FingerprintManager

def test_fetch_fp_from_mol():
    '''
    Test if the method can generate fingerprint from a molecule object without any errors.
    '''
    # create a test molecule object - caffeine :) 
    test_smiles = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    test_mol = Chem.MolFromSmiles(test_smiles)

    # create a FingerprintManager object
    fp_manager = FingerprintManager()

    # test if the method can generate fingerprint without any errors
    fp, bit_info = fp_manager.fetch_fp_from_mol(test_mol)

    # test if the generated fingerprint has the expected shape and data type
    assert isinstance(fp, np.ndarray)
    assert fp.shape == (2048,)

    # test if the generated fingerprint contains the expected number of on bits
    assert fp.sum() == 25

    # test if the generated fingerprint contains the expected bit info
    assert isinstance(bit_info, dict)
    assert len(bit_info) == 25

def test_fetch_fp_from_smiles():
    '''
    Test if the method can generate fingerprint from a SMILES without any errors.
    '''
    # create a test molecule object - caffeine again :) 
    test_smiles = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"

    # create a FingerprintManager object
    fp_manager = FingerprintManager()

    # test if the method can generate fingerprint without any errors
    fp, bit_info = fp_manager.fetch_fp_from_smiles(test_smiles)

    # test if the generated fingerprint has the expected shape and data type
    assert isinstance(fp, np.ndarray)
    assert fp.shape == (2048,)

    # test if the generated fingerprint contains the expected number of on bits
    assert fp.sum() == 25

    # test if the generated fingerprint contains the expected bit info
    assert isinstance(bit_info, dict)
    assert len(bit_info) == 25

def test_fetch_fp_from_df(example_data):
    '''
    Test if the method can generate fingerprints from a DataFrame without any errors.
    '''
    # create a FingerprintManager object
    fp_manager = FingerprintManager()

    # test if the method can generate fingerprints without any errors
    fp_manager.fetch_fp_from_df(example_data, smiles_col = "Smiles", n_bits= 4096, use_chirality = True)

    # test if the generated fingerprints have the expected shape and data type
    assert fp_manager.fp.shape == (1216, 4096), "Fingerprint array has the wrong shape." 

# to-do  : test further