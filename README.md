# X-FP - eXplainable FingerPrints

![x-fp logo](src/xfp/resources/X-FP_logo.png)

## Description
X-FP (eXplainable FingerPrints) is a software package which explains the model outputs for the QSAR models trained exclusively on Morgan Fingerprints. It is model-agnostic. 

## Installation
The tool can be installed using either `mamba` (recommended) or `pip`-only. First of all clone the repository and 
switch to the cloned directory:
```bash
git clone https://github.com/czodrowskilab/x-fp
cd x-fp
```

### Using Mamba
A conda-based installation, e.g. Anaconda, Miniconda, Miniforge or Mambaforge, must be available on your system and
`mamba` must be installed in the base environment. If no existing conda installation is available, visit 
[Mambaforge](https://github.com/conda-forge/miniforge#mambaforge) for more information.

To create the conda environment, activate it and install X-FP into the environment, run the following commands:
```bash
mamba env create -f environment_generic.yml
mamba activate xfp-env
pip install .
```
There are also environment files for specific operating systems available, e.g. [environment_linux64.yml](environment_linux64.yml) for Linux
and [environment_osx.yml](environment_osx.yml) for macOS. These files contain a fixed package set that was used for developing and testing.

### Using only Pip
To install the tool using only `pip`, it is recommended to create a virtual environment first. For this you have to use
a system Python interpreter of version 3.10 or higher. You can check this by calling `python3 -V`. Run the following
commands to create a virtual environment, activate it and install X-FP into the environment:
```bash
python3 -m venv xfp-env
source xfp-env/bin/activate
pip3 install rdkit .
```
If there is an error during the first command, you probably have to install the `venv` module first. This can be done
for example by running `sudo apt install python3-venv` on Debian/Ubuntu. You can also skip the first two commands and
just call `pip3 install rdkit .` directly, but this will install the tool and all dependencies into your system Python
environment which is not recommended.

## Usage

X-FP is fairly straight-forward to use. 

### Test Data and Model
hERG bioactivity data was accessed from CHEMBL 27. A [XGBoost binary classification model](test_data/herg_model-xgb.json) was trained on 6078 compounds. 
[Test data](test_data/df_test_herg.tsv) comprising of 1216 compounds is included. SHAP TreeExplainer object for the model is also included as a [pickled object](test_data/pickled_shap_explainer_fp4.pkl).

### SHAP TreeExplainer
```python

# Import X-FP for TreeExplainer
from xfp_tree_explainer import XFPTreeExplainer

# Read the pickled shap_explainer for your model
with open("test_data/pickled_shap_explainer_fp4.pkl", 'rb') as file:
    shap_explainer = pickle.load(file)

# Make a X-FP tree explainer object
xfp = XFPTreeExplainer(shap_explainer)

# Load the data (here, the test set - in case you want to visualize the important bits and the chemical substructures they encode in the test set)
df_test = pd.read_csv('test_data/df_test_herg.tsv', sep= '\t')

# Generate the fingerprints of the radius and number of bits on which the model was trained
xfp.fetch_fp_from_df(df_test, smiles_col = "Smiles", n_bits= 4096, use_chirality = True)

# The dataframe (now with extra columns for the RDKit molecule object and the Morgan fingerprints) can be checked by calling:
xfp.input_df

# This dataframe is a normal Pandas dataframe and can be used as such
xfp.input_df.info()
xfp.input_df.head()

# Generate the fragments. This is usually the time-consuming step, and can depend on the size of the dataset and parallelization.
xfp.generate_frags(n_jobs = 10, verbose = True)

# To visualize the number of substructures encoded per bit across the dataset, a box plot can be generated:
xfp.make_substruct_per_bit_plot()

# The SHAP values should be generated first
xfp.generate_shap_values()

# The SHAP summary plot can be generated by calling:
xfp.generate_shap_summary_plot()

# The report can be generated for any desired bits:
xfp.generate_bit_analysis_report([807, 3839, 2698, 3291, 2855, 800], report_title = "X-FP SHAP TreeExplainer Report: hERG Test Set")

```

### XGBoost Feature Importance
```python

# Import X-FP for XGBoost in-built feature importance methods
from xfp_xgb import XFPXGB

# Read the FP2 XGBoost classification model
model_fp2 = XGBClassifier()
fp2_model_path = "test_data/herg_model-xgb.json"

model_fp2.load_model(fp2_model_path)

# Make a X-FP XGBoost feature importance object
xfp = XFPXGB(model_fp2)

# Load the data (here, the test set - in case you want to visualize the important bits and the chemical substructures they encode in the test set)
df_test = pd.read_csv('test_data/df_test_herg.tsv', sep= '\t')

# Generate the fingerprints of the radius and number of bits on which the model was trained
xfp.fetch_fp_from_df(df_test, smiles_col = "Smiles", n_bits= 4096, use_chirality = True)

# The dataframe (now with extra columns for the RDKit molecule object and the Morgan fingerprints) can be checked by calling:
xfp.input_df

# This dataframe is a normal Pandas dataframe and can be used as such
xfp.input_df.info()
xfp.input_df.head()

# Generate the fragments. This is usually the time-consuming step, and can depend on the size of the dataset and parallelization.
xfp.generate_frags(n_jobs = 10, verbose = True)

# To visualize the number of substructures encoded per bit across the dataset, a box plot can be generated:
xfp.make_substruct_per_bit_plot()

# The feature importance should be generated first. Any of the in-built feature importance methods can be used.
xfp.make_feature_importance(importance_type = "gain")

# The feature importance plot can be generated by calling:
xfp.generate_feature_importance_plot()

# The report can be generated for any desired bits:
xfp.generate_bit_analysis_report([807, 3839, 2698, 3291, 2855, 800], report_title = "X-FP XGB FI Report: hERG Test Set")

```

## Cite us
Reference for citation will be added soon.

## Authors
- Marcel Baltruschat ([GitHub](https://github.com/mrcblt))
- Aishvarya Tandon ([GitHub](https://github.com/aish-tan))

## Acknowledgements
- Hanna Rieger: Validation case study
- Adam Adamczyk: Logo design

## License
[MIT License](LICENSE)
