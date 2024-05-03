# Reproducing the results / rerunning the code

The outputs of some scripts are liable to change due to f.ex. changed content in databases from when they were last run.
Files for the intermediate steps are therefore included.

## Environment setup
For some steps, external dependency versions can be relevant. We therefore include the `environment` folder that contains lists over all installed dependencies.
Python 3.10 was installed inside a conda environment, while dependencies inside the environment were installed via pip.

## Package installation
Once the general environment is set up, install the packages inside `code_packages`. This can be done from the packages' folder using `pip install .` .
Of the two packages provided, `cofactors` provides the relevant modeling code, while `plotting` only contains helper functions for visualization.


## External data

### Extracting proteomics averages
The proteomics data from mitoPARP cell lines and parental HEK293 used are included under `external_data`,
the raw data contained in `proteomics.xlsx`, the processed data integrated into the models are in `293parp_abundance_ratios_mapped.csv`.
`uniprot_acc_to_ensemble_gene.tsv` contains mappings to convert between the identifiers in the data and the model. 
They were downloaded from taken from [david.ncifcrf.gov](https://david.ncifcrf.gov/).
All the steps performed can be followed in `code_modeling/extract_proteomics.ipynb`

### Extracting Km
Kms can be downloaded to a default location using the script `./code_modeling/download_kms.py`.
This script relies on having the cofactors package installed.
Note that the database contents may change over time as more data is added.

### Mitocore model
`external_data/mitocore` contains the Mitocore model used, along with metadata and data extracted from the models.
See also the [official source](https://www.mrc-mbu.cam.ac.uk/research-resources-and-facilities/mitocore).


## Model building
Building the models happens both inside Matlab (proteomics integration) and Python (NAD integration).
The process therefore involves multiple scripts.
There are also models created for the supplemental data.
Relevant code is in `code_modeling`.

### Proteomics integration via GIMME
`code_modeling/integration_gimme.m` creates the 293mitoPARP and parental HEK293 models that are used in the main analyses.
`code_modeling/paramscan_gimme.m` creates the necessary series of models to estimate the effects of different objective fractions on the GIMME models.

### Integration of NAD levels
`code_modeling/integration_nad.ipynb` creates the models parameterized for NAD levels. 
Here, we also solve the models with FBA and pFBA and save the results to be used later.
Running this script relies on the `cofactors` package.

## Model analyses and visualizations
Some of the analyses (solving the models with pFBA and FBA) is already run along with the model creation.
The scripts here contain the concrete analyses performed with the solutions along with the visualizations.
The relevant scripts are contained in `code_analysis`.

### Main analyses
The core analyses are performed in `analyses_parp.ipynb`, where the illustrations related to the 293mitoPARP / HEK293 models are also created.
Running the analyses depends on having the `plotting` package installed.


### Supplemental analyses
`check_kms.ipynb` gives an overview over the Kms extracted from Brenda and SabioRK.
`comparison_brenda_sabiork.ipynb` shows how using the Brenda or SabioRK databases for obtaining Kms influences the models.
`gene_mapping.ipynb` analyzes how many of the IDs in the proteomics dataset could be mapped to the model and vice versa. 
`graphing_gimme_mapping.ipynb` gives an overview of how many of the reactions in Mitocore were ultimately covered by the proteomics integration.
`parameter_scan_c0.ipynb` tracks model results across a range of different C_0 parameters.
`scan_params_gimme.ipynb` tracks model results across a range of GIMME parameters.
`scan_decision_functions.ipynb` compares the results of different decision functions when choosing between multiple Km values.
`scan_fva_frac.ipynb` tracks model results across a range of FVA objective fractions.

# Result files
The `results` folder contains the direct analysis results from the models created.
In the folder itself the data are sorted by metric, as they are used for the analyses in the main text.

`results/supplements` contains similar data, sorted by model instead to give a better overview of the calculations.