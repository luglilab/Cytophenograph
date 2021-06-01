# Cytophenograph

Identifies subpopulations in high-dimensional single-cell data. Cytophenograph is a computational pipeline that was developed to avoid the disadvantages of manual gating. The pipeline is developed using Python3, the clustering method adopted is a custom version of Phenograph (https://github.com/jacoblevine/PhenoGraph) where we insert a blocked seed. Besides Phenograph pipeline needs the following package installed: Pandas,Numpy,Sklearn for data parsing and exploring and openTSNE,Seaborn,Matplotlib for data visualization. This method is adaptative both in terms of dimensionality and sample size, making it suitable in a range of settings for which single-cell population structure is of interest, including other cancers or healthy tissues, and for use with other emerging single-cell technologies.

## 1) Installation 

Install Miniconda
Miniconda is a Python distribution, package manager, and virtual environment solution. We recommend installing Miniconda with Python 3 (miniconda3), as many bioinformatics packages are now transitioning to Python 3. You can still install Python 2 software with miniconda3 by passing the python=2.7 flag when you create a new environment; otherwise the default Python version will be Python 3.

Begin by downloading Miniconda and following the associated installation instructions.

https://docs.conda.io/en/latest/miniconda.html

### Create your cytophenograph environment and install the dependences

Test if miniconda3 is installed


```python
which conda
```

Clone our repository


```python
git clone https://github.com/luglilab/Cytophenograph
```

Create a new environment

for installation on linux machine execute this command:
```python
conda env create -n cytophenograph2 -f ./Cytophenograph/environment_cytophenograph2_linux.yml
conda activate cytophenograph2
```
for installation on mac machine execute this command:
```python
conda env create -n cytophenograph2 -f ./Cytophenograph/environment_cytophenograph2_mac.yml
conda activate cytophenograph2
```

Install Phenograph


```python
pip install -e ./Cytophenograph/Phenograph_LugliLab --use-feature=2020-resolver
pip install scipy==1.4.1 --use-feature=2020-resolver
```

Move on Phenograph folder



```python
python ./Cytophenograph/cytophenograph.v2_0.py --help
```


Test Execution 
```python
abs_path=$(pwd)
mkdir $abs_path/Cytophenograph/output_test
# Run Phenograph
python ./Cytophenograph/cytophenograph.v2_0.py -i $abs_path/Cytophenograph/Test_dataset/CD8_Panel_II_channelvalues_GA_downSampled/ -o $abs_path/Cytophenograph/output_test -k 300 -m $abs_path/Cytophenograph/Test_dataset/CD8_bulk_markers_to_exclude.txt -n Test -t 10 -p $abs_path/Cytophenograph/Test_dataset/Info_file_bulk_Test.xlsx -c Phenograph
# Run PARC
python ./Cytophenograph/cytophenograph.v2_0.py -i $abs_path/Cytophenograph/Test_dataset/CD8_Panel_II_channelvalues_GA_downSampled/ -o $abs_path/Cytophenograph/output_test -k 300 -m $abs_path/Cytophenograph/Test_dataset/CD8_bulk_markers_to_exclude.txt -n Test -t 10 -p $abs_path/Cytophenograph/Test_dataset/Info_file_bulk_Test.xlsx -c Parc
# Run Phenograph and Parc
python ./Cytophenograph/cytophenograph.v2_0.py -i $abs_path/Cytophenograph/Test_dataset/CD8_Panel_II_channelvalues_GA_downSampled/ -o $abs_path/Cytophenograph/output_test -k 300 -m $abs_path/Cytophenograph/Test_dataset/CD8_bulk_markers_to_exclude.txt -n Test -t 10 -p $abs_path/Cytophenograph/Test_dataset/Info_file_bulk_Test.xlsx -c Both
```
# 

Pipeline has been testen on Linux and Mac OS. 
Know bug:  Scipy  version must <1.4.1, During the execution of "pip install scipy==1.4.1 --use-feature=2020-resolver". User could obtain this warning "ERROR: scanorama 1.6 requires intervaltree==2.1.0, but you'll have intervaltree 3.0.2 which is incompatible.
anndata 0.7.4 requires pandas>=1.0, but you'll have pandas 0.25.3 which is incompatible.
phenograph 1.5.7 requires scipy>=1.5.1, but you'll have scipy 1.4.1 which is incompatible."
Please ignore this warning. 

**File preparation:**

 - [ Input folder]: This folder must contains **only** csv file exported from flowjo after compensation and trasformation. CSV must have same identical header. 
 - [Marker list]: This file must contains the features (channel or marker) to exclude. Please check that channel is included in the header of all csv files. (Example here: https://github.com/luglilab/Cytophenograph/blob/master/Test_dataset/CD8_bulk_markers_to_exclude.txt )
 - [Pheno File]: Excel file (Example here: https://github.com/luglilab/Cytophenograph/blob/master/Test_dataset/Info_file_bulk_Test.xlsx ) with the following column "Sample Cell_type EXP ID Time_point Condition Count". Number of row should be the same of input CSV. 


**Output**
 - [Output Folder]: Empty folder where user will find .h5ad file ( ready to use for Cellxgene https://chanzuckerberg.github.io/cellxgene/), FCScluster[Phenograph or Parc] folder and FCSsample[Phenograph or Parc] folder with Tot_counts.txt and Tot_percentage.txtand with absolute and percentage frequency and log.txt with analysis execution information. 
 
 **Graphics output**
 h5ad file with UMAP and others graphical output could be open with Cellxgene ( https://chanzuckerberg.github.io/cellxgene/ ). 

### Please cite:
Alvisi G, Brummelman J, Puccio S, Mazza EM, Tomada EP, Losurdo A, Zanon V, Peano C, Colombo FS, Scarpa A, Alloisio M, Vasanthakumar A, Roychoudhuri R, Kallikourdis M, Pagani M, Lopci E, Novellis P, Blume J, Kallies A, Veronesi G, Lugli E. IRF4 instructs effector Treg differentiation and immune suppression in human cancer. J Clin Invest. 2020 Jun 1;130(6):3137-3150. doi: 10.1172/JCI130426. PMID: 32125291; PMCID: PMC7260038.
