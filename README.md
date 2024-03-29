# Cytophenograph

Identifies subpopulations in high-dimensional single-cell data. Cytophenograph is a computational pipeline that was developed to avoid the disadvantages of manual gating. The pipeline is developed using Python3, the clustering method adopted is a custom version of Phenograph (https://github.com/jacoblevine/PhenoGraph) where we insert a blocked seed. Besides Phenograph pipeline needs the following package installed: Pandas,Numpy,Sklearn for data parsing and exploring and UMAP,Seaborn,Matplotlib for data visualization. This method is adaptative both in terms of dimensionality and sample size, making it suitable in a range of settings for which single-cell population structure is of interest, including other cancers or healthy tissues, and for use with other emerging single-cell technologies.

## 1) Installation 

Install Miniconda
Miniconda is a Python distribution, package manager, and virtual environment solution. We recommend installing Miniconda with Python 3 (miniconda3), as many bioinformatics packages are now transitioning to Python 3. 
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


### Installation on LINUX machine
Tested on Debian GNU/Linux server

**Strategy 1 : Use YML file to clone environment** 
```python
conda env create -n Cytophenograph6 -f ./Cytophenograph/environment_cytophenograph6_linux.yml
conda activate Cytophenograph6
pip install -e ./Cytophenograph/FlowSOM_LugliLab
pip install phenograph==1.5.7
```
**Strategy 2 : Execute the following command** 
```python
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --add channels bpteague
conda create --name Cytophenograph6 python=3.9  r-base hnswlib
conda activate Cytophenograph6
pip install pyVIA
conda install -c bioconda scanorama -y
conda install -c conda-forge datashader -y
conda install -c anaconda openpyxl -y
pip install -e ./Cytophenograph/FlowSOM_LugliLab
pip install phenograph==1.5.7
pip install fcsy
conda install -c bioconda bioconductor-flowai -y
conda install scikit-image -y
```

### Installation on MAC machine
Tested on computer with ios 10.15.7 
**Strategy 1 : Use YML file to clone environment** 
```python
conda env create -n Cytophenograph6 -f ./Cytophenograph/environment_cytophenograph6_mac.yml
conda activate Cytophenograph6
pip install -e ./Cytophenograph/FlowSOM_LugliLab
pip install phenograph==1.5.7
```
**Strategy 2 : Execute the following command** 
```python
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --add channels bpteague
conda create --name Cytophenograph6 python=3.9  r-base hnswlib 
conda activate Cytophenograph6
pip install pyVIA
conda install -c bioconda scanorama -y
conda install -c conda-forge datashader -y
conda install -c anaconda openpyxl
pip install -e ./Cytophenograph/FlowSOM_LugliLab
pip install phenograph==1.5.7
pip install fcsy
conda install -c bioconda bioconductor-flowai
conda install scikit-image
```

### Docker 
```
Please visit Docker Hub:
https://hub.docker.com/r/sinnamone/cytophenograph6

or 

use Dockerfile

docker build  -t cytophenograph6 .

To test the execution:

abs_path=$(pwd)

mkdir -p $abs_path/Cytophenograph/output_test

docker run --entrypoint /bin/bash -v $abs_path:/data -w /data -p 8891:8891  -it cytophenograph6

Inside the docker execute:

python cytophenograph.v6.py -i /data/Example_One_Inputs/CSVFiles/ -o /data/output_test/ -k 60 -m /data/Example_One_Inputs/InfoFile/marker.txt -n Test -t 10 -p /data/Example_One_Inputs/InfoFile/Info_file_bulk_Test.xlsx -c Phenograph

```

###  Move on Phenograph folder



```python
python ./Cytophenograph/cytophenograph.v6.py --help
```


###  Test Execution 
```python
abs_path=$(pwd)
mkdir -p $abs_path/Cytophenograph/output_test
# Run Phenograph
python ./Cytophenograph/cytophenograph.v6.py -i $abs_path/Cytophenograph/Example_One_Inputs/CSVFiles/ -o $abs_path/Cytophenograph/output_test/ -k 60 -m $abs_path/Cytophenograph/Example_One_Inputs/InfoFile/marker.txt -n Test -t 10 -p $abs_path/Cytophenograph/Example_One_Inputs/InfoFile/Info_file_bulk_Test.xlsx -c Phenograph
# Run VIA
python ./Cytophenograph/cytophenograph.v6.py -i $abs_path/Cytophenograph/Example_One_Inputs/CSVFiles/ -o $abs_path/Cytophenograph/output_test/ -w 30 -z 1.0 -m $abs_path/Cytophenograph/Example_One_Inputs/InfoFile/marker.txt -n Test -t 10 -p $abs_path/Cytophenograph/Example_One_Inputs/InfoFile/Info_file_bulk_Test.xlsx -c VIA
# Run Flowsom
python ./Cytophenograph/cytophenograph.v6.py -i $abs_path/Cytophenograph/Example_One_Inputs/CSVFiles/ -o $abs_path/Cytophenograph/output_test/ -x 5 -y 31 -m $abs_path/Cytophenograph/Example_One_Inputs/InfoFile/marker.txt  -n Test -t 10 -p $abs_path/Cytophenograph/Example_One_Inputs/InfoFile/Info_file_bulk_Test.xlsx -c FlowSOM

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
 - [Output Folder]: Empty folder where user will find .h5ad file ( ready to use for Cellxgene https://chanzuckerberg.github.io/cellxgene/), CSVcluster folder and CSVsample folder with Tot_counts.txt and Tot_percentage.txtand with absolute and percentage frequency and log.txt with analysis execution information. 
 
 **Graphics output**
 h5ad file with UMAP and others graphical output could be open with Cellxgene ( https://chanzuckerberg.github.io/cellxgene/ ). 


