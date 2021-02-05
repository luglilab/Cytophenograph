# Cytophenograph

Identifies subpopulations in high-dimensional single-cell data. Cytophenograph is a computational pipeline that was developed to avoid the disadvantages of manual gating. The pipeline is developed using Python3, the clustering method adopted is a custom version of Phenograph (https://github.com/jacoblevine/PhenoGraph) where we insert a blocked seed. Besides Phenograph pipeline needs the following package installed: Pandas,Numpy,Sklearn for data parsing and exploring and openTSNE,Seaborn,Matplotlib for data visualization. This method is adaptative both in terms of dimensionality and sample size, making it suitable in a range of settings for which single-cell population structure is of interest, including other cancers or healthy tissues, and for use with other emerging single-cell technologies.

## 1) Installation 

Install Miniconda
Miniconda is a Python distribution, package manager, and virtual environment solution. We recommend installing Miniconda with Python 3 (miniconda3), as many bioinformatics packages are now transitioning to Python 3. You can still install Python 2 software with miniconda3 by passing the python=2.7 flag when you create a new environment; otherwise the default Python version will be Python 3.

Begin by downloading Miniconda and following the associated installation instructions.

https://docs.conda.io/en/latest/miniconda.htm

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


```python
conda env create -n cytophenograph -f ./Cytophenograph/environment.yml 
conda activate cytophenograph
```

Install Phenograph


```python
git clone https://github.com/luglilab/Phenograph_LugliLab
```

Move on Phenograph folder



```python
pip install -e Phenograph_LugliLab
```

```python
pip install scipy==1.4.1 --use-feature=2020-resolver
```

```python
conda install -c conda-forge umap-learn
```


```python
cd ../Cytophenograph/
python cytophenograph.py --help
```


Test Execution 
```python
mkdir output
python cytophenograph.py -i ./CD8_Panel_II_channelvalues_GA/ -o ./output/ -k 300 -m ./marker.txt -n TestCytophenograph -t 10
```
