from setuptools import setup

setup(
    name='FlowSom',
    url='https://github.com/Hatchin/FlowSOM',
    author='Sangyu Shen',
    author_email='sangyushen@gmail.com',
    packages=['flowsom'],
    install_requires=['FlowCytometryTools', 'matplotlib', 'networkx', 'minisom', 'scikit-learn', 'numpy', 'pandas'],
    version='0.1.1',
    license='MIT',
    description='A Python implementation of FlowSOM algorithm for clustering and visualizing a mass cytometry data set.',
    keywords=['flowsom', 'cytometry', 'machine learning', 'visualization', 'clustering', 'dimentionality reduction']
)