FROM continuumio/miniconda3
RUN apt-get -y install git
RUN apt-get install g++
RUN git clone https://github.com/luglilab/Cytophenograph
SHELL ["/bin/bash", "--login", "-c"]
RUN conda install -c conda-forge python=3.7.3
RUN conda install -c conda-forge hnswlib
RUN pip install igraph
RUN pip install parc
RUN pip install scanpy==1.9.1
RUN pip install phenograph==1.5.7
RUN pip install -e ./Cytophenograph/FlowSOM_LugliLab
RUN pip install openpyxl
RUN conda install -c conda-forge matplotlib=3.3.4
RUN conda install -c conda-forge gcc
RUN conda install -c conda-forge python-annoy
RUN pip install scanorama==1.6
