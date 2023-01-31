FROM continuumio/miniconda3
RUN apt-get -y install git
RUN apt-get install g++
RUN git clone https://github.com/luglilab/Cytophenograph
SHELL ["/bin/bash", "--login", "-c"]
RUN conda install -c conda-forge python=3.9
RUN conda install -c conda-forge hnswlib
RUN pip install igraph
RUN pip install pyVIA
RUN pip install scanpy==1.9.1
RUN pip install phenograph==1.5.7
RUN pip install -e ./Cytophenograph/FlowSOM_LugliLab
RUN pip install openpyxl
RUN conda install -c conda-forge matplotlib
RUN conda install -c conda-forge gcc
RUN conda install -c conda-forge python-annoy
RUN pip install scanorama==1.6
RUN pip install fcsy
#RUN conda install -c conda-forge r-base=4.2.1
RUN apt-get update && apt-get install -y r-base
RUN conda install -c anaconda libgcc-ng
RUN Rscript -e "install.packages('flowAI')"
#RUN conda install -c bioconda bioconductor-flowai=1.24.0
RUN conda install scikit-image
#ADD /Users/simonepuccio/Documents/GitHub/Cytophenograph /usr/local/bin
CMD python /Cytophenograph/cytophenograph.v6.py
 
