FROM continuumio/miniconda3:4.7.10
RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda
RUN conda install python=3.6 numpy theano keras  scikit-learn biopython
RUN git clone https://github.com/jessieren/DeepVirFinder
RUN chmod +x /DeepVirFinder/dvf.py
ENV PATH /DeepVirFinder:$PATH
ENTRYPOINT ["dvf.py"]
CMD ["-h"]
