FROM continuumio/miniconda3:4.7.10
RUN conda config --add channels anaconda
RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda
RUN apt install -y g++
RUN conda install python=3.6 numpy=1.17.0 theano=1.0.4 keras=2.2.5  scikit-learn=0.21.3 biopython=1.72  mkl=2019.4
RUN git clone https://github.com/jessieren/DeepVirFinder
RUN chmod +x /DeepVirFinder/dvf.py
RUN conda clean --all
ENV PATH /DeepVirFinder:$PATH
ENTRYPOINT ["dvf.py"]
CMD ["-h"]
