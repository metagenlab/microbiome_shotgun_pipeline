FROM ubuntu:18.04
RUN apt-get update && apt-get install -y wget bash zip
WORKDIR /usr/local/bin 
RUN wget https://github.com/arpcard/rgi/archive/5.1.0.zip && unzip 5.1.0.zip
RUN mv rgi-*/* .

RUN wget http://github.com/bbuchfink/diamond/releases/download/v0.8.36/diamond-linux64.tar.gz 
RUN tar xvf diamond-linux64.tar.gz 
RUN mv diamond /usr/bin 
RUN rm diamond-linux64.tar.gz

RUN apt-get install -y python3 python3-dev python3-pip ncbi-blast+ prodigal bowtie2 samtools bedtools bamtools

RUN pip3 install -r requirements.txt && \
    pip3 install .
 
RUN bash test.sh   

RUN apt-get clean

CMD ["rgi"]