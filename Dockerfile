FROM python:3.6-jessie

RUN apt-get update
RUN apt-get install unzip

# # get vcftools
RUN mkdir vcftools
RUN git clone https://github.com/vcftools/vcftools.git /vcftools
RUN (cd /vcftools; ./autogen.sh; ./configure; make; make install)

# get samtools
RUN git clone https://github.com/samtools/htslib
RUN git clone https://github.com/samtools/samtools
RUN (cd /samtools; autoheader; autoconf -Wno-syntax; ./configure; make; make install)

# install requirements
COPY ./requirements.txt /app/requirements.txt
RUN pip3 install -r /app/requirements.txt

COPY . /app
WORKDIR /app

CMD /bin/bash
