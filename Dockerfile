FROM ubuntu

RUN apt-get -qq update

ADD tzscript.sh /tzscript.sh
RUN chmod 777 /tzscript.sh; /tzscript.sh

RUN apt-get -qq -y install unzip vim tar sudo gcc g++ gfortran m4 python2.7 python-dev git wget bzip2 tar expect
RUN wget http://devtools.fg.oisin.rc-harwell.ac.uk/nightly/7.0/ccp4-7.0-linux64-latest.tar.bz2
RUN bunzip2 ccp4-7.0-linux64-latest.tar.bz2
RUN mkdir ./ccp4
RUN tar -xf ccp4-7.0-linux64-latest.tar -C ./ccp4 --strip-components=1

ADD ccp4.setup-sh ./ccp4/bin

ADD . /pandda
RUN chmod 775 /pandda/pandda-2_install.sh
RUN /pandda/pandda-2_install.sh

RUN mkdir /BAZ2B
RUN cd /BAZ2B; wget https://zenodo.org/record/48768/files/data.zip; /bin/bash -c "unzip /BAZ2B/data.zip"
RUN mkdir /BAZ2B_out
