FROM perl:latest

MAINTAINER Ken Youens-Clark <kyclark@email.arizona.edu>

COPY local /usr/local/

COPY scripts /usr/local/bin/

#COPY bin /usr/local/bin/

#COPY include /usr/local/include/

#COPY lib /usr/local/lib/

#COPY share /usr/local/share/

ENV LD_LIBRARY_PATH=/usr/local/lib

#RUN curl -L http://cpanmin.us | perl - App::cpanminus

#RUN cpanm Carton

#WORKDIR /usr/local/bin

#RUN pwd

#RUN carton install --deployment

#RUN perl Makefile.PL 
#RUN cpanm --installdeps .

ENV PERL5LIB=/usr/local/lib/perl5

CMD "run-fizkin.pl"
