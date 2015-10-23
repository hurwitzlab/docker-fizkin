FROM perl:latest

MAINTAINER Ken Youens-Clark <kyclark@email.arizona.edu>

COPY bin /usr/local/bin/

COPY scripts /usr/local/bin/

COPY include /usr/local/include/

COPY lib /usr/local/lib/

COPY share /usr/local/share/

ENV LD_LIBRARY_PATH=/usr/local/lib

ENV PERL5LIB=/usr/local/lib

CMD "run-fizkin"
