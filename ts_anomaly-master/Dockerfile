FROM rocker/r-ver:3.4.2

# Make sure g++ uses c++11
ENV PKG_CXXFLAGS=-std=c++11

RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    libssl-dev \
    git \
    libcurl4-openssl-dev \
    python python-dev python-distribute python-pip \
    supervisor uwsgi uwsgi-plugin-python uwsgi-plugin-graylog2 \
  &&  . /etc/environment

# install R packages from repo defined in rocker/r-ver
#RUN export MRAN=$(cat /etc/environment | cut -d '=' -f 2) \
RUN MRAN=$(cat /etc/environment | cut -d '=' -f 2) \
  && Rscript -e "install.packages(c('Rcpp','changepoint','zoo', 'jsonlite'), repo = '$MRAN')" \
  && Rscript -e "install.packages(c('Rserve', 'Matrix', 'MASS'), rep = '$MRAN')"

# Pull isolation forest package from appfigures github repository
RUN git clone -b master https://d8b0c0ce331e083d44fc26227601b5a2c2718c9b:x-oauth-basic@github.com/appfigures/isofor.git /isofor \
  && cd /isofor \ 
  && R CMD INSTALL .

ADD tsAnomaly /tsAnomaly
RUN cd /tsAnomaly \
  && R CMD INSTALL .


# Add rserve user
RUN groupadd -r rserve && useradd -r -g rserve rserve

ADD start.R /start.R
ADD Rserv.conf /Rserv.conf

# Files required for python http wrapper
ADD requirements.txt /requirements.txt
ADD httpwrapper.py /httpwrapper.py

# Install required python libraries.
RUN pip install -r /requirements.txt

# startup script
ADD uwsgi.ini /uwsgi.ini
COPY supervisor-app.conf /etc/supervisor/conf.d/

# Set default to run the R server and HTTP wrapper.
ENTRYPOINT ["supervisord", "-n"]

# Port for Rserve
EXPOSE 6312

# Port for HTTP wrapper
EXPOSE 6316
