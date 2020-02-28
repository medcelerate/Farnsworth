FROM continuumio/miniconda:latest
LABEL maintainer="medcelerate@gmail.com"

COPY environment.yml .
RUN apt-get update -y && apt-get install -y gcc
RUN conda update -n base -c defaults conda && \
    conda env update --name base --file /environment.yml && \
	rm /environment.yml

COPY farnsworth.py /bin/Farnsworth
RUN chmod +x /bin/Farnsworth
