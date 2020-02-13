FROM python:3.8.1-alpine3.10
LABEL maintainer="medcelerate@gmail.com"


RUN apk add --no-cache libstdc++
RUN apk add --no-cache python3-dev
RUN apk add --no-cache --virtual .build-deps g++
RUN ln -s /usr/include/locale.h /usr/include/xlocale.h

COPY requirements.txt requirements.txt

RUN apk add --no-cache python3 && \
    python3 -m ensurepip && \
    pip3 install --upgrade pip

RUN pip3 install -r requirements.txt
RUN apk del .build-deps
RUN apk update && apk add bash

COPY farnsworth.py /bin/Farnsworth
RUN chmod +x /bin/Farnsworth
