#
# docker run --rm \
#        --interactive \
#        --tty \
#        -v "$(pwd):$(pwd)" \
#        -w $(pwd) \
#        -it muccg/bpasubmit
#
FROM python:3.6-alpine
LABEL maintainer "https://github.com/muccg"

ENV VIRTUAL_ENV /env
ENV PIP_NO_CACHE_DIR="off"
ENV PYTHON_PIP_VERSION 9.0.1
ENV PYTHONIOENCODING=UTF-8

RUN apk --no-cache add \
    ca-certificates \
    git

RUN python3 -m venv $VIRTUAL_ENV \
    && $VIRTUAL_ENV/bin/pip install --upgrade \
    pip==$PYTHON_PIP_VERSION
ENV PATH $VIRTUAL_ENV/bin:$PATH

COPY requirements.txt /requirements.txt
RUN pip install -r /requirements.txt

COPY . /app
RUN pip install -e /app

RUN addgroup -g 1000 bpa \
  && adduser -D -h /data -H -S -u 1000 -G bpa bpa \
  && mkdir /data \
  && chown bpa:bpa /data

USER bpa

ENTRYPOINT ["/env/bin/bpa-submit"]
