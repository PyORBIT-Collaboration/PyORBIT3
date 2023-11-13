FROM ubuntu:latest
RUN apt-get update -y
RUN apt-get install -y  build-essential
RUN apt-get install -y python3 openmpi-bin openmpi-common libopenmpi-dev  zlib1g-dev libfftw3-dev\
    python3-distutils python3-venv
RUN apt-get install -y git
RUN apt-get install -y libpython3-dev
