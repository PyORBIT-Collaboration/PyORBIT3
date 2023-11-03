FROM quay.io/centos/centos:stream9
# RUN dnf -y update
RUN dnf group install -y "Development Tools"
RUN dnf install -y python3-devel openmpi-devel fftw3-devel

