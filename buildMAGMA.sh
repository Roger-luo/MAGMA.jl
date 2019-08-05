#!/bin/bash

# install OpenBLAS
cd ~/Downloads/OpenBLAS-0.3.6
make -j$(nproc)
sudo make install
# install MAGMA
cd ~/Downloads/magma-2.5.1-alpha1
cp make.inc-examples/make.inc.openblas make.inc
sed -i '$a export CUDADIR=/usr/local/cuda' ~/.bashrc
sed -i '$a export OPENBLASDIR=/opt/OpenBLAS/lib' ~/.bashrc
source ~/.bashrc
sudo OPENBLASDIR=/opt/OpenBLAS/lib CUDADIR=/usr/local/cuda make -j$(nproc)
sudo OPENBLASDIR=/opt/OpenBLAS/lib CUDADIR=/usr/local/cuda make install

