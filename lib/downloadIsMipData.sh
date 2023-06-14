#!/bin/bash

mkdir ~/Documents/
mkdir ~/Documents/MATLAB
mkdir ~/Documents/MATLAB/ISSM/

cd ~/Documents/MATLAB/ISSM/
wget https://zenodo.org/record/2651652/files/README.txt

wget https://zenodo.org/record/2651652/files/JPL1_ISSM_init.zip
wget https://zenodo.org/record/2651652/files/JPL1_ISSM_ctrl.zip

unzip *.zip
