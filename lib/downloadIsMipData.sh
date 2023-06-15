#!/bin/bash

mkdir ~/Documents/
mkdir ~/Documents/MATLAB
mkdir ~/Documents/MATLAB/ISSM

cd ~/Documents/MATLAB/ISSM/
curl https://zenodo.org/record/2651652/files/README.txt --output README.txt

curl https://zenodo.org/record/2651652/files/JPL1_ISSM_init.zip --output JPL1_ISSM_init.zip
curl https://zenodo.org/record/2651652/files/JPL1_ISSM_ctrl.zip --output JPL1_ISSM_ctrl.zip

unzip JPL1_ISSM_init.zip
unzip JPL1_ISSM_ctrl.zip

rm *.zip