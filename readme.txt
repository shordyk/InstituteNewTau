Model code developed in for work in progress on Institute Ice Stream and Shear Margin interactions

Requirements
You must download and install the following: 
Distmesh 
	https://popersson.github.io/distmesh/index.html
CVX 
	http://cvxr.com/cvx/download/
ALBMAP
	https://doi.pangaea.de/10.1594/PANGAEA.734145 (unzip this file in ALBMAP folder)

Matlab Add-ons required:
MEaSUREs https://www.mathworks.com/matlabcentral/fileexchange/47329-measures
    Follow Installation steps in An_Overview.m to get the data files
BedMachine https://www.mathworks.com/matlabcentral/fileexchange/69159-bedmachine
Antarctic Mapping Tools https://www.mathworks.com/matlabcentral/fileexchange/47638-antarctic-mapping-tools
Bedmap2 Toolbox for Matlab https://www.mathworks.com/matlabcentral/fileexchange/42353-bedmap2-toolbox-for-matlab
    Get the data from here https://nsidc.org/data/nsidc-0756/versions/3 
    Run bedmap2_download.m    
Curve Fitting Toolbox
Image Processing Toolbox
cbrewer is used for many colormaps, but not needed if you change the colormaps

Input:
This code can be in series using "MainHelper.m" or directly for 1 case using "ModelRunner.m".
First run distmesh.m to generate the mesh .mat file.  

Output:
Running "MainHelper.m" as configured now will produce basic figures for viewing the observed ice speed, computed ice speed, computed basal strength, and observed driving force. This will also save data files of the outputs in the data/ directory.

Limitations:
It currently is designed to run on a Mac/PC machine, or a linux server (like Sherlock) but plotting is disabled on linux. 
