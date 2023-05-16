Model code developed in for work in progress on Institute Ice Stream and Shear Margin interactions

Requirements
You must download and install the following: 
Distmesh 
	https://popersson.github.io/distmesh/index.html
CVX 
	http://cvxr.com/cvx/download/
MEaSUREs matlab plug in + data files
	https://www.mathworks.com/matlabcentral/fileexchange/47329-measures
BedMachine matlab plug in + data files
	https://www.mathworks.com/matlabcentral/fileexchange/69159-bedmachine
ALBMAP
	https://doi.pangaea.de/10.1594/PANGAEA.734145 (unzip this file in ALBMAP folder)
Input
This code can be in series using "MainHelper.m" or directly for 1 case using "ModelRunner.m".  

Output
Running "MainHelper.m" as configured now will produce basic figures for viewing the observed ice speed, computed ice speed, computed basal strength, and observed driving force. This will also save data files of the outputs in the data/ directory.
Once the data is saves, you can use "SaveFigs.m" to generate the figures as formatted in the final paper. 

