Model code developed in "Hydrology and Basal Composition Dominate the Variability of Basal Strength at Thwaites Glacier" by Paul T Summers, Cooper W Elsworth, Jenny Suckale.

This code was developed to understand the physical mechanism governing basal strength at Thwaites Glacier in Antarctica. Specifics are included in the original paper, and in comments within this code.

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
Smithe et al (2017) suppliment files
	https://doi.org/10.5194/tc-11-451-2017 (unzip in same directory as readme, replacing empty directory of the same name)

Input
This code can be in series using "MainHelper.m" or directly for 1 case using "ModelRunner.m".  

Output
Running "MainHelper.m" as configured now will produce basic figures for viewing the observed ice speed, computed ice speed, computed basal strength, and observed driving force. This will also save data files of the outputs in the data/ directory.
Once the data is saves, you can use "SaveFigs.m" to generate the figures as formatted in the final paper. 

References

Summers P T, Elsworth C W, and J Suckale (in review) Hydrology and Basal Composition Dominate the Variability of Basal Strength at Thwaites Glacier. Journal of Geophysical Research: Earth Surface.