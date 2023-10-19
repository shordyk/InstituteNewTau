
tautemp = ncread("~/Documents/MATLAB/ISSM/JPL1_ISSM_ctrl/SAMEDITstrbasemag_AIS_JPL1_ISSM_ctrl.nc","strbasemag")

importfile("~/Desktop/working_research/Institute/newbase.mat")
 



%tautemp(1,401:422,255:276) = newbase

for i=1:762
   tautemp(i,401:422,255:276) = newbase
end


%tautemp(:,401:422,255:276) = newbase

ncwrite("~/Documents/MATLAB/ISSM/JPL1_ISSM_ctrl/SAMEDITstrbasemag_AIS_JPL1_ISSM_ctrl.nc", "strbasemag", tautemp)