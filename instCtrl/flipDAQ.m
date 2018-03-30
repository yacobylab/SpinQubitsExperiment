function flipDAQ
% Change sign of DAQ. 
% function flipDAQ
global smdata; 
smdata.channels(chl('DAQ1')).rangeramp(4)=smdata.channels(chl('DAQ1')).rangeramp(4)*-1;
end