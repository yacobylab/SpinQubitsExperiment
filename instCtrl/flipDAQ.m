function flipDAQ(chanName)
% Change sign of DAQ. 
% function flipDAQ
global smdata; 
if ~exist('chanName','var'), chanName = 'DAQ1'; end
smdata.channels(chl(chanName)).rangeramp(4)=smdata.channels(chl(chanName)).rangeramp(4)*-1;
end