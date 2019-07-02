function plssetup
% Load updated pulses, awgdata, create awg inst. 
% Fixme: why is this still hardcoded... 

global awgdata; global plsdata; global smdata; 

if strcmp(computer, 'GLNX86') || strcmp(computer,'PCWIN64')
    plsdata.datafile = 'z:/qDots/awg_pulses/plsdata_2017_04_14.mat';    
else
    if exist('z:/qDots','file')
      plsdata.datafile = 'z:/qDots/awg_pulses/plsdata_2017_04_14.mat';
    else
      plsdata.datafile = 'y:/qDots/awg_pulses/plsdata_2017_04_14.mat';  
    end
end
plssync('load');

% Hack that only makes sense in our setup. 
awgList = sminstlookup('AWG5000'); 
awgloaddata;
 if exist('smdata', 'var') && isfield(smdata, 'inst')
     for i = 1 :1% length(awgList)
         awgdata(i).awg = smdata.inst(awgList(i)).data.inst;
     end
%     if ~isempty(sminstlookup('AWG7000')) && length(awgdata)>i
%         awgdata(i+1).awg = smdata.inst(sminstlookup('AWG7000')).data.inst;
%     end
% end

end