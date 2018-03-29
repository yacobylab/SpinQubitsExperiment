function config_script
global smdata
%% Load DLL and board handle 
addpath z:/qDots/daq/4chandaq 
addpath z:/qDots/daq/4chandaq/Include
if libisloaded('ATSApi')
  unloadlibrary('ATSApi');
end

if ~alazarLoadLibrary
   error('Unable to load library');
end
boardh = calllib('ATSApi', 'AlazarGetBoardBySystemID', 1, 1);
smdata.inst(sminstlookup('ATS660')).data.handle = boardh;

rmpath z:/qDots/daq/4chandaq 
rmpath z:/qDots/daq/4chandaq/Include

%% configure standard parameters
daqfn('SetExternalTrigger', boardh,2, 1); %AC/DC (1/2), 5V/1V (0/1)
daqfn('SetTriggerOperation', boardh, 0, 0, 2, 2, 140, 1, 3, 1, 1);
% TriggerOp, TriggerEng1, Source1, Slp1, Lvl1, TriggerEng2, Source2, Slp2, Level2
% J low to high ,engine J, trigin, negative, 
% source: 0=chan a, 1=chan b, 2=external, 3=off
% slope: 1 = rising, 2=falling.
% changed level from 146 to 140 on 11/29/10 to fix trigger issues.
% changed level back to 160 on  11/2/11 to try to fix trigger issues
% operation: J, K, J | K, J & K, J ^ K, J & !K !J & K = 0:6
% Engine1: J/K = 0/1, Sourve1: CH1, CH2, EXT, NONE = 0:3, Slope: +/-=1/2,
% level (0:255 = -100:100
% same for engine 2

daqfn('ConfigureAuxIO',boardh,14,0);
daqfn('SetTriggerDelay', boardh, 0);
daqfn('SetTriggerTimeOut', boardh, uint32(0)); 

daqfn('SetBWLimit', boardh, 1, 1);   % approx. 20 Mhz.
daqfn('SetBWLimit', boardh, 2, 1);

% Set DAQ range
nchans = 2; inst = sminstlookup('ATS660'); 
rngVals = [.2 .4 .8, 2, 5, 8, 16]; % range of the channel in V
rngRef =  [6, 7, 9, 11, 12, 14, 18]; % Alazar Ref for each V.
for ch = 1:nchans
    [~, rngInd] = min(abs(rngVals - smdata.inst(inst).data.rng(ch)));
    daqfn('InputControl', boardh, ch,2, rngRef(rngInd), 2-logical(smdata.inst(inst).data.highZ(ch)));
    smdata.inst(inst).data.rng(ch) = rngVals(rngInd);
end
end