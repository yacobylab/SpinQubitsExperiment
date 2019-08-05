function startupsm
% Run when start matlab on experiment computer. 
%% Turn off annoying warnings. 
warn{1} = 'MATLAB:legend:IgnoringExtraEntries';
warn{2} = 'MATLAB:Java:GenericException';
warn{3} = 'MATLAB:polyfit:PolyNotUnique'; 
warn{4} = 'MATLAB:legend:PlotEmpty'; 
warn{5} = 'MATLAB:nearlySingularMatrix';         
warn{6} = 'stats:nlinfit:ModelConstantWRTParam';
warn{7} = 'MATLAB:rankDeficientMatrix'; 
warn{8} = 'stats:nlinfit:Overparameterized'; 
warn{9} = 'stats:nlinfit:IterationLimitExceeded'; 
warn{10} = 'stats:nlinfit:IllConditionedJacobian'; 
warn{11} = 'MATLAB:loadlibrary:FunctionNotFound';
warn{12} = 'MATLAB:dispatcher:pathWarning'; 
for i = 1:length(warn) 
    warning('off',warn{i}); 
end
%% Load smdata and open instruments. 
global smdata; global qdata; 
tuning = 0; 
% This is supposed to make notifications better. 
%system_dependent('RemotePathPolicy', 'TimecheckDirFile');
%system_dependent('RemoteCWDPolicy', 'TimecheckDirFile');

load z:/qDots/data/mx50structs/smdata_2018_12_12; % Load the rack 
logsetfile(smdata.files.log);
smdata.inst(inl('AWG1')).data.inst.RemoteHost='140.247.189.243'; % Matlab seems to overwrite the IP of the awg. this will set it correctly.
olist={'DAC1','DAC2','AWG1','AWG2','MercuryIPS','stepAtten','Hittite'};%,'N5183'}; 
load(smdata.files.scandata);
global scandata;
if tuning
   olist = [olist {'LockinA','LockinB','DMM1'}];
   smdata.inst(inl('LockinA')).data.inst.InputBufferSize = 4e5; 
   smdata.inst(inl('LockinB')).data.inst.InputBufferSize = 4e5; 
end
olist = sminstlookup(olist); 

for i=1:length(olist)
    try
        smopen(olist(i));
    catch
        fprintf('Error opening device %d: %s\n',olist(i),smdata.inst(olist(i)).name);
    end
end
try 
    smadachandshake; % Will throw an error if the handshakes don't match. 
catch 
    warning('DAC handshake failure'); 
end
sminitdisp; % Initialize channel display
openLabBrick % Start lab bricks. 
%% Populate channel list, start DAQ, load remaining structs. Plot recent tuning data. 
try
  smget(1:19);  
  smget({'RFpow2','RFpow1'}); 
  smget({'Phase1','Phase2'}); 
catch
  warning('Error populating initial channel list');
end
try 
    plssetup
catch 
    warning('Error loading pulses and awgdata'); 
end
try
    configATS660
catch 
    warning('Error starting DAQ'); 
end
cd(smdata.files.dir); 
global fbdata; %#ok<*NUSED>
load(smdata.files.fbdata);
gatesList = {'1a','2a','1b','2b','N12','T12','3a','4a','3b','4b','N34','T34','SD1top','SD1mid','SD1bot','SD4top','SD4mid','SD4bot','VRes','VBias'};
readoutList = {'Phase1','Phase2','samprate','RFfreq2','RFpow2','RFfreq1','RFpow1','RFfreq3','RFpow3'}; 
if tuning
    tuningList = {'LockExcA','LockFreqA','LockinTauA','LockSensA','LockExcB','LockFreqB','LockinTauB','LockSensB'};
else
    tuningList = {};
end
smdata.configch = [gatesList, readoutList,tuningList,'Time']; 
smset('samprate',1e7); 
load(smdata.files.tunedata);
load(smdata.files.qdata); 
global tuneData
if ~tuning 
    tuneData.rePlot; 
end
tuneData.updateAll; 
end
