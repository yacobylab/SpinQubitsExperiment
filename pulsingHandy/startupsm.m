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
global smdata;
tuning = 1; 
load z:/qDots/sm_config/smdata_MX50_2016_06_13; % load the rack 
logsetfile(smdata.files.log);
smdata.inst(inl('AWG1')).data.inst.RemoteHost='140.247.189.243'; %Matlab seems to overwrite the IP of the awg. this will set it correctly.
olist={'DAC1','DAC2','AWG1','AWG2','MercuryIPS','stepAtten'};%,'N5183'}; 
load(smdata.files.scandata);
global scandata;
if tuning         
   olist = [olist {'LockinA','LockinB','Switch','DMM1'}];
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
    smadachandshake; % will throw an error if the handshakes don't match. 
catch 
    warning('DAC handshake failure'); 
end
sminitdisp; %initialize channel display
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

load(smdata.files.tunedata);
global tuneData
if ~tuning 
    tuneData.rePlot; 
end
end