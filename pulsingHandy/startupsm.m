function startupsm

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

global smdata;
tuning = 0; 
load z:/qDots/sm_config/smdata_MX50_2016_06_13; % load the rack 
logsetfile('z:/qDots/notes/log_2015_11_05.txt');
smdata.inst(inl('AWG1')).data.inst.RemoteHost='140.247.189.243'; %Matlab seems to overwrite the IP of the awg. this will set it correctly.
olist={'DAC1','DAC2','AWG1','AWG2','MercuryIPS','stepAtten'};%,'N5183'}; 

load('z:/qDots/data/data_2015_11_05/scandata_2016_02_12');
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
    warning('DAC handshake failure \n'); 
end
sminitdisp; %initialize channel display

try
    if ~libisloaded('vnx_fsynth')
        [success,warnings]=lbLoadLibrary2;
        if ~success
            error('Unable to load vnx_fsynth');
        end
    end
    labBricks = inl('LabBrick');
    calllib('vnx_fsynth','fnLSG_SetTestMode',false);
    brickfn('GetNumDevices'); % this needs to be run first. 
    [~,devIDs]=calllib('vnx_fsynth','fnLSG_GetDevInfo',uint32(zeros(1,length(labBricks))));
    
    for i = 1:length(labBricks)
        serialNum(i) = brickfn('GetSerialNumber',devIDs(i));      %#ok<AGROW>
    end
    for i = 1:length(labBricks)
        devIDCurr = devIDs(serialNum==smdata.inst(labBricks(i)).data.serial);
        smdata.inst(labBricks(i)).data.handle = devIDCurr;
        brickfn('InitDevice',devIDCurr);
    end
catch
    warning('Error initializing lab bricks');
end
try
  smget(1:19);  
  smget({'RFpow2','RFpow1'}); 
catch
  warning('Error populating initial channel list\n');
end
try 
    plssetup
catch 
    warning('Error loading pulses and awgdata'); 
end
try
    config_script
catch 
    warning('Error starting DAQ'); 
end
cd z:/qDots/data/data_2015_11_05
global fbdata; %#ok<*NUSED>
load z:/qDots/sm_config/fbdata_2017_04_28

load tunedata_2016_04_11.mat
global tuneData

tuneData.rePlot; 
end