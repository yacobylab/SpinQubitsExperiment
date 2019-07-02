function [scan, fbScan, measScanNoFb] = swfbConfig(scan)
% Configures a scan for software feedback, where the feedback and measurement pulses are separate groups.
% [scan, fbScan, measScanNoFb] = swfbConfig(scan)
% Creates dBz meas scan. 
% Before, prefns: (1) arm (2) set mask (3) set pulseline
% After: (1) arm for swfb (2) swfb (3) reconfigure DAQ (4) arm (5) set mask, (6) set pulse line       

global smdata; global fbdata;
nqubits=0;
for i=1:length(scan.loops(1).getchan) %decide which channels to feedback on based on the measurement scan.
    if any(strfind(scan.loops(1).getchan{i},'DAQ'))
        nqubits=nqubits+1;
    end
end
fbdata.times = []; fbdata.pumpHist = []; fbdata.gradHist = []; fbdata.set = []; fbdata.err = [];
datachans=scan.loops(1).getchan(1:nqubits);
ic=smchaninst(datachans(1));
daqInst=ic(1); daqChan=ic(2);
switch nqubits
    case 1
        ind = str2double(datachans{1}(end)); % determine side        
    case 2
        ind = 3; % Add this        
end
fbGroup = fbdata.params(ind).fbGroup; % This is the dbz measurement.
fbInit = fbdata.params(ind).fbInit;
measScanNoFb=scan; % Scan before swfb applied.

% Make the dBz measurement scan for feedback, part of prefn that sets gradient.
fbScan=makeFeedbackScan(fbGroup,fbdata.params(ind).nloopfb,{datachans});
fbScan.disp = [];
fbScanInit = fConfSeq(fbInit,struct('nloop',fbdata.params(ind).nloopfb,'nrep',1,'opts','raw','datachan',{datachans}));

noMask=~isempty(strfind(scan.loops(1).prefn(2).fn,'PulseLine')); % See if the scan uses a mask.
% remove arm from the config fn. 
scan.configfn(1).args{1} = 'fast pls'; 

% prefn 1: set gradient: 
prefn(1).fn=@setgradient;
prefn(1).args{1}=fbScan;
prefn(1).args{2}='';
prefn(1).args{3}=struct('figure',0,'tol',fbdata.params(ind).tol,'nloop',fbdata.params(ind).nloopfb);

if noMask % Need to clear mask after fb.    
    prefn(2).fn = @smaconfigwrap;
    prefn(2).args{1} = @clearMask;
    prefn(2).args{2} = datachans{1};          
end
% Need to reconfigure the card after feedback to accept different readout/nloop etc.
args=scan.configfn(1).args{3}; %configfn(1) should be smabufconfig. need element 1 in case of raw mode.
val=args(1); rate=args(2);
fnStr=func2str(smdata.inst(daqInst).cntrlfn);
prefn(end+1).fn=sprintf('@(x,daqInst,daqChan,val,rate,pls,npulse) %s([daqInst,daqChan, 5], val, rate,pls,npulse)',fnStr);
prefn(end).args{1}=daqInst;
prefn(end).args{2}=daqChan;
prefn(end).args{3}=val;
prefn(end).args{4}=rate;
prefn(end).args{5}='pls';
prefn(end).args{6}=scan.data.pulsegroups(1).npulse(1);

% Shift extant prefns back to after new swfb ones. 
nprefn = length(scan.loops(1).prefn); 
prefn(end+1:end+nprefn)=scan.loops(1).prefn;
scan.loops(1).prefn = prefn;     
% Add a configfn to set the gradient before the scan.
% This is useful for some scans, where all the loop variables are not 1 on the first loop.

scan.configfn(end+1).fn=@smaconfigwrap;
scan.configfn(end).args{1}=@setgradient;
scan.configfn(end).args{2}='dummy';
scan.configfn(end).args{3}=fbScanInit;
scan.data.setpt=[fbdata.params(ind).setpt];
scan.cleanupfn(end+1).fn = @saveFBdata; 
scan.cleanupfn(end).args = {}; 
end