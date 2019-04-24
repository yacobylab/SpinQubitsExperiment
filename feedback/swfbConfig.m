function [scan, fbScan, measScanNoFb] = swfbConfig(scan)
% configures a scan for software feedback, where the feedback and measurement pulses are separate groups.
% [ scan, fbScan, measScanNoFb ] = swfbConfig(scan)
% Uses fbdata.swfb to get information for the feedback.

global smdata; global fbdata;
nqubits=0;
for i=1:length(scan.loops(1).getchan) %decide which channels to feedback on based on the measurement scan.
    if any(strfind(scan.loops(1).getchan{i},'DAQ'))
        nqubits=nqubits+1;
    end
end
fbdata.times = []; fbdata.pumpHist = []; fbdata.gradHist = []; fbdata.set = [];
datachans=scan.loops(1).getchan(1:nqubits);
ic=smchaninst(datachans(1));
daqInst=ic(1); daqChan=ic(2);
switch nqubits
    case 1
        switch datachans{1}
            case 'DAQ1'
                ind = 1; 
                fbGroup=fbdata.params(1).fbGroup; % This is the dbz measurement. 
                fbInit = fbdata.params(1).fbInit; 
            case 'DAQ2'
                ind = 2; 
                fbGroup=fbdata.params(2).fbGroup;
        end
    case 2
        fbGroup=fbdata.fbGroup;
end
measScanNoFb=scan; %make the measurement scan
fbScan=fConfSeq(fbGroup,struct('nloop',fbdata.nloopfb,'nrep',1,'opts','raw','datachan',{datachans})); %make the feedback scan: this will run as a prefunction
fbScan.configch=[];
fbScan.loops(1).setchan='count2';
fbScan.consts(1) = [];
fbScan.figure = 1111; 
fbScan.disp = [];
xv=fbScan.data.pulsegroups(1).varpar(:,1)';
fbScan.xv=xv;
fbScanInit = fConfSeq(fbInit,struct('nloop',fbdata.nloopfb,'nrep',1,'opts','raw','datachan',{datachans})); %make the feedback scan: this will run as a prefunction

noMask=~isempty(strfind(scan.loops(1).prefn(2).fn,'PulseLine')); %See if the scan uses a mask.
% remove arm from the config fn. 
scan.configfn(1).args{1} = 'fast pls'; 
if ~noMask
    %incoming prefns are (1) arm (2) setmask (3) set pulseline
    %the outgoing prefunctions will be (1) arm (2) swfb, which will configure and set the mask, (3) re configure the card (4) rearm (5) setmask, (6) set pulse line
    if size(scan.loops(1).prefn,2)==3
        scan.loops(1).prefn(3:5)=scan.loops(1).prefn(1:3); %the scan came from fConfSeq2 or fConfSeq2_v2
    elseif size(scan.loops(1).prefn,2)==4
        scan.loops(1).prefn(3:6)=scan.loops(1).prefn(1:4); %the scan came from fConfSeq2_v4. There is an additional prefunction to tell the card how many pulses are in the group
    else
        error('Unexpected number of prefunctions');
    end        
    preNum = 2; 
else
    %some scans, like a T1 scan, have no mask and the incoming prefns are (1) arm (2) set pulseline
    %the outgoing prefunctions will be (1) arm (2) swfb, which will configure and set the mask, (3) get rid of the mask (4) re configure the card (5) rearm  (6) set pulse line
    scan.loops(1).prefn(4:5)=scan.loops(1).prefn(1:2);
    scan.loops(1).prefn(2).fn=@smaconfigwrap;
    scan.loops(1).prefn(2).args{1}=smdata.inst(daqInst).cntrlfn;
    scan.loops(1).prefn(2).args{2}=[daqInst daqChan 6];
    scan.loops(1).prefn(2).args{3}=[];            
    preNum=3;
end

% need to reconfigure the card after feedback to accept different readout/nloop etc.
args=scan.configfn(1).args{3}; %configfn(1) should be smabufconfig. need element 1 in case of raw mode.
val=args(1); rate=args(2);
fnStr=func2str(smdata.inst(daqInst).cntrlfn);
prefnStr=sprintf('@(x,daqInst,daqChan,val,rate,npulse,pls) %s([daqInst,daqChan, 5], val, rate,npulse,pls)',fnStr);    
% This is the prefn that sets the gradient. 
scan.loops(1).prefn(1).fn=@setgradient;
scan.loops(1).prefn(1).args={}; %get rid of the args that were there from the previous prefn.
scan.loops(1).prefn(1).args{1}=fbScan;
scan.loops(1).prefn(1).args{2}='';
scan.loops(1).prefn(1).args{3}=struct('figure',0,'tol',fbdata.params(ind).tol);

scan.loops(1).prefn(preNum).fn=prefnStr;
scan.loops(1).prefn(preNum).args{1}=daqInst;
scan.loops(1).prefn(preNum).args{2}=daqChan;
scan.loops(1).prefn(preNum).args{3}=val;
scan.loops(1).prefn(preNum).args{4}=rate;
scan.loops(1).prefn(preNum).args{6}=scan.data.pulsegroups(1).npulse(1);
scan.loops(1).prefn(preNum).args{5}='pls';
%add a configfn to set the gradient before the scan.
%This is useful for some scans, where all the loop variables are not 1 on the first loop.
scan.configfn(end+1).fn=@smaconfigwrap;
scan.configfn(end).args{1}=@setgradient;
scan.configfn(end).args{2}='dummy';
scan.configfn(end).args{3}=fbScanInit;
scan.data.setpt=[fbdata.params(:).setpt];
scan.cleanupfn(end+1).fn = @saveFBdata; 
scan.cleanupfn(end).args = {}; 
end