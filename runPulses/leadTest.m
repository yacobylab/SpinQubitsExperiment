
global tuneData
%% Make a lead scan based on load, then one on the other lead, then load scan. 
pg.pulses=19;
pg.dict=tuneData.activeSetName; 
side = upper(tuneData.activeSetName(1)); 
pg.chan=[getNum(tuneData.xyChan{1}),getNum(tuneData.xyChan{2})];
pg.ctrl='loop pack';  

npulse=200;
pulseLen=25;
rng = 3; 
% 
% if tuneData.sepDir(1)>0  % TL meas point. Set up direction of lead.
%     leadDir = -[1 tuneData.loadPos.slope];
%     loadDir = [1 1./leadDir(2)]; loadDir=-sqrt(2)*loadDir/(norm(loadDir));
% else  % BR meas point
%     leadDir = [1 tuneData.loadPos.slope];
%     loadDir = [1 -1./leadDir(2)]; loadDir=sqrt(2)*loadDir/(norm(loadDir));
% end
% leadDir = leadDir / norm(leadDir); % lead Direction, towards the right
% if isopt(opts,'vert') 
%     loadDir = [0,1];
% elseif isopt('opts','horz') 
%     loadDir = [1,0]; 
% end

juncDist = tuneData.chrg.trTriple(end,:)-tuneData.chrg.blTriple(end,:);

%pulse length, measure time, %measure location (2)
pg.params=[pulseLen pulseLen-1 0 0];
dir = {'x','y'}; 
loadDir = [0,1; 1,0]; 
dist = 3;
for i = 1:2
    leadDir = [1,tuneData.chrg.([dir{i} 'LeadSlope'])(end)]; leadDir = leadDir / norm(leadDir);
    if i ==1, leadDir = -leadDir; end
    loadCenter = -1e3*(1/2*juncDist+tuneData.chrg.defaultOffset) + leadDir * dist;    
    measLocs = loadCenter + linspace(-0.5,0.5,npulse)' * rng*loadDir(i,:)-tuneData.measPt;
    pg.varpar = measLocs; % Range of scan from the current reload val
    pg.name=sprintf('leadPos%s%s',upper(dir{i}),side);
    plsdefgrp(pg);
    leadGrp{i}=pg.name;    
end
awgadd(leadGrp);
%% Make a line scan 
dict = pdload(tuneData.activeSetName); 
rng = 2; center = 0; 
npulse=200;
pulseLen = 10; 
eps = linspace(-0.5,0.5,npulse)'.*rng+center * dict.exch.val-tuneData.measPt;
%pulse length, measure time, %measure location
pg.params=[pulseLen, pulseLen-1, 0, 0];
pg.name=sprintf('Line%s',side);
plsdefgrp(pg); 
awgadd(pg.name); 

scans{1}=fConfSeq(lineGrp,struct('nloop',1000,'nrep',15,'opts','','datachan',tuneData.dataChan));
scans{2}=fConfSeq(leadGrp{1},struct('nloop',200,'nrep',15,'opts','','datachan',tuneData.dataChan));
scans{3}=fConfSeq(leadGrp{2},struct('nloop',200,'nrep',15,'opts','','datachan',tuneData.dataChan));

namePat{1}='line%s_%1.1fmK'; 
namePat{2}='leadPosX%s_%1.1fmK'; 
namePat{3}='leadPosY%s_%1.1fmK';
awgcntrl('on start err');
%% Test scans before starting
mcTemp=cell2mat(smget('MC'));
fprintf('Temp is %1.1f \n',mcTemp*1e3);
for j = 1:length(scans)
    scanName = sprintf(namePat{j},side,mcTemp*1e3);
    smrun(scans{j},smnext(scanName));
    sleep('fast');
end
sleep('fast');
%%
temps = [60, 65, 70,80,90,100,120,140,170,200,250,300]*1e3;
for i = 1:length(temps)
    setMCSP(temps(i));
    MC = cell2mat(smget('MC'));
    err = (MC-temps(i))/temps(i);
    while err > 0.03 % Wait for things to settle.
        pause(60);
        MC = cell2mat(smget('MC'));
        err = (MC-temps(i))/temps(i);
    end
    fprintf('Temp is %1.1f \n',MC*1e3);
    for j = 1:length(scans)
        scanName = sprintf(namePat{j},side,mcTemp*1e3);
        smrun(scans{j},smnext(scanName));
        sleep('fast');
    end
end