function testSep
% Perform a ramsey type scan with no readout (i.e. a dBz scan but at variable epsilon). Sweep epsilon.  
% If J is still large, 
global tuneData
pg.pulses=15;
pg.ctrl='loop pack';
pg.chan=[str2double(char(regexp(tuneData.xyChan{1},'\d+','match'))),str2double(char(regexp(tuneData.xyChan{2},'\d+','match')))];
pg.dict={tuneData.activeSetName};

%Parameters: pulselength, eps, evo
pg.dict={struct('prep',struct('type','@null'),'read',struct('type','@null')),pg.dict};
eps = linspace(2,7,10);
evo = 1:128;
pg.varpar = evo';
for i = 1:length(eps)
    pg.params=[4.5 eps(i) 0];
    pg.name = sprintf('sepTest%d',i);
    plsdefgrp(pg);
    sep{i}=pg.name;
end
awgadd(sep);
awgcntrl('on start wait err');

scan=fConfSeq(sep,struct('nloop',200,'nrep',20,'datachan','DAQ1'));
data = smrun(scan,smnext('sepTestL'));
sleep;
figure(1120); clf; 
meanD = squeeze(nanmean(data{1})); 
imagesc(evo,eps,meanD); 
xlabel('Time'); ylabel('Epsilon'); 
set(gca,'YDir','Normal'); 
end