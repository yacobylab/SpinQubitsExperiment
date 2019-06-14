function testSep(config)
% Perform a ramsey type scan with no readout (i.e. a dBz scan but at variable epsilon). 
% Sweep epsilon.
% If J is still large, we will have a singlet state rotating around Z, so
% nothing will happen. As J turns off, we will have a singlet state
% rotating around X. Because dBz not set, will just create random mixed
% state, sig ~= 0.5.

global tuneData
pg.pulses=122;
pg.ctrl='loop pack';

pg.chan=[getNum(tuneData.xyChan{1}),getNum(tuneData.xyChan{2})];
pg.dict={tuneData.activeSetName};

%Parameters: pulselength, eps, evo
pg.dict={struct('prep',struct('type','@null'),'read',struct('type','@null')),pg.dict};
eps = linspace(0,4,12);
evo = 1:128;
pg.varpar = evo';
pg.trafofn.func=@rc_trafofn; pg.trafofn.args=0;
for i = 1:length(eps)
    pg.params=[5 eps(i) 0];
    pg.name = sprintf('sepTest%s%d',upper(tuneData.activeSetName(1)),i);
    plsdefgrp(pg);
    sep{i}=pg.name; %#ok<*AGROW>
end
awgadd(sep);
awgcntrl('on start wait err');

scan=fConfSeq(sep,struct('nloop',200,'nrep',8,'datachan',tuneData.dataChan));
scan = measAmp(scan); 
data = smrun(scan,smnext(sprintf('sepTest%s',upper(tuneData.activeSetName(1)))));
sleep('fast');
figure(1120); clf;
meanD = squeeze(nanmean(data{1}));
imagesc(evo,eps,meanD);
xlabel('Time'); ylabel('Epsilon');
set(gca,'YDir','Normal');
end