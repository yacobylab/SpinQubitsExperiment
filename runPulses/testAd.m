function testAd
% Check that adiabatic prep and read works. Vary voltage or time, perform just adprep or
% adprep/adread (combines to 4 scans). 
% Run with ampok, so that can use with nonzero measurement point. 
% For voltage/only prep, will end up with mixed state once we reach sep, singlet too small. 
% For time/meas, will end up with singlet with large time, mixed with small

global tuneData;
%% Define pulses 
 
pg.ctrl='loop pack';
pg.chan=[str2double(char(regexp(tuneData.xyChan{1},'\d+','match'))),str2double(char(regexp(tuneData.xyChan{2},'\d+','match')))];
pg.dict=tuneData.activeSetName;

%Parameters: pulselength, eps, evo
pg.pulses=32;
t=linspace(0,1,50);
pg.name = sprintf('adTestTime%s',upper(pg.dict(1)));
pg.params=[6 0];
pg.varpar =t';
plsdefgrp(pg);

name = {pg.name};
pg.pulses=33;
pg.name = sprintf('adTestTimeMeas%s',upper(pg.dict(1)));
plsdefgrp(pg);
name{2} = pg.name;

pg.pulses=115;
eps=linspace(1,8,100);
%Parameters: pulselength, eps, evo
pg.name = sprintf('adTestPos%s',upper(pg.dict(1)));
pg.params=[6 0.6 0];
pg.varpar =eps';
plsdefgrp(pg);

name{3} = pg.name;
pg.pulses=112;
pg.name = sprintf('adTestPosMeas%s',upper(pg.dict(1)));
plsdefgrp(pg);
name{4} = pg.name;
%% Run scan
awgadd(name);
awgcntrl('on start wait err');
nloop = 400; nrep = 10;
fignum = 400; figure(fignum+1); clf; hold on;
figure(fignum+2); clf; hold on;
pars ={t,t,eps,eps};
for i = 1:4
    scan = fConfSeq(name{i},{'nloop',nloop,'nrep',nrep, 'datachan',tuneData.dataChan,'opts','ampok'});
    scan.consts(end+1).setchan=tuneData.xyChan{1};
    scan.consts(end).val=tuneData.measPt(1);
    scan.consts(end+1).setchan=tuneData.xyChan{2};
    scan.consts(end).val=tuneData.measPt(2);
    data=smrun(scan,smnext(sprintf('%s_L',name{i})));
    meanData = nanmean(data{1});
    figure(fignum+ceil(i/2));
    plot(pars{i},meanData)
    title(name{i});
end
sleep
end