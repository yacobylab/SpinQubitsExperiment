function testAd
% Check that adiabatic prep and read works, and at what voltages one must
% prepare to, what times.
% Run 4 scans. One with adiabatic prep and no read, one with adprep and
% adread. 

% 
global tuneData;

pg.pulses=32;
pg.ctrl='loop pack';
pg.chan=[str2double(char(regexp(tuneData.xyChan{1},'\d+','match'))),str2double(char(regexp(tuneData.xyChan{2},'\d+','match')))];
pg.dict=tuneData.activeSetName;
t=linspace(0,1,50);
%Parameters: pulselength, eps, evo
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