function testAd(opts)
% Check that adiabatic prep and read works. Vary epsilon (voltage) or time, 
% perform just adprep or adprep/adread (combines to 4 scans). 
% function testAd(opts)
% Run with ampok, so that can use with nonzero measurement point.
% For eps and only prep, will end up with mixed state once we reach sep, but more singlets when eps too small.  
% For eps and prep/meas, should always get singlet. 
% For time and prep/meas, will end up with singlet when time large enough,
% mixed state when it's too fast. 
% For time and just prep, I think we'll mostly get mixed. 
% The epsilon of the time scan and time of the epsilon scan should be in
% safe region, but check them if things not working. 
% Need to make sure that ad skips STP. Edit start point in dictionary so
% that 
global tuneData;
if ~exist('opts','var'), opts = ''; end
% Define pulses
t = linspace(0,1,50);
eps = linspace(1,8,100);
    
pg.ctrl='loop pack';
pg.chan=[getNum(tuneData.xyChan{1}),getNum(tuneData.xyChan{2})];
pg.dict=tuneData.activeSetName;

pg.pulses=32;
pg.name = sprintf('adTestTimeMeas%s',upper(pg.dict(1)));
%Parameters: pulselength, ramptime. 
pg.params=[9 0];
pg.varpar =t';
plsdefgrp(pg);    
name = {pg.name};

pg.pulses=121;
pg.name = sprintf('adTestTime%s',upper(pg.dict(1)));
plsdefgrp(pg);
name{2} = pg.name;

pg.pulses=115;
%Parameters: pulselength, eps, evo
pg.name = sprintf('adTestPosMeas%s',upper(pg.dict(1)));
pg.params=[7 1 0];
pg.varpar =eps';
plsdefgrp(pg);

name{3} = pg.name;
pg.pulses=112;
pg.name = sprintf('adTestPos%s',upper(pg.dict(1)));
plsdefgrp(pg);
name{4} = pg.name;
if ~isopt(opts,'run') % Add to AWG
    awgadd(name);
    awgcntrl('on start wait err');
end
%% Run scan
nloop = 1000; nrep = 10;
fignum = 400; figure(fignum+1); clf; hold on;
figure(fignum+2); clf; hold on;
pars ={t,t,eps,eps};
for i = 1:4
    scan = fConfSeq(name{i},{'nloop',nloop,'nrep',nrep, 'datachan',tuneData.dataChan,'opts','ampok'});
    scan = measAmp(scan); 
    data=smrun(scan,smnext(sprintf('%s_L',name{i})));
    meanData = nanmean(data{1});
    figure(fignum+ceil(i/2));
    plot(pars{i},meanData)
    title(name{i});
end
sleep
end