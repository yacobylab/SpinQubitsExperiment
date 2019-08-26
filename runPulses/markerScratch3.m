% Create qdata
qdata(1).name = 'left'; 
qdata(1).chan = [2,1,5,6]; 
qdata(1).eps = 1; 

qdata(2).name = 'right'; 
qdata(2).chan= [3,4,7,8]; 
qdata(2).eps = 1; 
%% First IQ ramsey, varying phase
on = 2; off = 1;
name = qdata(on).name; 
clear pg; clear pgOff; clear pgOn; 
npulses=64;
plen=6;
ngrps = 20;

% Other side
offName = qdata(off).name; 
pgOff.pulses=132;
pgOff.ctrl='loop pack';
pgOff.chan = qdata(off).chan; 
offSide = upper(offName(1)); 
pgOff.name=sprintf('RabiToy%02d_%d_%s',plen,npulses,offSide);
pgOff.dict={struct('prep','@dbzprep','read','@dbzprep','measLoc','@wait'),offName};
pgOff.params=[plen, 0,qdata(off).eps, 1, 0, 10,5,0];
pgOff.varpar = (1:npulses)';
plsdefgrp(pgOff);
pg2=pgOff.name;

phaseVals = linspace(0,pi,ngrps); 
pgOn = pgOff; 
side = upper(name(1)); 
pgOn.dict={struct('prep','@dbzprep','read','@dbzread','measLoc','@meas'),name};
pgOn.chan = qdata(on).chan; 

pg.ctrl='grp loop pack';
pg.pulseind(2,:)=1:npulses;
pg.pulseind(1,:)=1:npulses;
pg.chan=[pgOn.chan, pgOff.chan];
pg.matrix=eye(8);

grps = {};
for i = 1:ngrps
    pgOn.name=sprintf('RabiToy%02d_%02d_%s',plen,i,side);
    
    iVal = cos(phaseVals(i));
    qVal = sin(phaseVals(i));
    %          PlsLen, wait time, exch, I, Q,  marker pre, marker delay, time
    pgOn.params=[plen, 0, qdata(on).eps, iVal, qVal, 10,5,0]; 
    pg1=pgOn.name;
    plsdefgrp(pgOn);
            
    pg.pulses.groups={pg1,pg2};
    pg.name=sprintf('IQToyLR_%02d',i);
    grps{end+1}=pg.name;      
    plsdefgrp(pg);
end
awgrm(fbdata.lastInd,'after'); 
awgadd(grps);
awgcntrl('on start wait err raw');
fprintf('ready!\n');
%% Run scan, single Ramsey
rabiScan=fConfSeq(grps,struct('nloop',128,'nrep',64,'opts','swfb','datachan',tuneData.dataChan));
smset('RFfreq3',1.45e9); 
smset('RFpow3',26); 

smrun(rabiScan,smnext(sprintf('RamseyIQ_%s',side)))
sleep
% best phase 0.6348 
%% Single ramsey scan at optimal phase
on = 1; 
name = 'left'; 
clear pg; clear pgOff; clear pgOn; 
npulses=64;
plen=6;

% Other side
off = 2;
offName = 'right';
pgOff.pulses=132;
pgOff.ctrl='loop pack';
pgOff.chan = qdata(off).chan;
pgOff.name=sprintf('RabiToy%02d_%d_%s',plen,npulses,offName);
pgOff.dict={struct('prep','@dbzprep','read','@dbzprep','measLoc','@wait'),offName};
pgOff.params=[plen, 0.11,qdata(off).eps, 1, 0, 10,5,0];
pgOff.varpar = (1:npulses)';
plsdefgrp(pgOff);
pg2=pgOff.name;

ngrps = 20;
phaseVals = linspace(0,pi,ngrps);
pgOn = pgOff;
pgOn.dict={struct('prep','@dbzprep','read','@dbzread','measLoc','@meas'),name};
pgOn.chan = qdata(on).chan;
pg.ctrl='grp loop pack';
pg.pulseind(2,:)=1:npulses;
pg.pulseind(1,:)=1:npulses;
pg.chan=[chan offChan];
pg.matrix=eye(8);
phs = 0.9;
pgOn.name=sprintf('RabiToy%02d_%s',plen,name);

iVal = cos(phs);
qVal = sin(phs);
%          PlsLen, wait time, exch, I, Q,  marker pre, marker delay, time
pgOn.params=[plen, 0, qdata(on).eps, iVal, qVal, 10,5,0];
pg1=pgOn.name;
plsdefgrp(pgOn);

pg.pulses.groups={pg1,pg2};
pg.name='IQToyLR';
plsdefgrp(pg);

awgadd(pg.name);
awgcntrl('on start wait err raw');
fprintf('ready!\n');
%% Ramsey, vary frequency 
scan=fConfSeq(pg.name,struct('nloop',128,'nrep',32,'opts','swfb','datachan',tuneData.dataChan));
scan.saveloop = [1 35];
side = upper(tuneData.activeSetName(1));

scan.loops(1).setchan = {'RFfreq3'}; 
scan.loops(1).npoints = 200; 
scan.loops(1).rng = [1.4e9, 1.5e9]; 

scan.loops(2).setchan = 'count'; 
scan.loops(2).npoints = 5; % Improves noise, try to do at least 2. 
smrun(scan,smnext('RamseyIQFreqL'))  
%% Single Echo group, L on, R on
on = 2; off = 1;

clear pg; clear pgOff; clear pgOn; 
npulses=64;
plen=6;
name = qdata(on).name; 
tau = 0.15; 
dt = (-33:30)';
phs = 0.63;

% Off side
offName = qdata(off).name; 
pgOff.pulses=131;
pgOff.ctrl='loop pack';

pgOff.chan = qdata(off).chan;
pgOff.name=sprintf('EchoToy%02d_%d_%s',plen,npulses,offName);
pgOff.dict={struct('prep','@dbzprep','read','@dbzprep','measLoc','@wait','evo','@IQburst'),offName};
dict = pdload(name); 
pgOff.params=[plen, qdata(off).eps, 1, 0, 10,5,tau,0,0];
waitTime = (npulses:-1:1)+dict.meas.time(1)*1e3-49; 
pgOff.varpar = [waitTime', dt]; 
plsdefgrp(pgOff);
pg2=pgOff.name;

pgOn = pgOff;
pgOn.dict={struct('prep','@dbzprep','read','@dbzread','measLoc','@meas','evo','@IQburst'),name};
pgOn.chan = qdata(on).chan;

% Put them together.
pg.ctrl='grp loop pack';
pg.pulseind(2,:)=1:npulses;
pg.pulseind(1,:)=1:npulses;
pg.chan=[pgOn.chan, pgOff.chan];
pg.matrix=eye(8);

pgOn.name=sprintf('EchoToy%02d_%s',plen,name);
waitTime = npulses:-1:1;
pgOn.varpar = [waitTime', dt]; 

iVal = cos(phs);
qVal = sin(phs);
%          PlsLen, exch,           I, Q,  marker pre, marker delay, time
pgOn.params=[plen, qdata(on).eps, iVal, qVal, 10,5,tau,0,0];
pg1=pgOn.name;
plsdefgrp(pgOn);

pg.pulses.groups={pg1,pg2};
pg.name='IQEchoLR';
plsdefgrp(pg);
awgrm(fbdata.lastInd,'after'); 
awgadd(pg.name);
awgcntrl('on start wait err raw');
fprintf('ready!\n');

%% Single Echo group, L off, R on
on = 2; off = 1;

clear pg; clear pgOff; clear pgOn; 
npulses=64;
plen=6;
name = qdata(on).name; 
tau = 0.15; 
dt = (-33:30)';
phs = 0.63;

% Off side
offName = qdata(off).name; 
pgOff.pulses=131;
pgOff.ctrl='loop pack';

pgOff.chan = qdata(off).chan;
pgOff.name=sprintf('EchoToy%02d_%d_%s',plen,npulses,offName);
pgOff.dict={struct('prep','@dbzprep','read','@dbzprep','measLoc','@wait','evo','@IQburst'),offName};
dict = pdload(name); 
pgOff.params=[plen, qdata(off).eps, 0, 0, 10,5,tau,0,0];
waitTime = (npulses:-1:1)+dict.meas.time(1)*1e3-49; 
pgOff.varpar = [waitTime', dt]; 
plsdefgrp(pgOff);
pg2=pgOff.name;

pgOn = pgOff;
pgOn.dict={struct('prep','@dbzprep','read','@dbzread','measLoc','@meas','evo','@IQburst'),name};
pgOn.chan = qdata(on).chan;

% Put them together.
pg.ctrl='grp loop pack';
pg.pulseind(2,:)=1:npulses;
pg.pulseind(1,:)=1:npulses;
pg.chan=[pgOn.chan, pgOff.chan];
pg.matrix=eye(8);

pgOn.name=sprintf('EchoToy%02d_%s',plen,name);
waitTime = npulses:-1:1;
pgOn.varpar = [waitTime', dt]; 

iVal = cos(phs);
qVal = sin(phs);
%          PlsLen, exch,           I, Q,  marker pre, marker delay, time
pgOn.params=[plen, qdata(on).eps, iVal, qVal, 10,5,tau,0,0];
pg1=pgOn.name;
plsdefgrp(pgOn);

pg.pulses.groups={pg1,pg2};
pg.name='IQEchoLR';
plsdefgrp(pg);
awgrm(fbdata.lastInd,'after'); 
awgadd(pg.name);
awgcntrl('on start wait err raw');
fprintf('ready!\n');
%% Single Echo group, L on, R off
on = 2; off = 1;

clear pg; clear pgOff; clear pgOn; 
npulses=64;
plen=6;
name = qdata(on).name; 
tau = 0.15; 
dt = (-33:30)';
phs = 0.63;

% Off side
offName = qdata(off).name; 
pgOff.pulses=131;
pgOff.ctrl='loop pack';

pgOff.chan = qdata(off).chan;
pgOff.name=sprintf('EchoToy%02d_%d_%s',plen,npulses,offName);
pgOff.dict={struct('prep','@dbzprep','read','@dbzprep','measLoc','@wait','evo','@IQburst'),offName};
dict = pdload(name); 
pgOff.params=[plen, qdata(off).eps, 1, 0, 10,5,tau,0,0];
waitTime = (npulses:-1:1)+dict.meas.time(1)*1e3-49; 
pgOff.varpar = [waitTime', dt]; 
plsdefgrp(pgOff);
pg2=pgOff.name;

pgOn = pgOff;
pgOn.dict={struct('prep','@dbzprep','read','@dbzread','measLoc','@meas','evo','@IQburst'),name};
pgOn.chan = qdata(on).chan;

% Put them together.
pg.ctrl='grp loop pack';
pg.pulseind(2,:)=1:npulses;
pg.pulseind(1,:)=1:npulses;
pg.chan=[pgOn.chan, pgOff.chan];
pg.matrix=eye(8);

pgOn.name=sprintf('EchoToy%02d_%s',plen,name);
waitTime = npulses:-1:1;
pgOn.varpar = [waitTime', dt]; 

iVal = cos(phs);
qVal = sin(phs);
%          PlsLen, exch,           I, Q,  marker pre, marker delay, time
pgOn.params=[plen, qdata(on).eps, 0, 0, 10,5,tau,0,0];
pg1=pgOn.name;
plsdefgrp(pgOn);

pg.pulses.groups={pg1,pg2};
pg.name='IQEchoLR';
plsdefgrp(pg);
awgrm(fbdata.lastInd,'after'); 
awgadd(pg.name);
awgcntrl('on start wait err raw');
fprintf('ready!\n');

%% Run scan, single Echo 
smset('RFfreq3',1.45e9); 
scan=fConfSeq(pg.name,struct('nloop',128,'nrep',32,'opts','swfb','datachan',tuneData.dataChan));
smrun(scan,smnext('EchoIQL')); 
sleep
%% Run scan, Echo vary phase
smset('RFfreq3',1.62e9); 
scan=fConfSeq(grps,struct('nloop',128,'nrep',32,'opts','swfb','datachan',tuneData.dataChan));
smrun(scan,smnext('EchoPhaseR')); 
sleep
%% Run scan, Echo vary frequency
smset('RFpow3',26); 
scan=fConfSeq(pg.name,struct('nloop',64,'nrep',32,'opts','swfb','datachan','DAQ2'));
scan.saveloop = [1 35];
side = upper(tuneData.activeSetName(1));

scan.loops(1).setchan = {'RFfreq3'}; 
scan.loops(1).npoints = 300; 
scan.loops(1).rng = [1.5e9, 0.5e9]; 

scan.loops(2).setchan = 'count'; 
scan.loops(2).npoints = 14; % Improves noise, try to do at least 2. 
smrun(scan,smnext('EchoIQFreqR'))  
sleep
%% Echo, vary phase
on = 2; off = 1;

clear pg; clear pgOff; clear pgOn; 
npulses=64;
plen=6;
name = qdata(on).name; 
tau = 0.15; 
dt = (-33:30)';
grps = {};

% Off side
offName = qdata(off).name; 
pgOff.pulses=131;
pgOff.ctrl='loop pack';

pgOff.chan = qdata(off).chan;
pgOff.name=sprintf('IQEcho_%s',offSide);
pgOff.dict={struct('prep','@dbzprep','read','@dbzprep','measLoc','@wait','evo','@IQburst'),offName};
dict = pdload(name); 
pgOff.params=[plen, qdata(off).eps, 1, 0, 10,5,tau,0,0];
waitTime = (npulses:-1:1)+dict.meas.time(1)*1e3-49; 
pgOff.varpar = [waitTime', dt]; 
plsdefgrp(pgOff);
pg2=pgOff.name;

pgOn = pgOff;
pgOn.dict={struct('prep','@dbzprep','read','@dbzread','measLoc','@meas','evo','@IQburst'),name};
pgOn.chan = qdata(on).chan;

% Put them together.
pg.ctrl='grp loop pack';
pg.pulseind(2,:)=1:npulses;
pg.pulseind(1,:)=1:npulses;
pg.chan=[pgOn.chan, pgOff.chan];
pg.matrix=eye(8);

pgOn.name=sprintf('IQEcho_%s',side);
waitTime = npulses:-1:1;
pgOn.varpar = [waitTime', dt]; 

for i = 1:ngrps
    pgOn.name=sprintf('EchoPhase%s_%d',side,i);
    
    iVal = cos(phaseVals(i));
    qVal = sin(phaseVals(i));    
    pgOn.params=[plen, qdata(on).eps, iVal, qVal, 10,5,tau,0,0];    
    pg1=pgOn.name;
    plsdefgrp(pgOn);
            
    pg.pulses.groups={pg1,pg2};
    pg.name=sprintf('IQToyLR_%02d',i);
    grps{end+1}=pg.name;      
    plsdefgrp(pg);
end

awgrm(fbdata.lastInd,'after'); 
awgadd(grps);
awgcntrl('on start wait err raw');
fprintf('ready!\n');
%% Echo, 4 groups, L on, R on, both on, both off. 
on = 2; off = 1;

clear pg; clear pgOff; clear pgOn; 
npulses=64;
plen=6;
name = qdata(on).name; 
tau = 0.15; 
dt = (-33:30)';
grps = {};

% Off side
offName = qdata(off).name; 
pgOff.pulses=131;
pgOff.ctrl='loop pack';

pgOff.chan = qdata(off).chan;

pgOff.dict={struct('prep','@dbzprep','read','@dbzprep','measLoc','@wait','evo','@IQburst'),offName};
dict = pdload(name); 

waitTime = (npulses:-1:1)+dict.meas.time(1)*1e3-49; 
pgOff.varpar = [waitTime', dt]; 

pgOn = pgOff;
pgOn.dict={struct('prep','@dbzprep','read','@dbzread','measLoc','@meas','evo','@IQburst'),name};
pgOn.chan = qdata(on).chan;

% Put them together.
pg.ctrl='grp loop pack';
pg.pulseind(2,:)=1:npulses;
pg.pulseind(1,:)=1:npulses;
pg.chan=[pgOn.chan, pgOff.chan];
pg.matrix=eye(8);

pgOn.name=sprintf('IQEcho_%s',side);
waitTime = npulses:-1:1;
pgOn.varpar = [waitTime', dt]; 

iVal = cos(phs);
qVal = sin(phs);    
ampOff = [0,1,0,1]; 
ampOn = [0,0,1,1]; 
for i = 1:length(ampOn)
    pgOff.params=[plen, qdata(off).eps, ampOff(i), 0, 10,5,tau,0,0];
    pgOff.name=sprintf('IQEcho%s_%d',offSide,i);
    plsdefgrp(pgOff);
    pgOn.name=sprintf('IQEcho%s_%d',side,i);
    pgOn.params=[plen, qdata(on).eps, ampOn(i)*iVal, ampOn(i)*qVal, 10,5,tau,0,0];
    plsdefgrp(pgOn);
    
    pg.pulses.groups={pgOn.name,pgOff.name};
    pg.name=sprintf('EchoPhaseLR_%02d',i);
    grps{end+1}=pg.name;
    plsdefgrp(pg);
end
awgrm(fbdata.lastInd,'after'); 
awgadd(grps);
awgcntrl('on start wait err raw');
fprintf('ready!\n');
%% Run scan, 4 groups, L on, R on, both on, both off. 
smset('RFfreq3',1.62e9); 
scan=fConfSeq(grps,struct('nloop',128,'nrep',64,'opts','swfb','datachan',tuneData.dataChan));
smrun(scan,smnext('EchoOnOff_R')); 
%%
% 1.413 seems to be best freq. 