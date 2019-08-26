%% Make groups for condEvo. 
name = 'left'; 
nGrps = 17; 
plen=8;
tau = 0.15; 
dt = (-30:55)';
phs = 0.63;
namePat = 'condEvoIQ_%02d_%s';

npulses = length(dt);
grps = {}; pg = [];
for i=1:nGrps
    for j = 1:length(qdata)
        dbzTime = i-1;
        if ~strcmp(qdata(j).name,name)
            pg(j).pulses=145;
            nullPrep = struct('type','wait','time',8e-3,'val',[0,0]);
            nullPi = struct('type','wait','time',16e-3,'val',[0,0]);
            customDict = struct('prep','@dbzprep','read',nullPrep,'measLoc','@meas','evo','@IQburst','pi',nullPi);
            pg(j).params=[plen, qdata(j).eps, 1, 0, 10,5,tau,dbzTime,0,0];
            waitTime = (npulses:-1:1)+1e3-49;
            pg(j).varpar = [waitTime', dt];
            pgFin(j) = make1grp(pg(j),{'n',i,'side',qdata(j).name,'namePat',namePat,'customDict',customDict,'opts','stag two'});
            else
            pg(j).pulses = 146;
            waitTime = (npulses:-1:1)+dbzTime; % dbz time.
            customDict = struct('prep','@dbzprep','read','@dbzprep','measLoc','@meas','evo','@exch');
            %          PlsLen, exch,         I, Q,  mk pre/delay, time, wait
            iVal = cos(phs); qVal = sin(phs);
            pg(j).params=[plen, qdata(j).eps, iVal, qVal, 10,5,tau,0,0];
            pg(j).varpar = [waitTime', dt];
            pgFin(j) = make1grp(pg(j),{'n',i,'side',qdata(j).name,'namePat',namePat,'customDict',customDict,'opts','stag'});
        end
    end
    grps{i}=sprintf(namePat,i,'LR');
    make2grp(pgFin(1),pgFin(2),grps{i});
end
awgrm(fbdata.lastInd,'after'); awgclear('unused');
awgadd(grps);
awgcntrl('on start wait err raw');
%% Run scan for condEvo
%for i = -0.014e9:0.007e9:0.014e9
smset('RFfreq3',0.8842e9); 
%smset('RFfreq3',0.8982e9); 
smset('RFpow3',-10);
scan=fConfSeq(grps,struct('nloop',64,'nrep',128,'opts','swfb','datachan','DAQ1'));
smrun(scan,smnext('EchoIQ_L')); 
sleep
%end
%% Run scan for condEvo at various power and frequencies
for j = 0.85e9:0.02e9:0.89e9
smset('RFfreq3',j);
for i=14:4:22
smset('RFpow3',i);
scan=fConfSeq(grps,struct('nloop',64,'nrep',64,'opts','swfb'));
smrun(scan,smnext('condEvo_LR')); 
sleep
end
end
for j = 1.05e9:0.05e9:1.15e9
smset('RFfreq3',j); 
for i=14:4:22
smset('RFpow3',i);
scan=fConfSeq(grps,struct('nloop',64,'nrep',64,'opts','swfb'));
smrun(scan,smnext('condEvo_LR')); 
sleep
end
end
for j = 1.4e9:0.05e9:1.5e9
smset('RFfreq3',j); 
for i=14:4:22
smset('RFpow3',i);
scan=fConfSeq(grps,struct('nloop',64,'nrep',64,'opts','swfb'));
smrun(scan,smnext('condEvo_LR')); 
sleep
end
end
%for i=10:2:22
%smset('RFfreq3',1.45e9); 
%smset('RFpow3',i);
%scan=fConfSeq(grps,struct('nloop',64,'nrep',128,'opts','swfb'));
%smrun(scan,smnext('condEvo_LR')); 
%sleep
%end

smset('RFpow3',10);
%% Run scan for echo at fixed frequency
%smset('RFfreq3',0.8842e9); 
smset('RFfreq3',0.5e9); 
smset('RFpow3',-35);
scan=fConfSeq(grps,struct('nloop',64,'nrep',64,'opts','swfb','datachan','DAQ1'));
smrun(scan,smnext('EchoIQFreq_L')); 
sleep
%% Make groups for echo on both sides
nGrps = 1; 
plen = 7;
tau = 0.150; 
dt = (-30:55)';
namePat = 'EchoIQ';

npulses = length(dt);
grps = {}; pg = [];
customDict = struct('prep','@dbzprep','read','@dbzprep','measLoc','@meas','evo','@IQburst');
for i=1:nGrps
    for j = 1:length(qdata)
        pg(j).pulses = 146;
        waitTime = npulses:-1:1;
        iVal = cos(qdata(j).phs); qVal = sin(qdata(j).phs);
        %          PlsLen, exch,         I, Q,  mk pre/delay, time, wait
        pg(j).params=[plen, qdata(j).eps, iVal, qVal, 10,5,tau,0,0];
        pg(j).varpar = [waitTime', dt];
        pgFin(j) = make1grp(pg(j),{'n',i,'side',qdata(j).name,'namePat',namePat,'customDict',customDict,'opts','stag two'});
    end
    grps{i}=sprintf([namePat '_%s'],'LR');
    make2grp(pgFin(1),pgFin(2),grps{i});
end
pg1 = pgFin(1).name;
pg2 = pgFin(2).name;
addGroups(grps); 
%% Run double echo scan while varying frequency. 
scan=fConfSeq(grps,struct('nloop',64,'nrep',32,'opts','swfb','datachan','DAQ1'));
scan.saveloop = [1 35];
scan.loops(1).setchan = {'RFfreq3'}; 
scan.loops(1).npoints = 900;

smset('RFpow3',5)
scan.loops(1).rng = [0.5e9, 2e9]; 
scan.loops(2).setchan = 'count';
scan.loops(2).npoints = 14; % Improves noise, try to do at least 2. 
smrun(scan,smnext([namePat '_L']));
sleep
%% Make groups to vary phase, measuring on both sides. 
name = 'right'; 
nGrps = 20; 
plen = 7;
tau = 0.15; 
dt = (-30:55)';
namePat = 'EchoPhaseIQ';

npulses = length(dt);
grps = {}; pg = [];
phaseVals = linspace(0,pi); 
customDict = struct('prep','@dbzprep','read','@dbzprep','measLoc','@meas','evo','@IQburst');
for i=1:nGrps
    for j = 1:length(qdata)
        pg(j).pulses = 146;
        waitTime = npulses:-1:1;
        if ~strcmp(qdata(j).name,name)
            iVal = 1; qVal = 0; 
        else
            iVal = cos(phaseVals(i)); qVal = sin(phaseVals(i)); 
        end
        %          PlsLen, exch,         I, Q,  mk pre/delay, time, wait
        pg(j).params=[plen, qdata(j).eps, iVal, qVal, 10,5,tau,0,0];
        pg(j).varpar = [waitTime', dt];
        pgFin(j) = make1grp(pg(j),{'n',i,'side',qdata(j).name,'namePat',namePat,'customDict',customDict,'opts','stag two'});
    end
    grps{i}=sprintf([namePat '_%02d_%s'],i,'LR');
    make2grp(pgFin(1),pgFin(2),grps{i});
end
pg1 = pgFin(1).name;
pg2 = pgFin(2).name;
addGroups(grps); 
%% Run phase varying scan for a variety of frequencies.
%freqVals = [0.9e9,0.92e9];
smset('RFpow3',5);
%freqVals = [0.8842e9,1.45e9,1.7684e9,1.8e9];
freqVals = 0.8842e9;
for i = 1:length(freqVals)
    smset('RFfreq3',freqVals(i));
    scan=fConfSeq(grps,struct('nloop',128,'nrep',32,'opts','swfb'));
    smrun(scan,smnext([namePat '_LR']));
    pause(30); 
    sleep
end
%% Make groups with RF on opposite qubit. 
name = 'left'; 
nGrps = 1; 
plen = 4;
tau = 0.15; 
dt = (-30:55)';
namePat = 'EchoIQ';
dict = pdload(name); 
npulses = length(dt);
grps = {}; pg = [];

for i=1:nGrps
    for j = 1:length(qdata)
        pg(j).pulses = 146;
        waitTime = npulses:-1:1;
        
        if ~strcmp(qdata(j).name,name)
            nullMeas = struct('type','wait','time',dict.meas.time(1),'val',[2,-2]);
            customDict = struct('prep','@dbzprep','read','@dbzprep','measLoc',nullMeas,'evo','@IQburst');
            optList = 'two'; 
        else
            customDict = struct('prep','@dbzprep','read','@dbzprep','measLoc','@meas','evo','@exch');
            optList = ''; 
        end
        iVal = cos(qdata(j).phs); qVal = sin(qdata(j).phs);
        %          PlsLen, exch,         I, Q,  mk pre/delay, time, wait
        pg(j).params=[plen, qdata(j).eps, iVal, qVal, 10,5,tau,0,0];
        %pg(j).params=[plen, qdata(j).eps, 0, 0, 10,5,tau,0,0]; % if you
        %want to turn RF off on the opposite dot
        pg(j).varpar = [waitTime', dt];
        pgFin(j) = make1grp(pg(j),{'n',i,'side',qdata(j).name,'namePat',namePat,'customDict',customDict,'opts',optList});
    end
    grps{i}=sprintf([namePat '_%s'],'LR');
    make2grp(pgFin(1),pgFin(2),grps{i});
end
pg1 = pgFin(1).name;
pg2 = pgFin(2).name;
addGroups(grps);
%%
scan=fConfSeq(grps,struct('nloop',128,'nrep',32,'opts','swfb','datachan','DAQ2'));
smrun(scan,smnext('EchoIQ_R'));