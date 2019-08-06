%% Make groups for condEvo. 
name = 'right'; 
nGrps = 17; 
plen=8;
tau = 0.15; 
dt = (-30:55)';
phs = 0.63;
namePat = 'condEvoIQ_%02d_%s';

clear pg; 
npulses = length(dt); 
grps = {};
for i=1:nGrps
    for j = 1:length(qdata)
        currName = qdata(j).name;
        
        pg(j).ctrl='loop pack';
        pg(j).chan = qdata(j).chan;
        pg(j).name=sprintf(namepat,i,upper(currName(1)));
        stagDict = sprintf('stag%s',lower(currName(1)));      
        if ~strcmp(currName,name)
            pg(j).pulses=145;
            nullPrep = struct('type','wait','time',8e-3,'val',[0,0]); 
            nullPi = struct('type','wait','time',16e-3,'val',[0,0]); 
            customDict = struct('prep','@dbzprep','read',nullPrep,'measLoc','@meas','evo','@IQburst','pi',nullPi);
            dbzTime = i-1;             
            pg(j).params=[plen, qdata(off).eps, 1, 0, 10,5,tau,dbzTime,0,0];
            waitTime = (npulses:-1:1)+1e3-49;            
        else
            pg(j).pulses = 146;
            waitTime = (npulses:-1:1)+dbzTime; % dbz time.             
            customDict = struct('prep','@dbzprep','read','@dbzprep','measLoc','@meas','evo','@IQburst');
            %          PlsLen, exch,         I, Q,  mk pre/delay, time, wait
            iVal = cos(phs); qVal = sin(phs);
            pg(j).params=[plen, qdata(on).eps, iVal, qVal, 10,5,tau,0,0];
        end
        pg(j).dict={customDict,stagDict,currName};
        pg(j).varpar = [waitTime', dt];
        plsdefgrp(pg(j));                
    end
    grps{i}=sprintf(namepat,i,'LR');
    make2grp(pg(1),pg(2),grps{i});
end
awgrm(fbdata.lastInd,'after'); awgclear('unused');
awgadd(grps);
awgcntrl('on start wait err raw');
%%
smset('RFfreq3',1.45e9); 
scan=fConfSeq(grps,struct('nloop',128,'nrep',32,'opts','swfb'));
smrun(scan,smnext('condEvo_LR')); 
sleep
%% Make groups for condEvo. 
name = 'right'; 
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
        else
            pg(j).pulses = 146;
            waitTime = (npulses:-1:1)+dbzTime; % dbz time.             
            customDict = struct('prep','@dbzprep','read','@dbzprep','measLoc','@meas','evo','@IQburst');
            %          PlsLen, exch,         I, Q,  mk pre/delay, time, wait
            iVal = cos(phs); qVal = sin(phs);
            pg(j).params=[plen, qdata(j).eps, iVal, qVal, 10,5,tau,0,0];
        end
        pg(j).varpar = [waitTime', dt];
        pgFin(j) = make1grp(pg(j),{'n',i,'side',qdata(j).name,'namePat',namePat,'customDict',customDict,'opts','stag two'});
    end
    grps{i}=sprintf(namePat,i,'LR');
    make2grp(pgFin(1),pgFin(2),grps{i});
end
awgrm(fbdata.lastInd,'after'); awgclear('unused');
awgadd(grps);
awgcntrl('on start wait err raw');
%% Echo on both sides
name = 'right'; 
nGrps = 1; 
plen = 7;
tau = 0.15; 
dt = (-30:55)';
phs = 0.63;
namePat = 'EchoIQ_%02d_%s';

npulses = length(dt); 
grps = {}; pg = [];
customDict = struct('prep','@dbzprep','read','@dbzprep','measLoc','@meas','evo','@IQburst');
for i=1:nGrps
    for j = 1:length(qdata)              
        pg(j).pulses = 146;        
        waitTime = npulses:-1:1;            
        if ~strcmp(qdata(j).name,name) % Off side
            pg(j).params=[plen, qdata(j).eps, 1, 0, 10,5,tau,0,0];            
        else
            %waitTime = (npulses:-1:1)+dbzTime; % dbz time.                         
            %          PlsLen, exch,         I, Q,  mk pre/delay, time, wait
            iVal = cos(phs); qVal = sin(phs);
            pg(j).params=[plen, qdata(j).eps, iVal, qVal, 10,5,tau,0,0];
        end
        pg(j).varpar = [waitTime', dt];
        pgFin(j) = make1grp(pg(j),{'n',i,'side',qdata(j).name,'namePat',namePat,'customDict',customDict,'opts','stag two'});
    end
    grps{i}=sprintf(namePat,i,'LR');
    make2grp(pgFin(1),pgFin(2),grps{i});
end
awgrm(fbdata.lastInd,'after'); awgclear('unused');
awgadd(grps);
awgcntrl('on start wait err raw');
%%
scan=fConfSeq(grps,struct('nloop',128,'nrep',32,'opts','swfb'));
scan.saveloop = [1 35];
side = upper(tuneData.activeSetName(1));

scan.loops(1).setchan = {'RFfreq3'}; 
scan.loops(1).npoints = 300; 
scan.loops(1).rng = [2e9, 1e9]; 

scan.loops(2).setchan = 'count'; 
scan.loops(2).npoints = 14; % Improves noise, try to do at least 2. 

smrun(scan,smnext('EchoIQ_LR')); 
sleep
%% Echo on both sides
name = 'right'; 
nGrps = 1; 
plen = 7;
tau = 0.15; 
dt = (-30:55)';
phs = 0.63;
namePat = 'EchoIQ_%02d_%s';

npulses = length(dt); 
grps = {}; pg = [];
customDict = struct('prep','@dbzprep','read','@dbzprep','measLoc','@meas','evo','@IQburst');
for i=1:nGrps
    for j = 1:length(qdata)              
        pg(j).pulses = 146;        
        waitTime = npulses:-1:1;            
        if ~strcmp(qdata(j).name,name) % Off side
            pg(j).params=[plen, qdata(j).eps, 1, 0, 10,5,tau,0,0];            
        else
            %waitTime = (npulses:-1:1)+dbzTime; % dbz time.                         
            %          PlsLen, exch,         I, Q,  mk pre/delay, time, wait
            iVal = cos(phs); qVal = sin(phs);
            pg(j).params=[plen, qdata(j).eps, iVal, qVal, 10,5,tau,0,0];
        end
        pg(j).varpar = [waitTime', dt];
        pgFin(j) = make1grp(pg(j),{'n',i,'side',qdata(j).name,'namePat',namePat,'customDict',customDict,'opts','stag two'});
    end
    grps{i}=sprintf(namePat,i,'LR');
    make2grp(pgFin(1),pgFin(2),grps{i});
end
awgrm(fbdata.lastInd,'after'); awgclear('unused');
awgadd(grps);
awgcntrl('on start wait err raw');
