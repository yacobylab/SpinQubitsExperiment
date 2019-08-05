name = 'right'; side = name(1); 
ind = strcmp({qdata.name},name); 
offQub = ~ind;
offName = qdata(offQub).name; offSide = offName(1); 

clear pgOff; clear pgOn; 
npulses=64;
plen=8;
tau = 0.15; 
dt = (-30:35)';
phs = 0.63;
customDict = struct('prep','@dbzprep','read','@dbzprep','measLoc','@meas','evo','@IQburst');  

pgOff.pulses=145;
pgOff.ctrl='loop pack';

% Off side
pgOff.chan = qdata(off).chan;
pgOff.name=sprintf('condEvoIQ_%s',upper(offSide(1)));
stagDict = sprintf('stag%s',offSide);
pgOff.dict={customDict,stagDict,offName};

dict = pdload(name); 
pgOff.params=[plen, qdata(off).eps, 1, 0, 10,5,tau,dbzTime,0,0];
waitTime = (npulses:-1:1)+dict.meas.time(1)*1e3-49; 
pgOff.varpar = [waitTime', dt]; 
plsdefgrp(pgOff);

pgOn = pgOff;
pgOn.chan = qdata(on).chan;
stagDict = sprintf('stag%s',side);
pgOn.dict={customDict,stagDict,name};
pgOn.name=sprintf('condEvoIQ_%s',upper(side(1)));
waitTime = npulses:-1:1;
pgOn.varpar = [waitTime', dt]; 

iVal = cos(phs); qVal = sin(phs);
%          PlsLen, exch,         I, Q,  mk pre/delay, time, wait
pgOn.params=[plen, qdata(on).eps, iVal, qVal, 10,5,tau,0,0];
plsdefgrp(pgOn);

grpName='condEvoIQ_LR';
make2grp(pgOn,pgOff,grpName); 
awgrm(fbdata.lastInd,'after'); awgclear('unused'); 
awgadd(grpName);
awgcntrl('on start wait err raw');
%%
grps = {};
name = 'right'; 
nGrps = 17; 
clear pg; 
plen=8;
tau = 0.15; 
dt = (-30:35)';
npulses = length(dt); 
phs = 0.63;
customDict = struct('prep','@dbzprep','read','@dbzprep','measLoc','@meas','evo','@IQburst');
for j=1:nGrps
    for i = 1:length(qdata)
        currName = qdata(i).name;
        
        pg(i).ctrl='loop pack';
        pg(i).chan = qdata(i).chan;
        pg(i).name=sprintf('condEvoIQ_%02d_%s',j,upper(currName(1)));
        stagDict = sprintf('stag%s',lower(currName(1)));
        pg(i).dict={customDict,stagDict,currName};
        if ~strcmp(currName,name)
            pg(i).pulses=145;
            dbzTime = j-1;             
            pg(i).params=[plen, qdata(off).eps, 1, 0, 10,5,tau,dbzTime,0,0];
            waitTime = (npulses:-1:1)+1e3-49;            
        else
            pg(i).pulses = 146;
            waitTime = (npulses:-1:1)+dbzTime; % dbz time.             
            %          PlsLen, exch,         I, Q,  mk pre/delay, time, wait
            iVal = cos(phs); qVal = sin(phs);
            pg(i).params=[plen, qdata(on).eps, iVal, qVal, 10,5,tau,0,0];
        end
        pg(i).varpar = [waitTime', dt];
        plsdefgrp(pg(i));                
    end
    grps{j}=sprintf('condEvoIQ_%02d_LR',j);
    make2grp(pg(1),pg(2),grps{j});
end
awgrm(fbdata.lastInd,'after'); awgclear('unused');
awgadd(grps);
awgcntrl('on start wait err raw');
%%
smset('RFfreq3',1.45e9); 
scan=fConfSeq(grps,struct('nloop',128,'nrep',32,'opts','swfb'));
smrun(scan,smnext('condEvo_LR')); 
sleep