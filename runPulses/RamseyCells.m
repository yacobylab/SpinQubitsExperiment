%% Single scan 
scan.loops(1).rng = [1.35e9 1.6e9]; 
smset('RFpow3',pow(i));
smrun(scan,smnext(sprintf('RamseySwitch%s',side)));
%% Sweep power and frequency (2d scan)
scan = fConfSeq('Ramsey_L_1',{'nloop',75,'nrep',40, 'datachan',tuneData.dataChan,'opts','swfb'}); % with feedback
scan.saveloop = [1 35];
side = upper(tuneData.activeSetName(1));

scan.loops(1).setchan = {'RFfreq3'}; 
scan.loops(1).npoints = 500; 

scan.loops(2).setchan = 'count'; 
scan.loops(2).npoints = 4; % Improves noise, try to do at least 2. 

pow = [-5,0];
%pow = [-23]; 
for j = 4:-1:4
    for i = 1:length(pow)
        freqMin = -5e8+j*1e9; 
        scan.loops(1).rng = [freqMin,freqMin+1e9];         
       % scan.loops(1).rng = [1.4e9 1.5e9]; 
        smset('RFpow3',pow(i));        
        smrun(scan,smnext(sprintf('RamseySwitch%s',side)));
        sleep;
        pause(i*20)
    end
end
%% Sweep power within a scan for the resonance frequency. 
scan = fConfSeq('Ramsey_L_1',{'nloop',75,'nrep',40, 'datachan',tuneData.dataChan,'opts','swfb'}); % with feedback
scan.loops(1).setchan = {'RFpow3'}; 
scan.loops(1).rng = [-35,-25]; 
scan.loops(1).npoints = range(scan.loops(1).rng)*10+1; 
scan.saveloop = [1 35];
side = upper(tuneData.activeSetName(1));
scan.loops(2).setchan = 'count'; 
scan.loops(2).npoints = 40; 
smset('RFfreq3',1.127e9);
smrun(scan,smnext(sprintf('RamseyPow%s',side)));
sleep;
%% Make the group with the switch
clear pg 
pg.pulses=118;
pg.ctrl='loop pack'; 
pg.chan=[2 1]; pg.dict='left';
l=pdload('left');
pg.trafofn.func=@skinTraf; pg.trafofn.args=10;
pg.dict={struct('prep',struct('type','@dbzprep'),'read',struct('type','@dbzread')),pg.dict};
eps = 0.9;
evo = 1:50;
pg.varpar = evo';

Ramsey = {};
namepat= sprintf('Ramsey_L_%d',1);
%         pulse length, exch dir (2), exch eps,
exchVals = [l.exch.val(1), l.exch.val(2), eps];
%                  ring up time, end delay time,
%                            wait time (coupled with ringup) , meas time, varpar.
pg.params=[5, exchVals, 200, 20, 0.2, 0.5 , 2.3, 0];
pg.name = namepat;
plsdefgrp(pg);
%Ramsey{i} = pg.name;
%end

awgadd(pg.name);
awgcntrl('on start wait err raw');

%%
clear pg 
pg.pulses=131;
pg.ctrl='loop pack'; 
pg.chan=[2 1]; pg.dict='left';
l=pdload('left');
pg.trafofn.func=@skinTraf; pg.trafofn.args=10;
pg.dict={struct('prep',struct('type','@dbzprep'),'read',struct('type','@dbzread')),pg.dict};
eps = 0.9;
evo = 1:50;
pg.varpar = evo';

Ramsey = {};
namepat= sprintf('RamseyE_L_%d',1);
%         pulse length, exch eps, echo time
%                       ring up time, end delay time,
%                              wait time (coupled with ringup), varpar.
pg.params=[5, eps, 200, 10, 20, 0.3, 0.01, 0];
pg.name = namepat;

plsdefgrp(pg);
%Ramsey{i} = pg.name;

awgadd(pg.name);
awgcntrl('on start wait err raw');
%%
clear pg 
pg.pulses=132;
pg.ctrl='loop pack'; 
pg.chan=[2 1,5,6]; pg.dict='left';
%pg.trafofn.func=@skinTraf; pg.trafofn.args=10;
pg.dict={struct('prep',struct('type','@dbzprep'),'read',struct('type','@dbzread')),pg.dict};
eps = 0.9;
evo = 100*ones(1,100);
pg.varpar = evo';

Ramsey = {};
namepat= 'IQtoy';
%     PlsLen, exch, I,  Q,  marker pre, marker delay, time
pg.params=[5, eps, 0, 1, 20, 20, 0];
pg.name = namepat;

plsdefgrp(pg);
%Ramsey{i} = pg.name;

awgadd(pg.name);
awgcntrl('on start wait err raw');
%% Run the group by itself 
scan=fConfSeq('Ramsey_L_1',struct('nloop',75,'nrep',100,'datachan','DAQ1','opts','swfb'));
smrun(scan,smnext('RamseyL')); 
sleep