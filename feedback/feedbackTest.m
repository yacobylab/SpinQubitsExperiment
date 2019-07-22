function feedbackTest(ctrl)
% update feedback pulses based on current dictionary. To make bigger
% changes (trafofn, meas time, pulse spacing, etc.) must edit function directly. 
% function feedbackTest(ctrl)
% ctrl can be all, pump, dbz, dbz2.
%   all: updates all pulses
%   pump the feedback pulses
%   dbz/dbz2 the measuring dbz pulses (dbz2 has larger spacing for faster
%   measurement)

global tuneData; global fbdata; global smdata; 
side = tuneData.activeSetName;
switch side
    case 'left'
        ind = 1;
    case 'right'
        ind = 2;
end
dict = pdload(side);
pg.chan=[getNum(tuneData.xyChan{1}),getNum(tuneData.xyChan{2})];
pg.dict=side;
samprate = smdata.inst(inl('ATS')).data.defSamprate; 
pulseInc = 1/samprate*1e9; 

switch ctrl
    case 'all'
        feedbackTest('pump');
        feedbackTest('dbz');
        feedbackTest('dbz2');
    case 'pump'
        pg.name = sprintf('pump_%s',upper(side(1)));
        pg.ctrl='loop';
        %pg.trafofn.func=@rc_trafofn;  pg.trafofn.args=-1;
        pg.varpar = 0; % Fully hardcoded?
        pg.params = 0;
        pg.pulses = [114,113];  % tpol, then spol.
        %   singlet goes from S to T
        pg.nrep = [0,0];
        pg.pulseind = [1 2];
        pg.jump = [];
        plsdefgrp(pg);
        awgadd(pg.name);
        awgcntrl('start on err');
        seq = awgseqind(pg.name);
        fbdata.params(ind).singletPulse=seq+2;
        fbdata.params(ind).tripletPulse=seq+1;
    case 'dbz' % dBz group for software feedback
        pg.pulses=12;
        pulses=128; % number of pulses
        scale=1; % time spacing betwen pulses 
        waitTime = ceil(pulses*scale/pulseInc)*pulseInc; 
        measTime = dict.meas.time(1); 
        minTime=measTime+sum(dict.reload.time)+pulses*scale/1000; % meas time, can reduce from dict val to speed up. 
        plsLen = ceil(minTime*4)/4+0.25; % Round to closest 250 ns. 
        namepat=sprintf('dBz_swfb_%d_%s',pulses,upper(side(1)));
        pg.ctrl='loop pack';
        pg.trafofn.func=@rc_trafofn; pg.trafofn.args=-1;
        pg.params=[plsLen, waitTime, measTime, 0]; %pulselength, max dt, meas time, septime
        %pg.params=[pulses*scale+1, mtime, 0]; %pulselength, max dt, meas time, septime
        pg.varpar=(1:pulses)'*scale;
        pg.name=namepat;
        plsdefgrp(pg);
        awgadd(pg.name);
        fbdata.params(ind).fbInit=namepat;
    case 'dbz2'
        pg.pulses=12;
        pulses=128; % number of pulses
        scale=4;                
        measTime = dict.meas.time(1); 
        minTime=measTime+sum(dict.reload.time)+pulses*scale/1000; % meas time, can reduce from dict val to speed up. 
        plsLen = ceil(minTime*4)/4+0.25; % Round to closest 250 ns. 
                
        namepat=sprintf('dBz_swfb_%d_%s',pulses*scale,upper(side(1)));
        pg.ctrl='loop pack';
        %pg.trafofn.func=@rc_trafofn; pg.trafofn.args=.5;
        pg.params=[plsLen, pulses*scale+1, measTime, 0]; %pulselength, max dt, meas time, septime
        pg.varpar=((0:(pulses-1))'*scale);
        pg.name=namepat;
        plsdefgrp(pg);
        awgadd(pg.name);
        %fbdata.params(ind).fbGroup=namepat;
end
end