function feedbackTest(ctrl)
% update feedback pulses. 
%function out=feedbackTest(ctrl)
% ctrl can be all, pump, dbz, dbz2. 
% all updates all pulses, pump the feedback pulses, dbz/dbz2 the measuring
% dbz pulses. 

global tuneData; global fbdata;
side = tuneData.activeSetName;
switch side 
    case 'left' 
        ind = 1;         
    case 'right' 
        ind = 2;         
end
dict = pdload(side); 
pg.chan=[str2double(char(regexp(tuneData.xyChan{1},'\d+','match'))),str2double(char(regexp(tuneData.xyChan{2},'\d+','match')))];
pg.dict=side;
switch ctrl
    case 'all'
        feedbackTest('pump'); 
        feedbackTest('dbz'); 
        feedbackTest('dbz2'); 
    case 'pump'
        pg.name = sprintf('pump_%s',upper(side(1))); %will get changed to 'R' below
        pg.ctrl='loop';             
        %pg.trafofn.func=@rc_trafofn;         pg.trafofn.args=-1;
        pg.varpar = 0; 
        pg.params = 0;
        pg.pulses = [114,113];  % tpol, then spol. 
        %   singlet goes from S to T 
        pg.nrep = [0,0];
        pg.pulseind = [1 2];
        pg.jump = [];
        plsdefgrp(pg); 
        awgadd(pg.name); 
        awgcntrl('start on err'); 
        seq= awgseqind(pg.name); 
        fbdata.params(ind).singletPulse=seq+2; 
        fbdata.params(ind).tripletPulse=seq+1;     
    case 'dbz' % dBz group for software feedback
        pg.pulses=12;
        pulses=128; %number of pulses
        scale=1;
                
        mtime=dict.meas.time(1); %meas time
        pulseLen=1+ceil(mtime);
        namepat=sprintf('dBz_swfb_%d_%s',pulses,upper(side(1)));
        pg.ctrl='loop pack';
        pg.trafofn.func=@rc_trafofn; pg.trafofn.args=-.5;
        pg.params=[pulseLen, pulses*scale+1, mtime, 0]; %pulselength, max dt, meas time, septime
        %pg.params=[pulses*scale+1, mtime, 0]; %pulselength, max dt, meas time, septime
        pg.varpar=((0:(pulses-1))'*scale);
        pg.name=namepat;
        plsdefgrp(pg);
        awgadd(pg.name);
        fbdata.params(ind).fbInit=namepat;                
    case 'dbz2' 
        pg.pulses=12;
        pulses=128; %number of pulses
        scale=4;
        mtime=dict.meas.time(1); %meas time
        pulseLen=1.5+ceil(mtime);        
        namepat=sprintf('dBz_swfb_%d_%s',pulses*scale,upper(side(1)));
        pg.ctrl='loop pack';
        pg.trafofn.func=@rc_trafofn; pg.trafofn.args=-.5;
        pg.params=[pulseLen, pulses*scale+1, mtime, 0]; %pulselength, max dt, meas time, septime
        pg.varpar=((0:(pulses-1))'*scale);
        pg.name=namepat;
        plsdefgrp(pg);
        awgadd(pg.name);
        %fbdata.params(ind).fbGroup=namepat;          
end
end