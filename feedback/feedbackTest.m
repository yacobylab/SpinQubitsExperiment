function feedbackTest(ctrl,side)
% update feedback pulses based on current dictionary. To make bigger
% changes (trafofn, meas time, pulse spacing, etc.) must edit function directly. 
% function feedbackTest(ctrl)
% ctrl can be all, pump, dbz, dbz2.
%   all: updates all pulses
%   pump the feedback pulses
%   dbz/dbz2 the measuring dbz pulses (dbz2 has larger spacing for faster
%   measurement)

global tuneData; global fbdata; global smdata; 
if ~exist('side','var'), side = tuneData.activeSetName; end

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
        if fbdata.activeQub == 3 && ind == 2 % Set up fb on both sides.
            feedbackTest('pumpBoth');
        elseif fbdata.activeQub == 3 && ind == 1
        else
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
        end
    case 'pumpBoth'
        % This is all pretty hard coded right now.
        
        fbpls = [140:142; 140:142];
        pg.ctrl='loop';
        pg.varpar = 0;
        dicts={'left','right'};
        
        % may be needed.  In order, the 8 pulses are: TT, oT, To, oF, Fo, TF, FT, FF
        % (these are indices into fbpls)
        % switch S for T
        pi1 = 1+[0 0 2 0 1 2 1 1; ...
                 0 2 0 1 0 1 2 1]; % type of pulse, 0,1,2 = [T+,S, off];
        
        pg.nrep = zeros(1,8);
        pg.jump = [];
        for i = 1:2
            pg.dict=dicts{i};
            pg.pulses = fbpls(i, :);
            pg.params = 0;
            pg.pulseind = pi1(i,:);
            pg.name = sprintf('pump%s',upper(dicts{i}(1)));
            if i == 2
                pg.chan = [3 4];
            else
                pg.chan = [2 1];
            end
            plsdefgrp(pg);
        end
        clear pg2;
        %This is the double.
        pg2.name = 'pumpBoth';
        pg2.jump= pg.jump;
        pg2.nrep = pg.nrep;
        pg2.matrix=eye(4);
        pg2.ctrl='';
        pg2.pulses.groups = {'pumpL','pumpR'};
        pg2.chan=[2 1 3 4];
        plsdefgrp(pg2);
        
        awgadd(pg2.name);
        
                % switch S for T
        pi1 = 1+[0 0 2 0 1 2 1 1; ...
                 0 2 0 1 0 1 2 1]; % type of pulse, 0,1,2 = [S,T+, off];       
        seq = awgseqind(pg2.name);
        fbdata.multi.SS=seq+1;
        fbdata.multi.ST=seq+4;
        fbdata.multi.TS=seq+5;
        fbdata.multi.TT=seq+8;
        fbdata.params(1).singletPulse=seq+2;
        fbdata.params(2).singletPulse=seq+3;
        fbdata.params(2).tripletPulse=seq+6;
        fbdata.params(1).tripletPulse=seq+7;
    case 'dbz' % dBz group for software feedback
        if fbdata.activeQub == 3 && ind == 2 % Set up fb on both sides.
            feedbackTest('dbzBoth')
        end
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
    case 'dbzBoth'
        pg.params = [7,130,0];        
        pulses = 128; 
        scale = 1; 
        pg.varpar=((1:(pulses))'*scale);
        pg.pulses = 134;
        pg.ctrl='loop pack';
        
        side='L'; pg.dict={'stagl','left'}; pg.chan=[2 1]; 
        nameL = sprintf('dBz_stag_%d%s',scale*pulses, side);
        pg.name = nameL; 
        plsdefgrp(pg);
        side='R'; pg.dict={'stagr','right'}; pg.chan=[3 4];
        nameR = sprintf('dBz_stag_%d%s',scale*pulses, side);
        pg.name = nameR; 
        plsdefgrp(pg);        
        
        % Combine the two groups above        
        pg2.matrix=eye(4);
        pg2.pulses.groups={nameL,nameR};
        pg2.chan=[2 1 3 4];
        pg2.name=sprintf('dBz_stag_%d_LR',pulses*scale);
        pg2.ctrl='grp loop pack';
        pg2.pulseind(2,:)=1:length(pg.varpar);
        pg2.pulseind(1,:)=1:length(pg.varpar);
        plsdefgrp(pg2);
        awgadd(pg2.name);
        fbdata.params(3).fbInit = pg2.name; 
        
        % Now the long one. 
        scale = 4; 
        pg.params = [8,525,0];        
        pg.varpar=((1:(pulses))'*scale);
        side='L'; pg.dict={'stagl','left'}; pg.chan=[2 1]; 
        nameL = sprintf('dBz_stag_%d%s',scale*pulses, side);
        pg.name = nameL; 
        plsdefgrp(pg);
        side='R'; pg.dict={'stagr','right'}; pg.chan=[3 4];
        nameR = sprintf('dBz_stag_%d%s',scale*pulses, side);
        pg.name = nameR; 
        plsdefgrp(pg);        
        
        pg2.name=sprintf('dBz_stag_%d_LR',pulses*scale);
        pg2.pulses.groups={nameL,nameR};
        plsdefgrp(pg2);
        awgadd(pg2.name);
        
        fbdata.params(3).fbGroup = pg2.name; 
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