%% List all pulses that have been defined 
plslist
%% List details about a particular pulse. 
plsprint(7)
%% adprep and then either adread or read.
%n = length(plsdata.pulses)+1; 
n = 121; 
clear pinf
pinf.name='adtest';
%       1          2       3        4         5        6       7       8        9%         
els={'@start','@wait','@fill','@wait','@reload','@wait','@adprep','@adread','@wait','@meas'};
pinf.data=struct('type',els,'time',[],'val',[]);

% Parameters: [pulselength, adprep/adread time]
pinf.pardef=[3 -1 ; 4 -1 ; 5 -1]; 
pinf.trafofn=@(x) [x(1),x(2),x(2)];
%plsplot(pinf)
plsreg(pinf,n); 
plssync('save'); 
%% Fully configurable dBz pulse
%n = length(plsdata.pulses)+1; 
plsnum = n;
                           %1          2      3         4        5      6   
pinf.data=struct('type',{'@start', '@fill','@wait','@reload','@exch','@meas'},...
                 'time',{  []      ,[]    ,   []      ,  []   , []     ,  [] },...
                 'val' ,{  []      ,[]    ,   []      ,  []   , []     ,  [] });

%params = [pulse length, max sep time, readdout time, sep time]
pinf.pardef = [2 -1; 5 3; 5 -1;];
pinf.trafofn = @(x)[x(1), x(2), x(3)*1e-3];

%plsplot(pinf, 'right')
plsreg(pinf, plsnum);
plssync('save'); 
%% dict based ST+ search STP
plsnum = 8;

clear pinf;
pinf.name = 'ST+';
%       1          2       3          4         5       6
els ={'@start', '@fill', '@wait', '@reload', '@exch', '@meas'};
pinf.data = struct('type', els, 'time', [], 'val', []);
pulselength = 4; %hardcode 4us pulselength
pinf.data(2).time = pulselength;

%params = [loiter time(ns) ,ST+ location (mV)]; 
pinf.pardef = [6,-1; 5, -1; 5, 3];
pinf.trafofn = @(x)[x(1),x(2)*1e-3, x(3)];

%plsplot(pinf, 'right');
plsreg(pinf, plsnum);
%%
clear pls
pinf.name='adtest2';
plsnum=33;
%       1          2       3        4         5       6      7      8  
els={'@start','@reload','@wait','@adprep','@wait','@fill','@wait','@meas'};
pinf.data=struct('type',els,'time',[],'val',[]);

% Parameters: [pulselength, adprep/adread time]
pinf.pardef=[6 -1 ; 4 -1]; 
pinf.trafofn=@(x) [x(1),x(2)];
%plsplot(pinf)
plsreg(pinf,plsnum); 
%% Pulsed Zoom
n = length(plsdata.pulses)+1; 
clear pls 
pls = plsdata.pulses(5); % start by copying old zoom pulse
pls.name = 'PulsedZoom'; 
pls.pardef= [6,-1; 7,2; 7,3];
pls.trafofn=@(x)[x(1), x(3),x(2)];
plsreg(pls,n); 
plssync('save'); 
%% 2D Load 
n = length(plsdata.pulses)+1; 
clear pls 
pls = plsdata.pulses(6); % start by copying old load pulse
pls.name = 'TwoDLoad'; 
pls.pardef= [6,-1; 6,-2; 6,1;6,2];
pls.trafofn = @(x)[x(1)*1e-3,x(2)*1e-3,x(3:4)+[x(5),x(6)]];
plsreg(pls,n); 
plssync('save'); 
%%
n = 117; 
clear pls 
pls = plsdata.pulses(15); 
pls.name = 'RamseyMarker'; 
pls.data(5).type = 'markerburst';
pls.data(5).val = [NaN, NaN, NaN, 1 0 0 0]; 
pls.pardef=     [9,-1; 5,1; 5,2; 5,3; 5,-1; 5 -2; 5 -3; 7 -1]; 
pls.trafofn=@(x)[x(1), x(2),x(3),x(4),x(7)*1e-3,x(5)*1e-3,x(6)*1e-3,.4-x(7)*1e-3];
plsreg(pls,n); 
plssync('save'); 

%% plsInfo 
n = 118; 
clear pls 
pls = plsdata.pulses(15); 
pls.name = 'RamseyMarker2'; 
pls.data(5).type = 'markerburst';
pls.data(5).val = [NaN, NaN, NaN, 0 1 0 0]; 
pls.pardef=     [9,-1; 5,1; 5,2; 5,3; 5,-1; 5 -2; 5 -3; 7 -1; 3 -1; 8 -1]; 
% PlsLen, 3 exch vals, exch time, marker pre, marker delay, wait time. 
pls.trafofn=@(x)[x(1), x(2),x(3),x(4),x(10)*1e-3,x(5)*1e-3,x(6)*1e-3,x(7)-x(10)*1e-3,x(8),x(9)];
plsreg(pls,n); 
plssync('save');

%%
n = 119; 
clear pls 
pls = plsdata.pulses(12); 
pls.name = 'dbzNoFill'; 
pls.data = pls.data(1:5); 
pls.pardef(4,:)=[]; 
% PlsLen, 3 exch vals, exch time, marker pre, marker delay, wait time. 
pls.trafofn=@(x)[(x(1)-x(3))*1e-3,x(3)*1e-3,x(2)];
plsreg(pls,n); 
plssync('save');
%%
dt=40e-3;
pls.data(1).type = '@start'; 
pls.data(2).type = 'markerburst'; 
pls.data(2).time = [dt 0 0]; 
pls.data(2).val = [0 0 0 1 1 0 0]; 
pls.data(3).type = '@wait'; 
pls.data(3).time = dt; 
pls.data(4).type = 'markerburst'; 
pls.data(4).time = [dt 0 0]; 
pls.data(4).val = [0 0 0 1 1 0 0]; 
pls.data(5).type = '@wait'; 
pls.data(5).time = dt; 
pls.data(6).type = 'markerburst'; 
pls.data(6).time = [dt 0 0]; 
pls.data(6).val = [0 0 0 1 1 0 0]; 
pls.data(7).type = '@wait'; 
pls.data(7).time = dt; 
pls.data(8).type = 'markerburst'; 
pls.data(8).time = [dt 0 0]; 
pls.data(8).val = [0 0 0 1 1 0 0]; 
pls.data(9).type = '@wait'; 
pls.data(9).time = dt; 

%plsdata.pulses(116) = pls; 
plsreg(pls,116); 
plssync('save'); 
 awgrm(13,'after'); awgclear('unused');

clear pg 
pg.pulses=116;
pg.ctrl='loop';
pg.chan=[str2double(char(regexp(tuneData.xyChan{1},'\d+','match'))),str2double(char(regexp(tuneData.xyChan{2},'\d+','match')))];
pg.dict={tuneData.activeSetName};

%Parameters: pulselength, eps, evo
%pg.varpar = evo';
%    pg.params=[8 eps(i) 0];
pg.name = 'markerDelay';
plsdefgrp(pg);
awgadd(pg.name);
awgcntrl('on start wait err');