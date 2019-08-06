%% List all pulses that have been defined 
plslist
%% List details about a particular pulse. 
plsprint(7)
% Started this file 2012/08/22
% Brand new everything.
%% 1 Everything off pulse
clear pinf;
plsnum = 1;
els = {'@start', '@fill', '@wait'};
pinf.data= struct('type', els, 'time', [], 'val', []);
pinf.data(2).time = 1;
plsplot(pinf, 'right');
plsreg(pinf, plsnum);
%% 2 Ch. 2 Marker 1 on pulse
clear pinf;
plsnum = 2;
els = {'@start', 'mark','@fill', '@wait'};

pinf.data= struct('type', els, 'time', [], 'val', []);
pinf.data(2).time=[0 0 0 0 1];
pinf.data(3).time = 1;
plsplot(pinf, 'right');
plsreg(pinf, plsnum);
%% 6 Dict based load, with varpars that help fitting, for load pos and time check.
plsnum = 6;
clear pinf;
pinf.name = 'reload';
%       1           2       3        4         5         6        7
els = {'@start', '@fill', '@wait', '@rand', '@wait', '@reload', '@meas'};
pinf.data = struct('type', els, 'time', [], 'val', []);
pulselength = 6; %hardcode 6us pulselength
measT = 1; %hard code measurment time;

pinf.data(2).time = pulselength;
pinf.data(5).time(1)=0.01;
pinf.data(7).time(1) = measT;

%params=[ramp to/from load (ns), loadTime (ns), load cntr loadpos offset(mV)]
pinf.pardef = [6 -1; 6 -2; 6 1; 6 2];
pinf.trafofn = @(x) [x(1)*1e-3, x(2)*1e-3, x(3:4)+x(5)*[1 1]];

plsplot(pinf, 'right');
plsreg(pinf, plsnum);
%% 7 top lead search (TL) dict based
plsnum = 7;
clear pinf;
pinf.name = 'TL';
%       1           2       3        4          5         6
els = {'@start', '@fill', '@wait', '@reload', '@sep', '@meas'};
pinf.data = struct('type', els, 'time', [], 'val', []);
pulselength = 4; %hardcode 4us pulselength
measT = 1; %hard code measurment time;

pinf.data(2).time = pulselength;
pinf.data(6).time(1) = measT;
pinf.data(5).time = 1;  %hard code 1 us loiter time

%params = [cntrX, cntrY, dirX, dirY, eps(uV)];
pinf.pardef = [5 1; 5 2];
pinf.trafofn = @(x) x(1:2)+(x(3:4)/norm(x(3:4)))*x(5)*1e-3;

%plsplot(pinf, 'right');
plsreg(pinf, plsnum);
%% 8 dict based ST+ search STP
plsnum = 8;

clear pinf;
pinf.name = 'ST+';
%       1          2       3          4         5       6
els ={'@start', '@fill', '@wait', '@reload', '@exch', '@meas'};
pinf.data = struct('type', els, 'time', [], 'val', []);
pulselength = 4; %hardcode 4us pulselength
measT = 1; %hard code measurment time;
pinf.data(2).time = pulselength;
pinf.data(6).time(1) = measT;

%params = [loiter time(ns) ,ST+ location (mV)]; 
pinf.pardef = [5 -1; 5 3];
pinf.trafofn = @(x)[x(1)*1e-3, x(2)];

%plsplot(pinf, 'right');
plsreg(pinf, plsnum);
%% 8 dict based ST+ search STP
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
%% 9 Singlet feedback polarizer, dict based, brave new world feedback.
% dict-based feedback elements.
% expects supplemental dictionary to define @stpsweep
plsnum = 9;
clear pinf;
pinf.name = 'Spol';
%       1          2       3         4       5
els ={'@start', '@fill', '@reload', '@sep', '@stpsweep'};
pinf.data = struct('type', els, 'time', [], 'val', []);
pinf.data(2).time = .25; %250ns pulse length
pinf.data(4).time = .001; %1ns loiter at far eps
pinf.data(3).time = [.005 .08 .005]; % hard code short load; this will get filled out to longer.
pinf.pardef = [4 -1];
pinf.xval=1;
%pinf.trafofn=@(x) [x(1)*1e-3];
%Parameters; sep time in ns
plsplot(pinf, {'right'});
plsreg(pinf,plsnum);
%% 10 Triplet unconditional polarizer, dict based, brave new world feedback
% dict-based feedback elements.
% expects supplemental dictionary to define @stpsweep, @tlsweep
% ignorse supplied parameters.
plsnum = 10;
clear pinf;
pinf.name = 'Tpol';
%       1          2       3         4       5
els ={'@start', '@fill', '@wait', '@tlsweep','@stpsweep'};
pinf.data = struct('type', els, 'time', [], 'val', []);
pinf.data(2).time = .250; %250ns pulse length
pinf.data(3).time=1e-3;
pinf.pardef = [];
pinf.trafofn=@(x) [];
pinf.xval=0;
%Parameters; sep time in ns
plsplot(pinf, {'right'});
plsreg(pinf,plsnum);
%% 11 250ns filler pulse for feedback
clear pinf;

pinf.data.pulsetab = [0 .25;0 0;0 0];
pinf.name = sprintf('fill_%d', round(plsdata.tbase * pinf.data.pulsetab(1, 2)));
plsreg(pinf, 11);
%% 12 Fully configurable dBz pulse
plsnum = 12;
                            %1          2       3         4        5
pinf.data=struct('type',{'@start', '@wait',  '@reload', '@sep', '@meas','fill','@wait'},...
                 'time',{  []        ,[]  ,   []      ,  []   , []     ,  []  , []},...
                 'val' ,{  []        , [] ,   []      ,  []   ,[]      ,  []  , []});
%params = [pulse length, max sep time, readdout time, sep time]
pinf.pardef = [2 -1; 4 -1; 5 -1; 6 -1];
pinf.trafofn = @(x)[(x(2)-x(4))*1e-3, x(4)*1e-3, x(3), x(1)];
%plsplot(pinf, 'right')
plsreg(pinf, plsnum);
%% 15 Plain old ramsey, adjustable prep
plsnum = 15;

clear pinf;
pinf.name='Ramsey';
%     1         2         3        4       5         6      7      8       9
els={'@start','@wait','@reload','@wait','@prep','@exch','@read','@wait', '@meas','@fill','@wait'};
pinf.data=struct('type',els,'time',[],'val',[]);
% Parameters: p(1) = pulse length p(2) = epsilon (mv) ; P(3) = evo time
pinf.pardef = [9 -1; 5 3; 5 -1; 7 -1]; 
pinf.data(7).val(3)=0; %enforce wait at 0;

pinf.trafofn = @(x) [x(1) x(2) x(3)*1e-3 .4-x(3)*1e-3];
%plsplot(pinf, {struct('read','@adread','prep','@adprep'),'left'});
plsreg(pinf, plsnum);
%% 16 autotune(read) scan

measT = 7; %read for 7 mus. 

plsnum = 16;

clear pinf;
pinf.name='Read';
%     1         2         3        4       5        6      7         8       9
els={'@start','@fill', '@wait','@reload','@wait','@rand','@wait', '@meas'};
pinf.data=struct('type',els,'time',[],'val',[]);
pinf.data(8).time(1) = measT;
% Parameters: p(1) = pulse length p(2) = wait time at 0,0 ;
pinf.pardef = [2 -1; 7 -1]; 

pinf.trafofn = @(x) [x(1) x(2)];
%plsplot(pinf, {struct('prep','@rand'),'left'});
plsreg(pinf, plsnum);

plsnum = 17;

clear pinf;
pinf.name='Read2';
%     1         2         3        4       5        6      7         8       9
els={'@start','@fill', '@wait','@reload','@wait', '@meas'};
pinf.data=struct('type',els,'time',[],'val',[]);
pinf.data(6).time(1) = measT;
% Parameters: p(1) = pulse length p(2) = wait time at 0,0 ;
pinf.pardef = [2 -1; 5 -1]; 

pinf.trafofn = @(x) [x(1) x(2)];
%plsplot(pinf, {struct('prep','@rand'),'left'});
plsreg(pinf, plsnum);
%% 18 Adiabatic-style stp pulse

plsnum = 18;

clear pinf;
pinf.name='AdSTP';
%     1         2         3        4         5       6      7      8     
els={'@start','@reload','@wait','@adprep','@exch','@meas','@fill','@wait'};
pinf.data=struct('type',els,'time',[],'val',[]);
% Parameters: p(1)=pulse length, p(2)=adprep time, p(3)=start_eps, p(4)=loiter time after adprep, p(5)=eps
pinf.pardef = [7 -1; 4 -1; 4 1; 4 2; 5 -1; 5 3]; 

pinf.trafofn = @(x) [x(1), x(2), x(3), x(5), x(4), x(5)];
plsplot(pinf,'left');
plsreg(pinf, plsnum);
%% 19 Line-scan pulse for determining tunnel coupling
plsnum = 19;
clear pinf;
pinf.name='Line';
%       1         2        3        4       5      6
els={'@start','@reload','@wait','@meas','@wait','@fill','@wait'};
% Parameters: p(1)=pulse length, p(2)=meas time p(3:4) = measure location
pinf.data=struct('type',els,'time',[],'val',[]);
pinf.data(6).time(1)=2;
pinf.data(3).time=0.001;
pinf.data(5).time=0.001;
pinf.data(4).time(1)=1;
pinf.data(4).val=[nan 1 -1];
pinf.pardef=[6 -1; 4 -1; 4 2 ; 4 3];
pinf.trafofn=@(x) x;
plsplot(pinf,'left');
plsreg(pinf,plsnum);

% Good pardef is [2 1 0 0], 2xn varpar
%% 20 The good T1 (for use with autotune)

% dBz, 15us readout time for T1 checks, parameterized by rotation time.
% Uses new pulse dictionary.
% Parameterized by measp, dBz rotation time.
clear pinf;
pinf.name='dBz';
pulsenum = 20;
%                           1        2       3       4          5      6                              
pinf.data=struct('type',{'@start','@fill','@wait','@reload','@sep','@meas'},...
                 'time',{[]      ,[17]   ,[]     ,[]       ,[.1]  ,[15]   },...
                 'val' ,{[1 1]      ,[]     ,[1 1]     ,[]       ,[]     ,[nan 1 1]     });
pinf.pardef = [5, -1; 6, 2; 6,3; 3,1 ; 3 ,2; 1,1 ; 1,2];
pinf.trafofn = @(x) [x(3)*1e-3 x(1) x(2) x(1) x(2) x(1) x(2)];
plsreg(pinf,pulsenum);
%plsplot(pinf,'right');
%% 21 Adiabatic-style st+ pulse, remove background.

plsnum = 21;

clear pinf;
pinf.name='AdSTP';
%     1         2         3        4         5       6      7      8     
els={'@start','@reload','@wait','@adprep','@adread','@meas','@fill','@wait'};
pinf.data=struct('type',els,'time',[],'val',[]);
% Parameters: p(1)=pulse length, p(2)=adprep time, p(3)=start_eps, p(4)=eps
pinf.pardef = [7 -1; 4 -1; 4 1; 4 2; 5 -1; 5 1; 5 2]; 

pinf.trafofn = @(x) [x(1), x(2), x(3), x(4), x(2), x(3), x(4)];
plsplot(pinf,'left');
plsreg(pinf, plsnum);
%% 22 Plain old ramseyE, either pi or dbz readout and prep. short pulse for awg7, parameterized length.
plsnum = 22;

clear pinf;
pinf.name='RamseyE';
%     1       2         3        4       5         6       7      8     9       10       11      12  
els={'@start','@reload','@wait','@prep','@exch','@pi','@exch','@read','@wait','@meas','@fill','@wait'};
pinf.data=struct('type',els,'time',[],'val',[]);
pinf.data(11).time = 4;
% Parameters: p(2) = epsilon (mv) ; p(3)=total evo time ; p(4) = dt in nsec
% p(1) is pulse length in us
pinf.pardef = [5 -1; 7 -1; 9 -1; 5 3 ; 7 3; 11 -1]; 

pinf.trafofn = @(x) [x(3)/2, x(3)/2+x(4)*1e-3, .3-x(4)*1e-3, x(2), x(2), x(1)];
%plsplot(pinf, {struct('read','@dbzprep','prep','@dbzprep','pi','@dbzpi'),'right'});
plsreg(pinf, plsnum);
%% 23 Dict based load test, for zoom scans
plsnum = 23;
clear pinf;
pinf.name = 'zoom';
%       1           2       3        4         5         6       7
els = {'@start', '@fill','@wait','@rand', '@wait', '@reload', '@meas'};
pinf.data = struct('type', els, 'time', [], 'val', []);

pinf.data(2).time = 2; %pulse length
pint.data(5).time=0.01;
pinf.data(6).time=[0.001 0.5 0.001];
pinf.data(7).time(1)=0.9;

%params=[ramp to/from load (ns), loadTime (ns), load cntr loadpos offset(mV)]
pinf.pardef = [6 -2];
pinf.trafofn = @(x) x;

plsplot(pinf, 'left');
plsreg(pinf, plsnum);
%% 24 dict based ST+ search STP
plsnum = 24;

clear pinf;
pinf.name = 'ST+';
%       1          2       3          4         5       6
els ={'@start', '@fill', '@wait', '@reload', '@exch', '@meas'};
pinf.data = struct('type', els, 'time', [], 'val', []);
pulselength = 4; %hardcode 4us pulselength
measT = 1; %hard code measurment time;
pinf.data(2).time = pulselength;
pinf.data(6).time(1) = measT;

%params = [pulse length, loiter time(ns) ,ST+ location (mV)]; 
pinf.pardef = [2 -1 ; 5 -1; 5 3];
pinf.trafofn = @(x)[x(1), x(2)*1e-3, x(3)];

%plsplot(pinf, 'right');
plsreg(pinf, plsnum);
%% 25 Pulse full of zeros

plsnum = 25;
clear pinf;
els= {'@start','@fill','@wait'};
% one parameter is fill time (is us)
pinf.data = struct('type',els,'time',[],'val',[]);

pinf.pardef = [2 -1];
pinf.trafofn = @(x) x;
plsreg(pinf,plsnum);
%% 27/28 Adprep ra msey designed to sync for controlled evution
plsnum=28;   %adprep and adread
plsnum2=27;  %dbz
clear pinf;
clear pinf2;
pinf.name='Ramsey_adprep_controlledevo';
pinf2.name='dBz_controlledevo';
%      1          2        3         4       5         6      7       8     9         10
els={'@start','@reload','@adprep','@sep','@exch','@adread','@wait','@fill','@wait','@meas','@wait'};
els2={'@start','@reload','@wait','@sep','@wait','@fill','@wait','@meas'};
pinf.data=struct('type',els,'time',[],'val',[]);
pinf2.data=struct('type',els2,'time',[],'val',[]);
pinf2.data(3).time=0.61;
pinf.data(4).time=0.1;
pinf.data(7).time=0.001;
pinf2.data(5).val = [0 0]; % to make the wait happen at measp
pinf.pardef=[8 -1 ; 5 3 ; 5 -1];
pinf2.pardef=[6 -1 ; 4 -1];
pinf.trafofn=@(x) x.*[1 1 1e-3];  % pulse time, eps, exch time
pinf2.trafofn=@(x) x.*[1 1e-3];   % pulse time, sep time.
plsreg(pinf,plsnum);
plsreg(pinf2,plsnum2);
%plsplot(pinf2,'right');
%% 29 RamseyE with tomo. not extra dtau rotation
plsnum = 29;

clear pinf;
pinf.name='RamseyE';
%     1       2         3        4       5         6       7      8     9       10       11      12  
els={'@start','@reload','@wait','@prep','@exch','@pi','@exch','@tomo','@wait','@fill','@wait','@meas','@wait'};
pinf.data=struct('type',els,'time',[],'val',[]);

%params = [pulse length, eps, T_evo(us)]

pinf.pardef = [5 -1; 7 -1; 5 3 ; 7 3; 10 -1]; 

pinf.trafofn = @(x) [x(3)/2, x(3)/2, x(2), x(2), x(1)];
%plsplot(pinf, {struct('read','@dbzprep','prep','@dbzprep','pi','@dbzpi'),'right'});
plsreg(pinf, plsnum);
%% 30 controlled evo pulse that allows wait at an epsilon
%this will replace pulse 27
plsnum2=30;  %dbz
clear pinf;
clear pinf2;

pinf2.name='dBz_controlledevo';
%      1          2        3       4       5       6      7       8     9         10

els2={'@start','@reload','@wait','@sep','@exch','@fill','@exch','@meas'};
%pinf.data=struct('type',els,'time',[],'val',[]);
pinf2.data=struct('type',els2,'time',[],'val',[]);
pinf2.data(3).time=0.61;
%pinf.data(4).time=0.1;
%pinf.data(7).time=0.001;
%
%pinf.pardef=[8 -1 ; 5 3 ; 5 -1];
pinf2.pardef=[6 -1 ; 5 3; 7 3; 4 -1];
%pinf.trafofn=@(x) x.*[1 1 1e-3];  % pulse time, eps, exch time
pinf2.trafofn=@(x) [x(1), x(2), x(2), x(3)*1e-3];   % pulse time, epsilon, sep time.
%plsreg(pinf,plsnum);
plsreg(pinf2,plsnum2);
%plsplot(pinf2,'right');
%% 31 state prep, exch. evolve, then tomo
plsnum = 31;

clear pinf;
pinf.name='Tomo';
%     1         2         3        4       5       6       7      
els={'@start','@reload','@wait','@prep','@exch', '@tomo','@meas','@fill','@wait'};
pinf.data=struct('type',els,'time',[],'val',[]);
pinf.data(8).time = 10.0;
% Parameters: p(1) J, p(2) wait time
pinf.pardef = [3 -1 ; 5 -1; 5 3]; 

pinf.trafofn = @(x) [.1-x(2)*1e-3, x(2)*1e-3, x(1)];
%plsplot(pinf, {struct('prep','@dbzpi', 'tomo', '@UDread'),'right'});
plsreg(pinf, plsnum);
%% 32 adprep and then either adread or read.
clear pinf
pinf.name='adtest';
plsnum=32;
%       1          2       3        4         5        6       7       8%         
els={'@start','@reload','@wait','@adprep','@adread','@wait','@fill','@wait','@meas'};
pinf.data=struct('type',els,'time',[],'val',[]);

% Parameters: [pulselength, adprep/adread time]
pinf.pardef=[7 -1 ; 4 -1 ; 5 -1]; 
pinf.trafofn=@(x) [x(1),x(2),x(2)];
%plsplot(pinf)
plsreg(pinf,plsnum); 

clear pinf
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
%% 33
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
%% 34 state prep, exch. evolve, then tomo, adding extra wait at the end to help with stag and combining groups
plsnum = 34;

clear pinf;
pinf.name='Tomo';
%     1         2         3        4       5       6       7      
els={'@start','@reload','@wait','@prep','@exch', '@tomo','@meas','@wait','@fill','@wait'};
pinf.data=struct('type',els,'time',[],'val',[]);
pinf.data(9).time = 10.0;
pinf.data(8).time = .002;
pinf.data(8).val = [0 0];
% Parameters: p(1) J, p(2) wait time
pinf.pardef = [3 -1 ; 5 -1; 5 3]; 

pinf.trafofn = @(x) [.1-x(2)*1e-3, x(2)*1e-3, x(1)];
%plsplot(pinf, {struct('prep','@dbzpi', 'tomo', '@UDread'),'right'});
plsreg(pinf, plsnum);
%% 35 Testing readout time. 
clear pinf
pinf.name='readouttest';
plsnum=35;
%       1          2       3        4        5      6       7       8%         
els={'@start','@wait','@fill','@wait','@reload','@prep','@wait','@meas'};
pinf.data=struct('type',els,'time',[],'val',[]);

% Parameters: [pulselength, wait time]
pinf.pardef=[3 -1 ; 7 -1]; 
pinf.trafofn=@(x) [x(1),x(2)];
%plsplot(pinf)
plsreg(pinf,plsnum); 
%% 36 RamseyE with tomo. fixing 29 which has too long wait bef rdout
plsnum = 36;

clear pinf;
pinf.name='RamseyE';
%     1       2         3        4       5         6       7      8     9       10       11      12  
els={'@start','@reload','@wait','@prep','@exch','@pi','@exch','@tomo','@meas','@wait','@fill','@wait'};
pinf.data=struct('type',els,'time',[],'val',[]);
pinf.data(10).time = .005;
%params = [pulse length, eps, T_evo(us)]

pinf.pardef = [5 -1; 7 -1; 5 3 ; 7 3; 11 -1]; 

pinf.trafofn = @(x) [x(3)/2, x(3)/2, x(2), x(2), x(1)];
%plsplot(pinf, {struct('read','@dbzprep','prep','@dbzprep','pi','@dbzpi'),'right'});
plsreg(pinf, plsnum);
%% 37 RamseyE with tomo. fixing 29 which has too long wait bef rdout
plsnum = 37;

clear pinf;
pinf.name='RamseyE';
%     1       2         3        4       5         6       7      8     9       10       11      12  
els={'@start','@wait','@reload','@wait','@prep','@exch','@pi','@exch','@tomo','@meas','@wait','@fill','@wait'};
pinf.data=struct('type',els,'time',[],'val',[]);
pinf.data(11).time = .005;
%params = [pulse length, max evo, eps, T_evo(us)]

pinf.pardef = [12 -1; 2 -1; 6 -1; 8 -1; 6 3 ; 8 3]; 

pinf.trafofn = @(x) [x(1), x(2)-x(4), x(4)/2, x(4)/2, x(3), x(3)];
%plsplot(pinf, {struct('read','@dbzprep','prep','@dbzprep','pi','@dbzpi'),'right'});
plsreg(pinf, plsnum);
%% 38 Generic prep exchange read
clear pinf
pinf.name='generic';
plsnum=38;
%       1          2       3       4        5       6       7       8%         
els={'@start','@wait','@fill','@wait',  '@reload','@prep','@exch', '@read','@meas'};
pinf.data=struct('type',els,'time',[],'val',[]);

% Parameters: [pulselength, exch eps, exch time]
pinf.pardef=[3 -1 ; 7 3; 7 -1]; 
pinf.trafofn=@(x) [x(1),x(2), x(3)*1e-3];
%plsplot(pinf)
plsreg(pinf,plsnum);
%% 41 Load, then hold at some point. 
plsnum = 41;

clear pinf;
pinf.name='hold';
%     1         2         3        4       5 
els={'@start','@reload','@wait','@fill','@wait'};
pinf.data=struct('type',els,'time',[],'val',[]);
% Parameters: p(1) = pulse length p(2) = epsilon (mv) ; P(3) = evo time
pinf.pardef = [4 -1; 3 3; 5 3]; 

pinf.trafofn = @(x) [x(1) x(2) x(2)];
plsreg(pinf, plsnum);
%% 42 Pulse to hold qubit at some eps. 

plsnum = 42;
clear pinf;
        %1       2        3       4
els= {'@start','@exch','@fill','@exch'};
% params = fill time, eps
pinf.data = struct('type',els,'time',[],'val',[]);

pinf.pardef = [3 -1; 2 3; 4 3];
pinf.trafofn = @(x) [x(1) x(2) x(2)];
plsreg(pinf,plsnum);
%% 43 Pulse to hold qubit at some eps. 

plsnum = 43;
clear pinf;
        %1       2        3       4
els= {'@start','@fill','@exch'};
% one parameter is fill time (is us)
pinf.data = struct('type',els,'time',[],'val',[]);

pinf.pardef = [2 -1; 3 3];
pinf.trafofn = @(x) x;
plsreg(pinf,plsnum);
%% 43 Fully configurable dBz pulse with wait at the end. 

plsnum = 43;

                            %1         2       3         4        5    6
pinf.data=struct('type',{'@start','@fill','@wait',  '@reload', '@sep', '@meas'},...
                 'time',{  []        ,[]  ,    []   , []     ,  []  , []},...
                 'val' ,{  []        , [] ,    []   ,[]      ,  []  , []});

%params = [pulse length, max sep time, readdout time, sep time]
pinf.pardef = [3 -1; 5 -1; 6 -1; 2 -1];
pinf.trafofn = @(x)[(x(2)-x(4))*1e-3, x(4)*1e-3, x(3), x(1)];

%plsplot(pinf, 'right')
plsreg(pinf, plsnum);
%% 44 exchange with parameterizable x and y (not just eps). 

plsnum = 44;

%     1         2         3        4       5         6      7      8       9
els={'@start','@fill','@reload' '@wait','@prep','@exch','@read','@wait', '@meas'};
pinf.data=struct('type',els,'time',[],'val',[]);

%params = [pulse length, exch x, exch y, exch time]
pinf.pardef = [2 -1; 6 -1; 6 1; 6 2];
pinf.trafofn = @(x)[x(1), x(4)*1e-3, x(2) x(3)];

%plsplot(pinf, 'right')
plsreg(pinf, plsnum);
%% 45 dict based ST+ search STP
plsnum = 45;

clear pinf;
pinf.name = 'ST+';
%       1          2       3          4         5       6
els ={'@start', '@fill', '@wait', '@reload', '@exch45', '@meas'};
pinf.data = struct('type', els, 'time', [], 'val', []);
pulselength = 4; %hardcode 4us pulselength
measT = 1; %hard code measurment time;
pinf.data(2).time = pulselength;
pinf.data(6).time(1) = measT;

%params = [loiter time(ns) ,ST+ location (mV)]; 
pinf.pardef = [5 -1; 5 3];
pinf.trafofn = @(x)[x(1)*1e-3, x(2)];

%plsplot(pinf, 'right');
plsreg(pinf, plsnum);
%% 46 Singlet feedback polarizer, dict based, brave new world feedback.
% dict-based feedback elements.
% expects supplemental dictionary to define @stpsweep
plsnum = 46;
clear pinf;
pinf.name = 'Spol2';
%       1          2       3         4       5
els ={'@start', '@fill', '@reload', '@sep45', '@stpsweep'};
pinf.data = struct('type', els, 'time', [], 'val', []);
pinf.data(2).time = .25; %250ns pulse length
pinf.data(4).time = .001; %1ns loiter at far eps
pinf.data(3).time = [.005 .08 .005]; % hard code short load; this will get filled out to longer.
pinf.pardef = [4 -1];
pinf.xval=1;
%pinf.trafofn=@(x) [x(1)*1e-3];
%Parameters; sep time in ns
plsplot(pinf, {'right'});
plsreg(pinf,plsnum);
%% 47 Fully configurable dBz pulse w/ sep 45 to check feedback. 
plsnum = 47;
                           %1          2       3         4        5
pinf.data=struct('type',{'@start', '@wait',  '@reload', '@sep45', '@meas','fill','@wait'},...
                 'time',{  []        ,[]  ,   []      ,  []   , []     ,  []  , []},...
                 'val' ,{  []        , [] ,   []      ,  []   ,[]      ,  []  , []});

%params = [pulse length, max sep time, readdout time, sep time]
pinf.pardef = [2 -1; 4 -1; 5 -1; 6 -1];
pinf.trafofn = @(x)[(x(2)-x(4))*1e-3, x(4)*1e-3, x(3), x(1)];

%plsplot(pinf, 'right')
plsreg(pinf, plsnum);
%% 50 Dict based load test, for zoom scans (parameterizable meas time)
plsnum = 50;
clear pinf;
pinf.name = 'zoom';
%       1           2       3        4         5         6       7
els = {'@start', '@fill','@wait','@rand', '@wait', '@reload', '@meas'};
pinf.data = struct('type', els, 'time', [], 'val', []);

pinf.data(2).time = 2; %pulse length
pint.data(5).time=0.01;
pinf.data(6).time=[0.001 0.5 0.001];


% params: pulse len, meas T, reload T (so can set to zero);
pinf.pardef = [2 -1; 7 -1; 6 -2];
pinf.trafofn = @(x) x;

%plsplot(pinf, 'left');
plsreg(pinf, plsnum);
%% 51 Fully configurable dBz pulse with forced 1,2 readout (looks like pls 12)

plsnum = 51;

                            %1          2       3         4        5
pinf.data=struct('type',{'@start', '@wait',  '@reload', '@sep', '@meas12','fill','@wait'},...
                 'time',{  []        ,[]  ,   []      ,  []   , []     ,  []  , []},...
                 'val' ,{  []        , [] ,   []      ,  []   ,[]      ,  []  , []});

%params = [pulse length, max sep time, readdout time, sep time]
pinf.pardef = [2 -1; 4 -1; 5 -1; 6 -1];
pinf.trafofn = @(x)[(x(2)-x(4))*1e-3, x(4)*1e-3, x(3), x(1)];

%plsplot(pinf, 'right')
plsreg(pinf, plsnum);
%% 69 Singlet feedback polarizer, dict based, brave new world feedback.
% dict-based feedback elements.
% expects supplemental dictionary to define @stpsweep
% this will shorten load at all costs
plsnum = 69;
clear pinf;
pinf.name = 'Spol2';
%       1          2       3         4       5
els ={'@start', '@fill', '@reload', '@sep45', '@stpsweep'};
pinf.data = struct('type', els, 'time', [], 'val', []);
pinf.data(2).time = .25; %250ns pulse length
pinf.data(4).time = .001; %1ns loiter at far eps
pinf.data(3).time = [.001 .03 .001]; % hard code short load; this will get filled out to longer.
pinf.pardef = [4 -1];
pinf.xval=1;
%pinf.trafofn=@(x) [x(1)*1e-3];
%Parameters; sep time in ns
%plsplot(pinf, {'right'});
plsreg(pinf,plsnum);
%% 68 fb crazy
plsnum = 68;
clear pinf;
pinf.name = 'FB_crazy';
%       1          2       3         4       5
els ={'@start', '@fill', '@reload', '@sep', '@tlsweep','@stpsweep'};
pinf.data = struct('type', els, 'time', [], 'val', []);
pinf.data(3).time=1e-3;
pinf.pardef = [2 -1; 4 -1; 5 -1];
pinf.trafofn=@(x) [x];

%Parameters; pulse length, sep time, tlsweep time 
plsreg(pinf,plsnum);
%% 71 Fully configurable dBz pulse, with extra markers on readout

plsnum = 71;
          %1          2       3         4        5
els={ 'mark',   '@start', '@wait',  '@reload', '@sep', '@meas','fill','@wait'};
pinf.data=struct('type',els,'time',[],'val',[]);
%params = [pulse length, max sep time, readdout time, marker_fire_start, marker_fire_duration, sep time]
pinf.pardef = [1 -1; 1 -2; 1 -3; 1 -4; 1 -5; 3 -1; 5 -1; 6 -1; 7 -1];
pinf.trafofn = @(x)[x(4),0,0,x(5),0,(x(2)-x(6))*1e-3, x(6)*1e-3, x(3), x(1)];

%plsplot(pinf, 'right')
plsreg(pinf, plsnum);
%% 76 check T1 of the T+ state (or try)
plsnum = 76;

clear pinf;
pinf.name='T1T+';
%     1       2         3        4       5         6       7      8     9       10       11      12  
els={'@start','@fill','@wait','@wait','@reload','@prep','@read','@meas'};
pinf.data=struct('type',els,'time',[],'val',[]);
pinf.data(2).time = 17;
pinf.data(8).time = 15;
pinf.data(1).val = [1 1];
pinf.data(8).val = [NaN, 1, 1];
% Parameters: 
pinf.pardef = [4 -1]; 

pinf.trafofn = @(x) x;
%plsplot(pinf, {struct('read','@dbzprep','prep','@dbzprep','pi','@dbzpi'),'right'});
plsreg(pinf, plsnum);
%% 89 Dict based load test, for zoom scans, with measurement time control. 
plsnum = 89;
clear pinf;
pinf.name = 'zoom_time';
%       1           2       3        4         5         6       7
els = {'@start', '@fill','@wait','@rand', '@wait', '@reload', '@meas'};
pinf.data = struct('type', els, 'time', [], 'val', []);

pinf.data(2).time = 5; %pulse length
pint.data(5).time=0.01;
pinf.data(6).time=[0.001 0.5 0.001];
%pinf.data(7).time(1)=0.9;

%params=[ramp to/from load (ns), loadTime (ns), load cntr loadpos offset(mV)]
pinf.pardef = [7 -1; 6 -2];
pinf.trafofn = @(x) x;

plsplot(pinf, 'left');
plsreg(pinf, plsnum);
%% 90 Fully configurable ramsey pulse, with extra markers on readout

plsnum = 90;
pinf.name='ramseyramsey'; 

%     1         2        3        4         5       6       7        8     9  
els={'mark','@start','@wait', '@reload','@wait','@prep','@prep2','@exch','@read', '@read2','@meas','@fill','@wait'};
pinf.data=struct('type',els,'time',[],'val',[]);

%               1             2            3   
% Parameters: [pulselength, wait_eps, wait_time(ns),] 
pinf.pardef = [1 -1; 1 -2; 1 -3; 1 -4; 1 -5; 12 -1; 8 3; 8 -1]; 

pinf.trafofn = @(x) [x(3)+x(5)*1e-3,0,0,x(4),0,x(1), x(2), x(5)*1e-3];
%plsplot(pinf, {struct('read','@adread','prep','@adprep'),'left'});
plsreg(pinf, plsnum);
%% 93 Fully configurable dBz pulse with marker during sep (to add to sep)

plsnum = 93;

                            %1          2       3         4        5
pinf.data=struct('type',{'mark', '@start', '@wait',  '@reload', '@sep', '@meas','fill','@wait'},...
                 'time',{[],      []        ,[]  ,   []      ,  []   , []     ,  []  , []},...
                 'val' ,{[],      []        , [] ,   []      ,  []   ,[]      ,  []  , []});

%params = [pulse length, max sep time, readdout time, mark_start, sep time]
pinf.pardef = [3 -1; 5 -1; 6 -1; 7 -1; 1 -1; 1 -2; 1 -3; 1 -4; 1 -5];
pinf.trafofn = @(x)[(x(2)-x(5))*1e-3, x(5)*1e-3, x(3), x(1), x(4), x(5)*1e-3*[0, 1, 1, 0]];

%plsplot(pinf, 'right')
plsreg(pinf, plsnum);
%% 95 Plain old ramsey, adjustable prep, adjustable exch voltage

plsnum = 95;

clear pinf;
pinf.name='Ramsey';
%     1         2         3        4       5         6      7      8       9
els={'@start','@reload','@wait','@prep','@exch','@read','@wait', '@meas','@fill','@wait'};
pinf.data=struct('type',els,'time',[],'val',[]);
% Parameters: p(1) = pulse length p(2) = voltage1(mv) p(3)=voltage2 mv  ; P(4) = evo time
pinf.pardef = [9 -1; 5 1; 5 2; 5 -1; 7 -1]; 

pinf.trafofn = @(x) [x(1) x(2) x(3) x(4)*1e-3 .4-x(4)*1e-3];
plsplot(pinf, {struct('read','@adread','prep','@adprep'),'left'});
plsreg(pinf, plsnum);
%% 96 Plain old ramsey, adjustable prep, adjustable exch voltage 1 flip

plsnum = 96;

clear pinf;
pinf.name='Ramsey';
%     1         2         3        4       5         6     7      8      9       10
els={'@start','@reload','@wait','@prep','@exch','@exch','@read','@wait', '@meas','@fill','@wait'};
pinf.data=struct('type',els,'time',[],'val',[]);
% Parameters: p(1) = pulse length p(2) = ex1_1(mv) p(3)=ex1_2 mv  p(4)ex2_1 (mv) p(5)ex2_2; P(6) = evo time
pinf.pardef = [10 -1; 5 1; 5 2; 5 -1; 6 1; 6 2; 6 -1; 3 -1]; 

pinf.trafofn = @(x) [x(1) x(2) x(3) x(6)*1e-3/2 x(4) x(5) x(6)*1e-3/2 .4-(x(6)*1e-3)];
plsplot(pinf, {struct('read','@adread','prep','@adprep'),'left'});
plsreg(pinf, plsnum);
%% 97 Ramsey with generalized evo and controllable gate voltages 8 flips.
plsnum = 97;

clear pinf;
pinf.name='Ramsey';
%     1         2         3        4       5       6     7      8      9       10
els={'@start','@fill','@wait','@reload','@wait','@prep',};
for j = 1:16
    els{end+1} = '@exch';
end
els = [els,{'@read','@wait', '@meas'}];
pinf.data=struct('type',els,'time',[],'val',[]);
% Parameters: p(1) = pulse length, eps_x1, eps_y1, eps_x2, epsy2, time 
pinf.pardef = [2 -1;];
pinf.pardef = zeros(3*16+1,2);
pinf.pardef(1,:) = [2 -1];
pinf.pardef(2:9,:) = [(7:2:21)',ones(8,1)];
pinf.pardef((2:9)+8,:) = [(7:2:21)',2*ones(8,1)];
pinf.pardef((2:9)+8*2,:) = [(8:2:22)',ones(8,1)];
pinf.pardef((2:9)+8*3,:) = [(8:2:22)',2*ones(8,1)];
pinf.pardef((2:17)+8*4,:) = [(7:22)',-1*ones(16,1)];

pinf.trafofn = @(x) [x(1), x(2)*ones(8,1)',x(3)*ones(8,1)',x(4)*ones(8,1)',x(5)*ones(8,1)',1e-3*(1/16)*x(6)*ones(16,1)',];

plsreg(pinf, plsnum);
%% 99 Ramsey, 1 flip, with echo
plsnum = 99;

clear pinf;
pinf.name='RamseyE';
%     1          2        3        4      5        6      7      8        9      10       11      12     13      14
els={'@start','@fill','@wait','@reload','@wait','@prep','@exch','@exch','@pi','@exch','@exch','@read','@wait','@meas'};
pinf.data=struct('type',els,'time',[],'val',[]);
% Parameters: p(2) = epsilon (mv) ; p(3)=total evo time ; p(4) = dt in nsec
% p(1) is pulse length in us
% Parameters: p(1) = pulse length p(2) = ex1_1(mv) p(3)=ex1_2 mv  p(4)ex2_1
% (mv) p(5)ex2_2; P(6) = evo time; p(7) dt in ns
pinf.pardef = [2 -1; 7 1; 11 1; 7 2; 11 2; 8 1; 10 1; 8 2; 10 2; 7 -1; 8 -1; 10 -1; 11 -1;]; 

pinf.trafofn = @(x) [x(1), x(2), x(2), x(3), x(3), x(4), x(4), x(5), x(5), x(6)/4*1e-3, x(6)/4*1e-3, (x(6)/4+x(7)/2)*1e-3, (x(6)/4+x(7)/2)*1e-3];
%plsplot(pinf, {struct('read','@dbzprep','prep','@dbzprep','pi','@dbzpi'),'right'});
plsreg(pinf, plsnum);
%% 98 Ramsey with generalized evo and controllable gate voltages 2 flips.
plsnum = 98;

clear pinf;
pinf.name='Ramsey';
%     1         2         3        4       5       6     7      8      9       10
els={'@start','@fill','@wait','@reload','@wait','@prep',};
for j = 1:4
    els{end+1} = '@exch';
end
els = [els,{'@read','@wait', '@meas'}];
pinf.data=struct('type',els,'time',[],'val',[]);
% Parameters: p(1) = pulse length, eps_x1, eps_y1, eps_x2, epsy2, time 
pinf.pardef = [2 -1;];
pinf.pardef = zeros(3*4+1,2);
pinf.pardef(1,:) = [2 -1];
pinf.pardef(2:3,:) = [(7:2:10)',ones(2,1)];
pinf.pardef((2:3)+2,:) = [(7:2:10)',2*ones(2,1)];
pinf.pardef((2:3)+2*2,:) = [(8:2:10)',ones(2,1)];
pinf.pardef((2:3)+2*3,:) = [(8:2:10)',2*ones(2,1)];
pinf.pardef((2:5)+2*4,:) = [(7:10)',-1*ones(4,1)];

pinf.trafofn = @(x) [x(1), x(2)*ones(2,1)',x(3)*ones(2,1)',x(4)*ones(2,1)',x(5)*ones(2,1)',1e-3*(1/4)*x(6)*ones(4,1)',];

plsreg(pinf, plsnum);
%% Pulsed Zoom
n = length(plsdata.pulses)+1; 
clear pls 
pls = plsdata.pulses(5); % start by copying old zoom pulse
pls.name = 'PulsedZoom'; 
pls.pardef= [6,-1; 7,2; 7,3];
pls.trafofn=@(x)[x(1), x(3),x(2)];
plsreg(pls,n); 
plssync('save'); 
%% 116 wacky marker wait
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

plsdata.pulses(116) = pls; 
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
%% 117 Ramsey with marker burst
n = 117; 
clear pls 
pls = plsdata.pulses(15); 
pls.name = 'RamseyMarker'; 
pls.data(5).type = 'markerburst';
pls.data(5).val = [NaN, NaN, NaN, 1 0 0 0]; 
%           PlsLen, 3 exch vals, exch time, marker pre.
pls.pardef= [9,-1; 5,1; 5,2; 5,3; 5,-1; 5 -2; 5 -3; 7 -1]; 
pls.trafofn=@(x)[x(1), x(2),x(3),x(4),x(7)*1e-3,x(5)*1e-3,x(6)*1e-3,.4-x(7)*1e-3];
plsreg(pls,n); 
plssync('save');
%% 118 Ramsey with marker 
n = 118; 
clear pls 

pls.name = 'RamseyMarker2'; 
%     1         2        3        4       5         6           7        8       9      10
els={'@start','@reload','@wait','@prep','@markerburst','@read','@wait', '@meas','@fill','@wait'};
pls.data=struct('type',els,'time',[],'val',[]);

pls.data(5).val = [NaN, NaN, NaN, 0 1 0 0]; 
pls.pardef=     [9,-1; 5,1; 5,2; 5,3; 5,-1; 5 -2; 5 -3; 7 -1; 3 -1; 8 -1]; 
%           PlsLen, 3 exch vals, exch time, marker pre, marker delay, wait time. 
pls.trafofn=@(x)[x(1), x(2),x(3),x(4),x(10)*1e-3,x(5)*1e-3,x(6)*1e-3,x(7)-x(10)*1e-3,x(8),x(9)];
plsreg(pls,n); 
plssync('save');
%% 119
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
%% 120 2D Load 
%n = length(plsdata.pulses)+1; 
n = 120;
clear pls 
pls = plsdata.pulses(6); % start by copying old load pulse
pls.data(end).time = [];
pls.name = 'TwoDLoad'; 
pls.pardef= [7 -1;6,-1; 6,-2; 6,1;6,2];
pls.trafofn = @(x)[x(1),x(2)*1e-3,x(3)*1e-3,x(4:5)+[x(6),x(7)]];
plsreg(pls,n); 
plssync('save');
%% 121 adprep and then either adread or read.
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
%% 131 Echo with marker, not symmetric
clear pls 
plsnum = 131;
pls.name='RamseyEMarker';

%     1        2       3        4        5       6       7      8       9       10       11      12
els={'@start','@fill','@wait','@reload','@wait','@prep','@evo','@dbzpi','@exch','@read','@wait','@measLoc'};
pls.data=struct('type',els,'time',[],'val',[]);

pls.pardef = [2,-1; 7,3; 9,3; 7,4; 7,5; 7,-2; 7,-3; 11 -1; 7 -1; 9, -1];
%           PlsLen, exch,exch, I, Q, marker pre/post echo time, wait time, offset, 
pls.trafofn=@(x)[x(1),x(2),x(2),x(3),x(4),x(5)*1e-3,x(6)*1e-3,x(8)*1e-3, x(7)/2, x(7)/2+x(9)*1e-3];
plsreg(pls,plsnum); 
plssync('save');
%% 132 IQ burst start. 
clear pls 
plsnum = 132;
pls.name='IQtoy';

%     1        2         3       4       5          6       7       8          9      10       
els={'@start','@reload','@wait','@prep','@IQburst','@read','@wait','@measLoc','@fill','@wait'};
pls.data=struct('type',els,'time',[],'val',[]);

pls.pardef=  [9,-1; 3,-1; 5,3; 5 4; 5 5; 5 -2; 5 -3; 5 -1]; 

%             PlsLen, wait time, exch, I, Q,   wait, marker pre, marker delay, time
pls.trafofn=@(x)[x(1),x(2),x(3),x(4),x(5),x(6)*1e-3,x(7)*1e-3,x(8)*1e-3];
plsreg(pls,plsnum); 
plssync('save');
%% 133 IQ Echo, symmetric, fully customizable 
clear pls 
plsnum = 133;
pls.name='RamseyESym';

%     1        2       3        4        5       6       7      8       9       10       11      12
els={'@start','@fill','@wait','@reload','@wait','@prep','@evo','@dbzpi','@evo','@read','@wait','@measLoc'};
pls.data=struct('type',els,'time',[],'val',[]);
%             1    2    3    4    5    6    7     8    9    10   11   12      13      14
pls.pardef = [2,-1; 7,3; 9,3; 7,4; 9,4; 7,5; 9,5; 7,-2; 9,-2; 7,-3; 9,-3; 11 -1; 7 -1; 9, -1];
%           PlsLen, exch,exch, I, Q, marker pre/post echo time, wait time, offset, 
pls.trafofn=@(x)[x(1),x(2),x(2),x(3),x(3),x(4),x(4),x(5)*1e-3,x(5)*1e-3,x(6)*1e-3,x(6)*1e-3,x(8)*1e-3, x(7)/2, x(7)/2+x(9)*1e-3];
plsreg(pls,plsnum); 
plssync('save');
%% 134 dBz pulse 
plsnum = 134;
                            %1          2       3         4        5
pinf.data=struct('type',{'@start', '@wait',  '@reload', '@sep', '@meas','fill','@wait'},...
                 'time',{  []        ,[]  ,   []      ,  []   , []     ,  []  , []},...
                 'val' ,{  []        , [] ,   []      ,  []   ,[]      ,  []  , []});
%params = [pulse length, max sep time, sep time]
pinf.pardef = [2 -1; 4 -1; 6 -1];
pinf.trafofn = @(x)[x(1), (x(2)-x(3))*1e-3, x(3)*1e-3];
%plsplot(pinf, 'right')
plsreg(pinf, plsnum);
plssync('save');
%% 135 Load, ramp, exchange, ramp,read 
plsnum = 135;

clear pinf;
pinf.name='STP_ad';
%     1         2         3        4     5      6      7      8       9     10
els={'@start','@fill','@wait','@load','@exch','ramp','@exch','@exch','ramp','@read' '@meas','@wait'};
pinf.data=struct('type',els,'time',[],'val',[]);
% Parameters: [pauseTime,start1, end1, time1, start2, end2, time2, wait_loc, pls_length, wait_time]
pinf.pardef = [2 -1; 5 3; 6 3; 6 -1; 8 3; 9 3; 9 -1; 7 3; 7 -1; 6 1; 6 2; 9 1; 9 2; 5 -1; 8 -1;];

pinf.trafofn = @(x) [x(9), x(2), x(3), x(4), x(5), x(6), x(7),x(8),x(10), -1, 1 ,-1,1,x(1),x(1)];

plsreg(pinf, plsnum);
%% 136 variable length everything off pulse

clear pinf;
pinf.name='All_off_var'
plsnum = 136;
els = {'@start', '@fill', '@wait'};
pinf.data= struct('type', els, 'time', [], 'val', []);
pinf.pardef=[2 -1];
pinf.trafofn=@(x)[x(1)];
plsreg(pinf, plsnum);
%% 137 Singlet feedback polarizer, dict based, brave new world feedback, 500 ns
% dict-based feedback elements.
% expects supplemental dictionary to define @stpsweep
plsnum = 137;
clear pinf;
pinf.name = 'Spol_500';
%       1          2       3         4       5
els ={'@start', '@fill', '@reload', '@sep', '@stpsweep'};
pinf.data = struct('type', els, 'time', [], 'val', []);
pinf.data(2).time = .5; %500ns pulse length
pinf.data(4).time = .001; %1ns loiter at far eps
pinf.data(3).time = [.005 .08 .005]; % hard code short load; this will get filled out to longer.
pinf.pardef = [4 -1];
pinf.xval=1;
%pinf.trafofn=@(x) [x(1)*1e-3];
%Parameters; sep time in ns
plsplot(pinf, {'right'});
plsreg(pinf,plsnum);
%% 138 Triplet unconditional polarizer, dict based, brave new world feedback, 500 ns
% dict-based feedback elements.
% expects supplemental dictionary to define @stpsweep, @tlsweep
% ignorse supplied parameters.
plsnum = 138;
clear pinf;
pinf.name = 'Tpol_500';
%       1          2       3         4       5
els ={'@start', '@fill', '@wait', '@tlsweep','@stpsweep'};
pinf.data = struct('type', els, 'time', [], 'val', []);
pinf.data(2).time = .5; %500ns pulse length
pinf.data(3).time=1e-3;
pinf.pardef = [];
pinf.trafofn=@(x) [];
pinf.xval=0;
%Parameters; sep time in ns
plsplot(pinf, {'right'});
plsreg(pinf,plsnum);
%% 139 Filler pulse for feedback, 500 ns
clear pinf;

pinf.data.pulsetab = [0 .5;0 0;0 0];
pinf.name = sprintf('fill_%d', round(plsdata.tbase * pinf.data.pulsetab(1, 2)));
plsreg(pinf, 139);
%% 140 Singlet feedback polarizer, dict based, brave new world feedback, 750 ns
% dict-based feedback elements.
% expects supplemental dictionary to define @stpsweep
plsnum = 140;
clear pinf;
pinf.name = 'Spol_750';
%       1          2       3         4       5
els ={'@start', '@fill', '@reload', '@sep', '@stpsweep'};
pinf.data = struct('type', els, 'time', [], 'val', []);
pinf.data(2).time = .75; %500ns pulse length
pinf.data(4).time = .001; %1ns loiter at far eps
pinf.data(3).time = [.005 .08 .005]; % hard code short load; this will get filled out to longer.
pinf.pardef = [4 -1];
pinf.xval=1;
%pinf.trafofn=@(x) [x(1)*1e-3];
%Parameters; sep time in ns
%plsplot(pinf, {'right'});
plsreg(pinf,plsnum);
plssync('save'); 
%% 141 Triplet unconditional polarizer, dict based, brave new world feedback, 750 ns
% dict-based feedback elements.
% expects supplemental dictionary to define @stpsweep, @tlsweep
% ignorse supplied parameters.
plsnum = 141;
clear pinf;
pinf.name = 'Tpol7500';
%       1          2       3         4       5
els ={'@start', '@fill', '@wait', '@tlsweep','@stpsweep'};
pinf.data = struct('type', els, 'time', [], 'val', []);
pinf.data(2).time = .75; %500ns pulse length
pinf.data(3).time=1e-3;
pinf.pardef = [];
pinf.trafofn=@(x) [];
pinf.xval=0;
%Parameters; sep time in ns
%plsplot(pinf, {'right'});
plsreg(pinf,plsnum);
plssync('save'); 
%% 142  Filler pulse for feedback, 750 ns
clear pinf;

pinf.data.pulsetab = [0 .75;0 0;0 0];
pinf.name = sprintf('fill_%d', round(plsdata.tbase * pinf.data.pulsetab(1, 2)));
plsreg(pinf, 142);
plssync('save');
%% 143 Fully configurable dBz pulse, fill at start
% looks exactly like pulse 12, except the fill is at the start

clear pinf;
plsnum = 143;
                           %1         2       3         4        5     6
pinf.data=struct('type',{'@start', 'fill','@wait', '@reload', '@sep', '@meas','@wait'},...
                 'time',{  []        ,[]  ,   []      ,  []   , []     ,  [] ,[]},...
                 'val' ,{  []        , [] ,   []      ,  []   ,[]      ,  [] ,[]});

%params = [pulse length, max sep time, readdout time, sep time]
pinf.pardef = [2 -1; 5 -1; 6 -1; 2 -1];
pinf.trafofn = @(x)[(x(2)-x(4))*1e-3, x(4)*1e-3, x(3), x(1)];
%plsplot(pinf, 'right')
plsreg(pinf, plsnum);
plssync('save'); 
%% 144 rabi in the rotating frame, markerburst, eps correction based on time.
% expects prep2 to have the form of markerburst
plsnum = 144;

clear pinf;
pinf.name='RabiVCO';
%     1         2        3        4         5       6       7        8        9
els={'@start','@fill','@wait', '@reload','@wait','@prep','@prep2','@read','@meas','@wait'};
pinf.data=struct('type',els,'time',[],'val',[]);

%               1             2          3             4
% Parameters: [pulselength, wait_eps, slope(mV/us), wait_time(ns)] 
pinf.pardef = [2 -1; 7 3; 7 -1]; 

pinf.trafofn = @(x) [x(1), x(2)+x(3)*x(4)*1e-3, x(4)*1e-3];
%plsplot(pinf, {struct('read','@adread','prep','@adprep'),'left'});
plsreg(pinf, plsnum);
plssync('save'); 
%% 145 CondEvo with marker, not symmetric
clear pls 
plsnum = 145;
pls.name='CondEvoIQ';

%     1        2       3        4        5       6       7      8       9       10       11      12
els={'@start','@fill','@wait','@reload','@wait','@sep','@evo','@pi','@exch','@read','@wait','@measLoc'};
pls.data=struct('type',els,'time',[],'val',[]);
            % 1     2    3    4    5    6     7     8      9     10    11    
pls.pardef = [2,-1; 7,3; 9,3; 7,4; 7,5; 7,-2; 7,-3; 6 -1; 5 -1; 7 -1; 9, -1];
%             1 PlsLen, 2 exch, 3 I 4 Q, 5/6 mk pre/post, 8 echo time, 9 dBz time, 10 offset
pls.trafofn=@(x)[x(1),x(2),x(2),x(3),x(4),x(5)*1e-3,x(6)*1e-3,x(8)*1e-3,x(9)*1e-3, x(7)/2, x(7)/2+x(10)*1e-3];
plsreg(pls,plsnum); 
plssync('save');
%% 146, Ramsey assymmetric echo with the wait moved 
clear pls 
plsnum = 146;
pls.name='RamseyEMarker2';

%     1        2       3        4        5       6       7      8       9       10       11      12
els={'@start','@fill','@wait','@reload','@wait','@prep','@evo','@dbzpi','@exch','@read','@wait','@measLoc'};
pls.data=struct('type',els,'time',[],'val',[]);

pls.pardef = [2,-1; 7,3; 9,3; 7,4; 7,5; 7,-2; 7,-3; 5 -1; 7 -1; 9, -1];
%           PlsLen, exch,exch, I, Q, marker pre/post echo time, wait time, offset, 
pls.trafofn=@(x)[x(1),x(2),x(2),x(3),x(4),x(5)*1e-3,x(6)*1e-3,x(8)*1e-3, x(7)/2, x(7)/2+x(9)*1e-3];
plsreg(pls,plsnum); 
plssync('save');