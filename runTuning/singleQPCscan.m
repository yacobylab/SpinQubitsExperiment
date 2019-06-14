%% Make scan 

smset('LockinTauA',0.03);
lockin = 'LockinA'; 

scan.loops(1).ramptime = .06;
scan.loops.rng = [-.0 -.91];
scan.loops(1).npoints= 100; 

scan.loops(1).getchan=lockin;
scan.saveloop= [2 1]; 
scan.disp(1).channel = 1; scan.disp(1).dim = 1; scan.disp(1).loop=1;  
scan.loops(1).testfn.fn = @atZero; scan.loops(1).testfn.args = {1e-10,10};
%scan.disp(2).channel = 2; scan.disp(2).dim = 1; scan.disp(2).loop=1;  
%% Close off extra ohmic channel
smset({'3b','4b','T34'},-.85); 
%% 
smset(1:20,0); 
smset({'1b','2b','T12'},-.85)
%% Run scan 
currSetChans = {'1b','T12'}; 
scan.loops(1).setchan = currSetChans; 
scanName=smnext(sprintf('qpc_4k_%s_%s',currSetChans{1},currSetChans{2}));
%scanName=smnext(sprintf('qpc_4k_%s_%s',currSetChans{1}));
scan.comment = '';
smrun(scan,scanName);
smset(currSetChans,0);
%% Scan close val
currSetChans = {'SD4top'}; 
scan.loops(1).setchan = currSetChans; 
scanName=smnext('qpc_leftClose');
%scanName=smnext(sprintf('qpc_4k_%s',currSetChans{1}));
%scanName=smnext(sprintf('qpc_4k_%s_%s',currSetChans{1}));
smrun(scan,scanName);

%%
scanCount = scan; 
scanCount.loops(1).npoints = 1000; 
scanCount.loops(1).setchan = {'count'}; 
scanCount.loops(1).ramptime = 0.001; 
smrun(scanCount,smnext('qpc_testNoise'))
%% yoko scan
scanYoko = scanCount; 
scanCount.loops(1).npoints = 250;
smset('LockinTauB',0.003)
scanCount.loops(1).setchan = {'yoko'}; 
scanCount.loops(1).rng = [0 -1]; 
scanCount.loops(1).ramptime = 0.001; 
smrun(scanCount,smnext('qpc_yoko_SD1bot_1b'))

%%

%smset({'1b','T12','2b'},-0.85)
%smset({'SD4top','4b','SD4bot','4a','T34','3b'},0)
scandac = scan; 
scandac.disp(1).channel = 1; scandac.disp(1).dim = 1; scandac.disp(1).loop=2;  
scandac.configfn.fn = @smabufconfig2; 
scandac.configfn.args = {'trig arm'}; 
currSetChans = {'SD4top','4b'}; 
scandac.loops(1).testfn =[];

scandac.loops(2).getchan = {'QPCbufA'}; 
scandac.loops(1).getchan = [];
scandac.loops(1).ramptime = -1/16; 
scandac.loops(2).npoints =1; 
scandac.loops(1).npoints = ceil(2350/16); 
scandac.loops(1).rng = [0 -.91]; 
%profile on;  

scandac.loops(1).setchan = currSetChans; 
scanName=smnext(sprintf('qpc_4k_%s_%s',currSetChans{1},currSetChans{2}));
smrun(scandac,scanName); 
%scandac.loops(1).setchan = currSetChans; 
%scanName=smnext(sprintf('qpc_4k_%s_%s',currSetChans{1},currSetChans{2}));
%profile viewer; 
smset(scandac.loops(1).setchan,0); 

%smrun(scandac,smnext('qpc_SD1bot_1b'))
%%

smset({'SD1top','1b','SD1bot','1a','T12','2b'},0)
smset({'T34','3b'},-.9)
%smset({'SD4top','SD4bot','4a','T34','2b','T34','3b','4b','1b'},-.9);
%smset({'SD4top','SD4bot','4b','4a'},0);
scandac = scanYoko; 
scandac.disp(1).channel = 1; scandac.disp(1).dim = 1; scandac.disp(1).loop=1;  
%scandac.configfn.fn = @smabufconfig2; 
%scandac.configfn.args = {'trig arm'}; 
currSetChans = {'SD1bot','1a'}; 
scandac.loops(1).testfn =[];

scandac.loops(1).getchan = 'LockinA';
scandac.loops(1).ramptime = 0.06; 
scandac.loops(1).npoints =100; 
scandac.loops(1).rng = [0 -.91]; 
%profile on;  

smrun(scandac,scanName); 
%profile viewer; 
smset(scandac.loops(1).setchan,0); 

%smrun(scandac,smnext('qpc_SD1bot_1b'))