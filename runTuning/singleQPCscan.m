%% Make scan 
smset('LockinTauA',0.1);
lockin = 'LockinA'; 

scan.loops(1).ramptime = .2;
scan.loops.rng = [0 -1.9];
scan.loops(1).npoints= 100; 

scan.loops(1).getchan=lockin;
scan.saveloop= [2 1]; 
scan.disp(1).channel = 1; scan.disp(1).dim = 1; scan.disp(1).loop=1;  
scan.loops(1).testfn.fn = @atZero; scan.loops(1).testfn.args = {1e-10,10};
scan.disp(2).channel = 2; scan.disp(2).dim = 1; scan.disp(2).loop=1;  
%% Close off extra ohmic channel
smset({'3b','4b','T34'},-1.5); 
%% Run scan 
currSetChans = {'T12','Res'}; 
scan.loops(1).setchan = currSetChans; 
scanName=smnext(sprintf('qpc_4k_%s_%s',currSetChans{1},currSetChans{2}));
%scanName=smnext(sprintf('qpc_4k_%s_%s',currSetChans{1}));
smrun(scan,scanName);
smset(currSetChans,0); 