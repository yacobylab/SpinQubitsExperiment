%% Change sides
autotune.swap('right'); 
%% Examine autotune pulses. 
atplschk('loadTime_1_R','clf',{'offset',-tuneData.measPt})

atplschk('loadPos_1_R','clf',{'offset',-tuneData.measPt,'pulses',[1 50 100]})
atplschk('STP_1_R','clf',{'offset',-tuneData.measPt,'pulses',[1 50 100]})
atplschk('topLead_2_R','clf',{'offset',-tuneData.measPt,'pulses',[1 50 100]})
%%
%% Show current dictionary values 
printDict('right'); 
%% Manually change pulse 
r=pdload('right'); 
r.reload.val = [-2.5,-11]; 
pdsave('right',r); 
tuneData.loadPos.updateGroup; 
%% This will update all pulses to have the load position of the minimum value for load Pos scan. 
tuneData.loadPos.updateGroup('target'); 
%% 
tuneData.updateGroups('both'); 
%% Choose direction of 1,1 vs. 0,3 
tuneData.sepDir = [-1,1]; 
%% Print groups on awg 
awggroups
%% Things to run without measurement point
tuneData.lead.run
tuneData.line.run
tuneData.zoom.run; 
%% For wider junction 
tuneData.lead.run
tuneData.line.run
tuneData.zoom.run('wide'); 
%% With meas pt
tuneData.loadPos.run;  % This doesn't need load point
% These do: 
tuneData.loadTime.run; 
tuneData.stp.run; 
tuneData.tl.run
%% Move scans
tuneData.twoSen.scan.loops.rng = [];
tuneData.sensor.scan.loops.rng = tuneData.twoSen.scan.loops.rng;
tuneData.chrg.loops(2).settle = 0.25; 

%% Make basis 
atxyfixAll('right','chrg')
atgradfix('all',2e-3,'right')
atgradfix('Lead3 VRes',-2e-3,'right') % Just do a couple of gates. 