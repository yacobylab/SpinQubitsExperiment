% Notes and def'ns 
% Green lines show range of scan. Starts at white dot and ends at black
% dot. 
% default offset is the cetner of zoom scan, defined as the 
%% Change sides
autotune.swap('right'); 
%% Examine autotune pulses. 
atplschk('loadTime_1_R','clf',{'offset',-tuneData.measPt})
atplschk('loadPos_1_R','clf',{'offset',-tuneData.measPt,'pulses',[1 50 100]})
atplschk('STP_1_R','clf',{'offset',-tuneData.measPt,'pulses',[1 50 100]})
atplschk('topLead_2_R','clf',{'offset',-tuneData.measPt,'pulses',[1 50 100]})
atplschk('dBz_swfb_128_R','clf',{'offset',-tuneData.measPt,'pulses',[1 50 100]})
%% Show current dictionary values 
printDict('right'); 
%% Manually change pulse 
r=pdload('right'); 
r.reload.val = [-2,-15]; 
pdsave('right',r); 
tuneData.loadPos.updateGroup; 
%% 
tuneData.updateGroups('both'); 
%% Choose direction of 1,1 vs. 0,3 
tuneData.sepDir = [-1,1]; 
%% Print groups on awg 
awggroups
%% This will update all pulses to have the load position of the minimum value for load Pos scan. 
tuneData.loadPos.updateGroup('target'); 
%% Retune sensore and run charge scan
tuneData.sensor.run; 
tuneData.chrg.run
%% Things to run without measurement point
tuneData.lead.run
tuneData.line.run
tuneData.zoom.run; 
%% For wider junction 
tuneData.lead.run
tuneData.line.run
tuneData.zoom.run('wide'); 
%% Choose new measurement point 
tuneData.zoom.ana('man last'); 
%% With meas pt
tuneData.loadPos.run;  % This doesn't need load point
tuneData.loadTime.run; 
tuneData.stp.run; 
tuneData.tl.run
%% Move scans
tuneData.twoSen.scan.loops.rng = [];
tuneData.sensor.scan.loops.rng = tuneData.twoSen.scan.loops.rng;
tuneData.chrg.loops(2).settle = 0.25; 
%% Change stp target 
tuneData.stp.target=800;
tuneData.stp.updateGroup;
tuneData.stp.run;
%% Make basis 
atxyfixAll('right','chrg')
atgradfix('all',-2e-3,'right')
%atgradfix('Lead3 VRes',-2e-3,'right') % Just do a couple of gates. 
%% Run t1 scan without working gradient. 
tuneData.t1.run('nograd') 
%% Remove set of groups 
awgrm(13,'after'); 
awgclear('unused'); 
%% Pulsed zoom scan
tuneData.zoom.pulsed('wide') %
%% Center load scan certain distance down loead 
tuneData.loadPos.dist = 5; % dist from triple point 
%tuneData.loadPos.dist = -7; % slope of lead
tuneData.loadPos.updateGroup('init'); 
%% Center dot 
tuneData.center