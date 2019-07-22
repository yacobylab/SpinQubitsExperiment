% Notes and def'ns 
% Green lines show range of scan. Starts at white dot and ends at black
% dot. 
% default offset is the cetner of zoom scan, defined as the 
%% Change sides
autotune.swap('right'); 
%% Examine autotune pulses. 
atplschk('zoom_1_R','clf',{'offset',-tuneData.measPt,'pulses',[1 2]})
atplschk('loadTime_1_R','clf',{'offset',-tuneData.measPt})
atplschk('loadPos_1_R','clf',{'offset',-tuneData.measPt,'pulses',[1 50 100]})
atplschk('STP_1_R','clf',{'offset',-tuneData.measPt,'pulses',[1 50 100]})
atplschk('topLead_2_R','clf',{'offset',-tuneData.measPt,'pulses',[1 50 100]})
atplschk('dBz_swfb_128_R','clf',{'offset',-tuneData.measPt,'pulses',[1 50 100]})
%% Show current dictionary values 
printDict;
%% Manually change load pulse (for instance when tuning up new device)
r=pdload('right'); 
r.reload.val = [-2,-15]; 
pdsave('right',r); 
%tuneData.loadPos.updateGroup; 
tuneData.updateAll('nodict'); 
%% Add groups on both sides
tuneData.updateGroups('both'); 
%% Create new charge template (for auto triple point identification) 
chrgTemplate; 
%% Choose direction of 1,1 vs. 0,3 
tuneData.sepDir = [-1,1]; 
%% Print groups on awg 
awggroups
%% sensor dot scan
tuneData.twoSen.run
%% Retune sensor and run charge scan
tuneData.sensor.run('fine'); 
tuneData.chrg.run
%% Reanalyze charge scan 
tuneData.chrg.ana('man mnslp last'); 
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
%% Load scan (need meas pt)
tuneData.loadPos.run; % This doesn't need load point
%% This will update all pulses to have the load position of the minimum value for load Pos scan. 
tuneData.loadPos.updateGroup('target');
%% New load pulse 
%tuneData.loadPos.slope = -2.2; 
tuneData.loadPos.dist = 3; 
tuneData.loadPos.rangeScale = 4;
tuneData.loadPos.updateGroup('init'); 
tuneData.loadPos.run
%% Autotune pulses with readout. 
tuneData.loadTime.run; 
tuneData.stp.run; 
tuneData.tl.run
%% rundBz, no feedback
rundBz
%% check Sep amp. 
testSep
%% Run t1 scan without working gradient. 
tuneData.t1.run('nograd') 
%% Center dot 
tuneData.center
%% Move scans
tuneData.twoSen.scan.loops(1).rng = [-.4 -.55];
tuneData.twoSen.scan.loops(2).rng = [-.38 -.28];
tuneData.sensor.scan.loops(1).rng = tuneData.twoSen.scan.loops(1).rng;
%% Change stp target 
%tuneData.stp.target=00;
tuneData.stp.updateGroup;
tuneData.stp.run;
%% Change tl point
tuneData.tl.target=0;
tuneData.tl.dist = 2; 
tuneData.tl.updateGroup;
tuneData.tl.run;
%% Remake tl group to fit point (only use if fit worked)
tuneData.tl.updateGroup('target'); 
tuneData.stp.updateGroup('target'); 
tuneData.tl.run
tuneData.stp.run
%% Make basis 
atxyfixAll('right','chrg')
atgradfix('all',-5e-3,'right')
%%
atxyfixAll('left');         
%atgradfix('Lead3 VRes',-2e-3,'right') % Just do a couple of gates. 
%% Remove set of groups 
awgrm(13,'after'); 
awgclear('unused'); 
%% Pulsed zoom scan
tuneData.zoom.pulsed('wide') %
%% Center load scan certain distance down loead 
%tuneData.loadPos.dist = 2; % dist from triple point 
%tuneData.loadPos.dist = -7; % slope of lead
tuneData.loadPos.updateGroup('init'); 
%%
smset('Bz',0.7)
smaMercury3axis('heaterOff')
%% Check phase
scandata.autoramp = 0; 
autoscan('RF'); 
scandata.autoramp = 1; 
%% Change charge scan
tuneData.chrg.scan.loops(2).rng = [-0.007 0.007];
tuneData.chrg.scan.loops(1).rng = [-0.007 0.007];