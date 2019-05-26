%% Running quantum dot / s-t0 qubit experiments 
%% Starting computer: 
startmx
%% Run QPC scans 
measureQPCauto('SL','quiet')
measureQPCauto('SR','quiet')
measureQPCauto('QL','quiet')
measureQPCauto('QR','quiet')
%% Move scan to the left and down by 20 mV
autoscan('last',{'move2',-0.02,'move1',-0.02})
autoscan('last'); 
%% Zoom in on part of previous scan you want to scan next, then run this to create new scan.
autoscan('last','gca') 
%% Change scan direction easily
autoscan('last','up') % Run scan from minimum value up (down is option for opposite)
%% Set up scan with a given point spacing
autoscan('last','gofast',{'pointSpacing1',1e-3});
%% Change to charge sensing with lock in 
autoscan('','startSens') 
% Copies ranges but reduces point spacing to 1 mV in X. 
%% Change gate used in sensor scan 
autoscan('',{'sensorGate','SD4bot'}); 
%% Set up DAQ
autoscan('','startDAQ'); 
autoscan('RF'); 
% Will set to minimum value. 
%% Charge Sensing DAQ
autoscan('SD');
autoscan('sensor'); 
autoscan('sens','trafa'); 
autoscan('sens'); 
%% New approach to starting charge sensing... 
% First, analyze data so we can get diff. Then, run set function click on the point you want
% to use. 
plotChrg('simple noppt'); % Simple gets rid of button down functions. 
center = autoscan('','set'); 
% Check what value trafofn had at center value. 
sensorVal = autoscan('',{'sensorVal',center}); 
smset(scandata.sens.loops(2).setchan{2},sensorVal); 
%% Reset to the center values. 
smset(scandata.sens.loops(1).setchan{1},center(1)); 
smset(scandata.sens.loops(2).setchan{1},center(2));
%% Set up junction scan. Trafofn will be automatically copied from sens. 
%center = [-.8, -.9]; % Here, fill in value you want to center scan at. 
autoscan('juncd',{'center',center,'diff1',0.03,'diff2',0.03}); 
autoscan('juncd'); 
%% Junc scan with no trafofn
%center = [-.8, -.9]; % Here, fill in value you want to center scan at. 
scandata.juncd.loops(2).trafofn=[];
scandata.autoramp=0; 
autoscan('juncd',{'center',center,'diff1',0.03,'diff2',0.03}); 
autoscan('juncd'); 
%% Change to charge scan. 
autoscan('','set'); % Click on value for the junction you are using. 
%% Go back to gate values for better tuning set: 
smrestore 
%% Print list of channels associated with instrument. 
