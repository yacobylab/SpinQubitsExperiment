function scan = defScan
scan.saveloop = [2 1];
scan.disp(1) = struct('loop',2,'channel',1,'dim',1);
scan.disp(1) = struct('loop',2,'channel',1,'dim',2);

scan.configfn.fn = @smabufconfig2;
scan.configfn.args = {'arm',1};

scan.cleanupfn.fn = @smaconfigwrap;
scan.cleanupfn.args = {@smset,{'PlsRamp1','PlsRamp2','PlsRamp3','PlsRamp4'},0};

scan.loops(1).trigfn.fn = @smatrigAWG; 
scan.loops(1).trigfn.args = {'AWG1'}; 

scan.consts(1) = struct('setchan','samprate','val',4e7); 
scan.consts(2).setchan = 'PulseLine'; 

scan.disp(1).channel = 1; 
scan.disp(1).dim = 1; 
scan.disp(1).loop = 2; 

scan.disp(2).channel = 1; 
scan.disp(2).dim = 2; 
scan.disp(2).loop = 2; 
end