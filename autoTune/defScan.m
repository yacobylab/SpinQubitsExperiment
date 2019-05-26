function scan = defScan(opts,side)
% Create simple 2D scan with DAQ sensing. 
% function scan = defScan
% chrg: autoscan type scans, not using PlsRamps. Otherwise autotune scan. 
% side: left or right. 
if ~exist('opts','var'), opts = ''; end

scan.saveloop = [2 1];
scan.disp(1) = struct('loop',2,'channel',1,'dim',1);
scan.disp(2) = struct('loop',2,'channel',1,'dim',2);

scan.configfn.fn = @smabufconfig2;
scan.configfn.args = {'arm',1};
if ~isopt(opts,'chrg')
    scan.cleanupfn.fn = @smaconfigwrap;
    scan.cleanupfn.args = {@smset,{'PlsRamp1','PlsRamp2','PlsRamp3','PlsRamp4'},0};
end
scan.loops(1).trigfn.fn = @smatrigAWG; 
scan.loops(1).trigfn.args = {'AWG1'}; 
scan.loops(1).stream = 1; scan.loops(1).waittime=-1; 
scan.consts(1) = struct('setchan','samprate','val',4e7); 
if isopt(opts,'chrg') 
    if ~exist('side','var'), side = 'left'; end    
    scan.configfn(end+1).fn = @pulseLineWrap;
    if strcmp(side,'left')        
        scan.configfn(end).args = {'chrg_1_L'};
        scan.loops(2).getchan = 'DAQ1'; 
    else
        scan.configfn(end).args = {'chrg_1_R'};
        scan.loops(2).getchan = 'DAQ2'; 
    end
else
    scan.consts(2).setchan = 'PulseLine'; 
end
end