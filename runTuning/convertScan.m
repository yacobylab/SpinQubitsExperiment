function scanChan = convertScan(scanChan,scanRng) 
% convert between DAQ and lock in scans. 
% function scanChan = convertScan(scanChan,scanRng) 

global smdata; 
    for i = 1:length(scanRng.loops) 
        scanChan.loops(i).rng = scanRng.loops(i).rng; 
        scanChan.loops(i).setchan = scanRng.loops(i).setchan; 
        if isfield(scanRng.loops,'trafofn') && ~isempty(scanRng.loops(i).trafofn)
            scanChan.loops(i).trafofn = scanRng.loops(i).trafofn; 
        end
    end
    getchan = [scanChan.loops.getchan]; 
    instchan = smdata.channels(chl(getchan)).instchan;  
    if isfield(smdata.inst(instchan(1)).data,'Tau') 
        if smdata.inst(instchan(1)).data.Tau* 2 > abs(scanChan.loops(1).ramptime) 
            scanChan.loops(1).ramptime = -smdata.inst(instchan(1)).data.Tau* 2;
        end
    end
end