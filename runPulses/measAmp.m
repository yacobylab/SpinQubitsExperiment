function scan = measAmp(scan)
global tuneData
scan.consts(end+1).setchan=tuneData.xyChan{1};
scan.consts(end).val=tuneData.measPt(1);
scan.consts(end+1).setchan=tuneData.xyChan{2};
scan.consts(end).val=tuneData.measPt(2);
end