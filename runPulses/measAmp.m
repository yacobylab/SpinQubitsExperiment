function scan = measAmp(scan,opts)
% Set up scan to measure in amp mode by setting the PlsRamps as consts. 
global tuneData
if ~exist('opts','var'), opts = ''; end
scan.consts(end+1).setchan=tuneData.xyChan{1};
scan.consts(end).val=tuneData.measPt(1);
scan.consts(end+1).setchan=tuneData.xyChan{2};
scan.consts(end).val=tuneData.measPt(2);
if isopt(opts,'both')
    autotune.swap('next'); 
    scan = measAmp(scan);     
end
end