function runTimerScans(~,~)
persistent scanInd; global scanInfo;
if isempty(scanInd)
    scanInd = 1;
end
chanVals = scanInfo.gateValsMat(:,scanInd);
smset(scanInfo.channels,chanVals);
scanInd = scanInd +1;
if scanInfo.runSensor
    autoscan('sensor'); autoscan('sens','trafa')
end
data = autoscan('sens');
if mean(nanstd(data{1},[],2))<0.01
    scanInfo.runSensor = true;
else
    scanInfo.runSensor = false;
end
end