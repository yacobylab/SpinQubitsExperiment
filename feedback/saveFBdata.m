function scan = saveFBdata(scan)
% helper function for saving FBdata with scan. 
% function scan = saveFBdata(scan)
global fbdata; 

scan.data.FBset = fbdata.set;
scan.data.pumpHist = fbdata.pumpHist; 
scan.data.gradHist = fbdata.gradHist; 
scan.data.err = fbdata.err; 
fbdata.set = []; fbdata.gradHist={}; fbdata.pumpHist = {}; fbdata.err = {};

scan.data.FBData = fbdata; 
end