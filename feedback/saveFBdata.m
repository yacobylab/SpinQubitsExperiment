function scan = saveFBdata(scan)
global fbdata; 

scan.data.FBset = fbdata.set;
scan.data.pumpHist = fbdata.pumpHist; 
scan.data.gradHist = fbdata.gradHist; 

fbdata.set = []; fbdata.gradHist={}; fbdata.pumpHist = {}; 
end