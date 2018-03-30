% Assume we are in current data directory 
cd .. 
newDir =['data_' datestr(now,'yyyy_mm_dd')];
mkdir(newDir); 
cd(newDir) 
scandata.file = ['Z:\qDots\data\mx50structs\scandata_' datestr(now,'yyyy_mm_dd')]; 
smdata.file = ['Z:\qDots\data\mx50structs\smdata_' datestr(now,'yyyy_mm_dd')]; 
fbdata.file = ['Z:\qDots\data\mx50structs\fbdata_' datestr(now,'yyyy_mm_dd')]; 

makeNewTuneData(newDir); 
global tuneData 
tuneData.file = ['Z:\qDots\data\mx50structs\tuneData_' datestr(now,'yyyy_mm_dd')]; 

sleep