function newDataDir(opts)
% It's a good policy to change folders once there are a few thousand items
% in it, as startup will slow down at that point. 
% You will want to create a new data folder, and copy the data structs to
% the new experiment. There are (approximately) three possibilities of why
% you are starting a new directory, and each should be dealt with
% differently. 
% This assumes you already have the structs loaded and that they are globals. We'll add a version
% that copies files as well later. 
% We will also add plsdata handling. 

% 1) Same experimental setup, same sample (e.g. if folder has too many
% files) 
%      -Copy all the structs with new names. 
% 2) Same experimental setup, new sample 
%       -Copy smdata, otherwise create new default structs. 
% 3) New experimental setup, new sample. 
%       -Create new smdata. This is not yet set up.
% 1 is default. For 2, give opt "newSample". 

% The location of each struct is stored in two places: smdata.files and in
% the struct itself. The smdata location is useful for loading all of the
% structs, and storing in the struct itself good for ease when saving
% manually. 
if ~exist('opts','var'), opts = ''; end 

baseDir = 'Z:/qDots/data/';
structFolder = 'mx50structs/'; 
dateInfo = datestr(now,'yyyy_mm_dd'); 
newDir =fullfile(baseDir,['data_' dateInfo]);
mkdir(newDir); 
cd(newDir)
mkdir(['qpc_' dateInfo]); 
global tuneData; global scandata; global fbdata; global smdata; 
if isopt(opts,'newSample') 
    makeNewTuneData(newDir);               
end

tuneData.file = [baseDir,structFolder,'tuneData_',dateInfo]; 
scandata.file = [baseDir,structFolder,'scandata_',dateInfo]; 
fbdata.file = [baseDir,structFolder,'fbdata_',dateInfo]; 
smdata.files.smdata = [baseDir,structFolder,'smdata_',dateInfo]; 
smdata.files.scandata = scandata.file; 
smdata.files.fbdata = fbdata.file; 
smdata.files.tunedata = tuneData.file;
smdata.files.log = ['Z:/qDots/notes/log_' dateInfo '.txt']; 
smdata.files.dir = newDir;
smdata.files.ppt = ['Z:/qDots/PPT/' dateInfo]; % This is where pptplot puts ppt data.  
mkdir(smdata.files.ppt); 
sleep % This will save all of the new racks. 
end