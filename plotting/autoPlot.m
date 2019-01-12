function out=autoPlot(devs,cooldowns,scans,sideList,opts,filename)
% Automatically plot all data from a given category. 
% function out=autoPlot(devs,cooldowns,scans,sideList,opts,filename)
% Use with configurePlotting 
% devs should be an array of numbers referring to devices.  
% cooldowns should be the index of the cooldowns, ..unless opts is recent? 
% give a negative number for qpc cooldowns. 
% sideList : set of sides to use 1 = L, 2 = R, 3 = none 
% scans: set of which scans to plot. 
% filename: 
% opts: 
%   save: save file at the end. 
%   end: close ppt at the end. 
%   color: do all 2D scans. 
%   noppt, show: args for plotChrgB. mostly for debugging. 
%   old
%   hyst: analyze hyst on 1d data
%   new: just plot new data since last one (which was logged... -- writing
%   to file???) 
% 1: qpc, 5: hyst, 10: NT, 12: Wall 
% Devs: 1: ResDev, Ozone, Ebeam, Old, EBeam Res1,___ 3qubit, EBeamRes2
% Defaults are to run for everything. 
global pptdata;
if ~exist('opts','var'),    opts = '';  end
if ~exist('devs','var') || isempty(devs) 
    devs = 1:length(pptdata.dev);     
end
if exist('cooldowns','var') && ~isempty(cooldowns) && length(devs)>1
    error('Can only do specific cooldowns for single device')
end
dev = pptdata.dev; cats = pptdata.cats; 
scanTypeList = [cats{:}];
if ~exist('scans','var') || isempty(scans) && ~isopt(opts,'color')
    scans = 1:length(scanTypeList); 
elseif isopt(opts,'color')
    scans = length([pptdata.cats{1:2}])+1:length([pptdata.cats{:}]);
end % just plot the 2d scans. 
sides = {'L','R',''}; 
if ~exist('sideList','var') || isempty(sideList)
    sideList = 1:3; 
elseif ischar(sideList) || iscell(sideList)
    sideList = find(strcmpi(sides,sideList)); 
end
if exist('cooldowns','var')  && length(cooldowns) == 1 && length(devs) == 1 && ~isopt(opts,'new') %figure out legend name for qpc
    if cooldowns > 0
        name = dev(devs(1)).CooldownDesc{cooldowns};
    else
        name = dev(devs(1)).qpcDesc{-cooldowns};
    end
else
    name = '';
end
if isopt(opts,'noppt')  % define options to call plotChrgB with. 
    defopt = 'noppt '; 
elseif isopt(opts,'show')
    defopt = 'autoword';
else
    defopt = 'autoword invis';
end
%%  
out = {};
for m = devs
    if exist('filename','var') && ~isempty(filename)
        pptname = filename;
    else
        pptname = sprintf('dev%d',m);
    end
    pptdata.filename = pptname; 
    if ~isopt(opts,'noppt') && ~isopt(opts,'old') && ~isopt(opts,'new'),	pptControl('start'),	end    % Create new PPT
    for k = sideList
        for i = scans
            namePat = sprintf(scanTypeList{i},sides{k});            
            if i > 1 && i ~= 5 % Most names follow basic format. just qpcs and hyst don't. 
                fpat =sprintf('sm_%s_(\\d{4}).mat', namePat);
            elseif i == 1 && k == 1 % only one side for qpc
                fpat = 'sm_qpc_.*_(\d{4}).mat';         
            elseif (k == 1 || k == 2) && i == 5 %only do sides for hyst
                 fpat = sprintf('sm_%s_(\\d{4}).mat', namePat);
            else
                continue
            end
            if isopt(opts,'new') % just plot the new files (since last ppt plotting). 
                finFiles = grabFiles(pptdata.fileNames, fpat, [],cooldowns(1));
            else
                finFiles = grabFiles(pptdata.fileNames, fpat, cooldowns,dev(m).IndRng); % Takes all the files in the Ind rngs, for that cooldown, etc.            
            end
            if isfield(dev,'qpcRng') && ~isempty(dev(m).qpcRng)
                qpcRng = dev(m).qpcRng;
                if isopt(opts,'new')
                    qpcFiles = grabFiles(pptdata.qpcfileNames, fpat, [],cooldowns(2));
                else
                    qpcFiles = grabFiles(pptdata.qpcfileNames, fpat, -cooldowns,qpcRng);
                end
                fileList=cellfun(@(x) fullfile(pptdata.qpcFolder,x),qpcFiles,'UniformOutput',false);
                finFiles = [finFiles fileList];
            end
            if ~isempty(finFiles)
                if any(strcmpi(scanTypeList{i},cats{1})) && k == 1 % qpc
                    if ~isopt(opts,'hyst')
                        [~,out]=qpcPlot([defopt 'filt'],finFiles,[dev(m).name ' ' name]);
                        out.name = [m,cooldowns]; 
                    else
                        [~,out]=qpcPlot([defopt 'filt hyst'],finFiles,[dev(m).name ' ' name]);                    
                        out.name = [m,cooldowns]; 
                    end
                elseif any(strcmpi(scanTypeList{i},cats{2})) % 
                        [~,out]=qpcPlot2([defopt, 'mean hyst'],finFiles);               
                         out.name = [m,cooldowns];                             
                elseif any(strcmpi(scanTypeList{i},cats{3})) % 
                    plotChrg([defopt, 'water stop'], finFiles);
                elseif any(strcmpi(scanTypeList{i},cats{4})) % 
                    cycleThrough(finFiles, 9,[defopt, 'cbar']);
                elseif any(strcmpi(scanTypeList{i},cats{5})) % 
                    cycleThrough(finFiles, 9,[defopt, 'cbar glitch difffirst']);
                elseif any(strcmpi(scanTypeList{i},cats{6})) % 
                    cycleThrough(finFiles, 9,[defopt, 'chrg sub stop']);
                else
                    fprintf('Category %s not described \n ', scanTypeList{i});
                end       
            end
        end
    end
    if ~isopt(opts,'noppt') && isopt(opts,'save'), pptControl('save');     end
    if isopt(opts,'end'),	pptControl('end');	end
end
end