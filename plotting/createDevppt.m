function createDevppt(side,opts,dev)
% Makes ppt with tuning data for device. 
% function createDevppt(side,opts,n)
% sides: L, R, empty string.
% dev: which number device (these shoudl be defined in configurePlotting/2)
% opts: update, addSide
%   update: Add new files to ppt instead of starting anew.
%   addSide: Add files from new qubit to ppt.
%   files: don't update files, just use current set. 
% dev is number for which device
global pptdata;

if ~exist('dev','var') || isempty(dev),  dev = length(pptdata.dev); end
if ~exist('opts','var'), opts = ''; end
if ~isopt(opts,'files'), updateFiles; end

pptdata.start = 0;
if isopt(opts,'update')
    pptdata.filename = sprintf('dev%d%s',dev,side);
    pptControl('load')
    autoPlot(dev,pptdata.next,'',side,'new');
else
    if ~isopt(opts,'addSide')
        pptControl('start')
        pptdata.start = 1;
    else
        pptControl('load')
    end
    for j = 1:size(pptdata.dev(dev).qpcRng,1)
        autoPlot(dev,-j,1,side,'old');
    end
    for j = 1:size(pptdata.dev(dev).IndRng,1)
        autoPlot(dev,j,1,side,'old');
    end
    for j = 1:size(pptdata.dev(dev).qpcRng,1)
        autoPlot(dev,-j,'',side,'old color');
    end
    for j = 1:size(pptdata.dev(dev).IndRng,1)
        autoPlot(dev,j,'',side,'old color');
    end
end

pptdata.filename = sprintf('dev%d',dev);
pptControl('save'); pptControl('end');
%% Save the file number start adding next set of data.
numsCell=regexp([pptdata.dir.name],'(\d{4}).mat','tokens');
nums = cellfun(@str2num,[numsCell{:}]);
lastnum = max(nums); nextnum = lastnum+1;
pptdata.next(1) = nextnum;
numsCell=regexp([pptdata.qpcDir.name],'(\d{4}).mat','tokens');
nums = cellfun(@str2num,[numsCell{:}]);
lastnum = max(nums); nextnum = lastnum+1;
pptdata.next(2) = nextnum;
end