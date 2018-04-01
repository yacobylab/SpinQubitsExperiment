function createDevppt(side,opts,n)
%function createDevppt(side,opts,n)
% sides: L, R, empty string. 
% opts: update, addSide 
% n is number for which device 
global pptdata; 

if ~exist('n','var') || isempty(n) 
    n = length(pptdata.dev); 
end
if ~exist('opts','var') 
    opts = ''; 
end
if ~isopt(opts,'files'), updateFiles;   end
pptdata.start = 0; 
if isopt(opts,'update')        
    pptdata.filename = sprintf('dev%d',n); 
    pptControl('load')
    autoPlotFunc(n,pptdata.next,'',side,'new');
else
    if ~isopt(opts,'addSide')
        pptControl('start')
        pptdata.start = 1;
    else
        pptControl('load')
    end
    for j = 1:size(pptdata.dev(n).qpcRng,1)
        autoPlotFunc(n,-j,1,side,'old');
    end
    for j = 1:size(pptdata.dev(n).IndRng,1)
        autoPlotFunc(n,j,1,side,'old');
    end        
    for j = 1:size(pptdata.dev(n).qpcRng,1)
        autoPlotFunc(n,-j,'',side,'old color');
    end
    for j = 1:size(pptdata.dev(n).IndRng,1)
        autoPlotFunc(n,j,'',side,'old color');
    end    
end

pptdata.filename = sprintf('dev%d',n); 
pptControl('save'); pptControl('end');
%% Save where to start for adding next set of data. 
numsCell=regexp([pptdata.dir.name],'(\d{4}).mat','tokens');
nums = cellfun(@str2num,[numsCell{:}]);
lastnum = max(nums); nextnum = lastnum+1;
pptdata.next(1) = nextnum; 
numsCell=regexp([pptdata.qpcDir.name],'(\d{4}).mat','tokens');
nums = cellfun(@str2num,[numsCell{:}]);
lastnum = max(nums); nextnum = lastnum+1;
pptdata.next(2) = nextnum; 
end