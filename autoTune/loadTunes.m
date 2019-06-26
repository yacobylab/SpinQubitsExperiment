function [data,out]=loadTunes(data,opts,pat)
global tuneData
if isempty(data) || ischar(data) || numel(data) == 1
    if (~exist('data','var') || isempty(data)) && ~isopt(opts,'last')
        [data,out.scan,~,out.time]=loadAna(['sm_',pat, '*']);        
    elseif exist('data','var') && ~isempty(data) && ischar(data)
        [data,out.scan,~,out.time]=loadAna(data);        
    else
        side = upper(tuneData.activeSetName(1));
        if isopt(opts,'last'), data = tuneData.runNumber; end        
        fpat = sprintf('sm_%s%s_%04i_(\\d{3}).mat', pat,side,data);
        tuneFiles = dir(tuneData.dir); tuneFiles = {tuneFiles.name};
        fnumsCell = regexp(tuneFiles,fpat,'tokens');
        fnumsCell = [fnumsCell{:}];
        fnums = str2double([fnumsCell{:}]);
        num = max(fnums);
        fileName = sprintf('sm_%s%s_%04i_%03i.mat', pat,side,data,num);
        [data,out.scan,~,out.time]=loadAna(fileName);        
        
    end
    out.anaData=1;
else
    out.anaData=0; out.time = now; 
end
end