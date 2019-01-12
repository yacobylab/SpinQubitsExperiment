function cycleThrough(fileList, nFiles,opts,varargin)
% call PlotChrg for a long list of files nFiles at a time. 
%function cycleThrough(fileList, nFiles,opts,varargin)
% nFiles default is 9 
% plotChrg called with opts and varargin. 

if ~exist('nFiles','var'), nFiles = 9; end
nplots = ceil(length(fileList)/nFiles);
for j = 1:nplots
    if j < nplots
        plotChrg(opts,fileList((j-1)*nFiles+1:(j*nFiles)),varargin{:});
    else
        plotChrg(opts,fileList((j-1)*nFiles+1:end),varargin{:});
    end
end
end