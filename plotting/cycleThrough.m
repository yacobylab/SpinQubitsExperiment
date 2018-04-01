function cycleThrough(fileList, nFiles,opts,varargin)
% call PlotChrgB for a long list of files only 9 at a time. 
%function cycleThrough(fileList, nFiles,opts,varargin)
nplots = ceil(length(fileList)/nFiles);
for j = 1:nplots
    if j < nplots
        plotChrgB(opts,fileList((j-1)*nFiles+1:(j*nFiles)),varargin{:});
    else
        plotChrgB(opts,fileList((j-1)*nFiles+1:end),varargin{:});
    end
end
end