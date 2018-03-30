function [out, fileDir] = getFiles(ffilter)
% Returns list of file names as cell. Can use filter. 
%function [out, fileDir] = getFiles(ffilter)

if ~exist('ffilter','var'), ffilter = ''; end
[out, fileDir] = uigetfile(ffilter,'MultiSelect','On');
if isempty(out), return; end
if ~iscell(out), out = {out}; end
end