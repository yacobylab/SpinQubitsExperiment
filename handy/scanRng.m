function vals = scanRng(scan,n)
% return linspace of scanRng. 
% vals = scanRng(scan)
if ~exist('n','var') 
    n = 1; 
end
if isfield(scan.loops,'rng') && ~isempty(scan.loops(n).rng)
    vals = linspace(scan.loops(n).rng(1),scan.loops(n).rng(2),scan.loops(n).npoints);
elseif isfield(scan.loops,'setchanranges')
    vals = linspace(scan.loops(n).setchanranges{1}(1),scan.loops(n).setchanranges{1}(2),scan.loops(n).npoints);
else
    vals = linspace(1,scan.loops(n).npoints,scan.loops(n).npoints);
end

end