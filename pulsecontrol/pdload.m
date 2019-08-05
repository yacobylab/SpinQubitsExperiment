function dict = pdload(name)
% Load most recent entry in a pulse dictionary.  
%function pd = pdload(name)

global plsdata;
if isstruct(name) % if you've passed it a dict, return. 
    dict=name;
    return;
end
load([plsdata.grpdir, 'pd_', name,'_last']);
if ~exist('dict','var'), dict = pd; end
end
