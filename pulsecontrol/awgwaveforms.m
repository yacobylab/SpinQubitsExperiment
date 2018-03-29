function waveforms = awgwaveforms(group,awg,opts)
% Give a list of waveforms that are known to be loaded, return list of names.
% function waveforms = awgwaveforms(group,awg,opts)
% no arguments: return a list of group names
% if group is specified, give the names of the waveforms in that group
% if awg is specified, work on awg awg
% if 'opts' is 'delete' and group is specifed, erase the waveforms from the
% table (does not do anything to AWG). 

global awgdata;
awgsyncwaveforms();
if ~exist('awg','var') || isempty(awg), awg=1; end
if ~exist('group','var') || isempty(group)
    waveforms = regexp(awgdata(awg).waveforms,'(.*)_\d+_\d+','tokens');
    waveforms = [waveforms{:}];
    waveforms = [waveforms{:}];
    waveforms=unique(sort(waveforms));
    i=strcmp('zero',waveforms);
    waveforms(i)=[];
    i=strcmp('trig',waveforms);
    waveforms(i)=[];
else
    waveforms = regexp(awgdata(awg).waveforms,sprintf('(%s)_\\d+_\\d+',group));
    ind = find(cellfun(@(x) ~isempty(x) && x == 1 , waveforms));    
    waveforms=awgdata(awg).waveforms(ind);
    if exist('opts','var') && strcmp(opts,'delete'), awgdata(awg).waveforms(ind)=[]; end
end
end