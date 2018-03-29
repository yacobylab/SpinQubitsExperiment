function awglist(ind,awg)
% List the waveforms present on an AWG.  If no arguments given, list all waveforms.
% function awglist([ind],[awg])
% If ind is negative, list the last (-ind) waveforms.
% awg specifies which awg in awgdata(i) to use, and defaults to 1.  Nominally the output should not depend on awg.

global awgdata;
if ~exist('awg','var') || isempty(awg)
    awg=1;
end
nwaves = query(awgdata(awg).awg, 'WLIS:SIZE?', '%s\n', '%i');
if ~exist('ind','var') || isempty(ind)    
    ind = 1:nwaves-1;
elseif ind < 0    
    ind = nwaves+(ind:-1);
end
for i = ind
    wf = query(awgdata(awg).awg, sprintf('WLIS:NAME? %d', i));
    if ~query(awgdata(awg).awg, sprintf('WLIS:WAV:PRED? %s', wf), '%s\n', '%i') % check this isn't a built in wave. 
        fprintf('%i: %s', i, wf);
    end
end
end
