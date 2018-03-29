function seqind = awgseqind(pulses,rep)
% Find the pulse line associated with a pulse group or pulse index.
% seqind = awgseqind(pulses, rep)
% Returns negative number for groups, positive for pulse index.
% Pulse can be given as a struct with fields pulses, rep. 
% Retuns nan in sequence not loaded. 
% Seems to be something wonky where can't provide rep. 
global awgdata;
if isstruct(pulses)
    rep=[pulses.rep];
    pulses=[pulses.pulses];
elseif ischar(pulses)
    pulses = {pulses};
end
seqind = nan(1, length(pulses));
for i = 1:length(pulses)
    if iscell(pulses) % could allow ints as well
        ind = strcmp(pulses{i}, {awgdata(1).pulsegroups.name});
        if isempty(ind) || ~any(ind) % no such group
            seqind(i) = nan;
        else
            seqind(i) = awgdata(1).pulsegroups(ind).seqind;
        end
    elseif pulses(i) > 0
        if exist('rep', 'var')
            ind = find(pulses(i) == awgdata(1).seqpulses);
            ind=ind(rep(i));
        else
            ind = find(pulses(i) == awgdata(1).seqpulses, 1);
        end
        if ~isempty(ind)
            seqind(i) = ind;
        end
    else
        seqind(i) = awgdata(1).pulsegroups(-pulses(i)).seqind;
    end
end
if any(isnan(seqind))
    fprintf('WARNING: Some pulses not present in sequence.\nHit Ctrl-C to abort, or any key to continue.\n');
    pause;
end
end