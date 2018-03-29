function awgrm(grp, ctrl)
% Remove a set of groups currently loaded on awg. 
% awgrm(grp, ctrl)
% grp: 'all' or group name
% ctrl: 'after' remove all following groups
% Otherwise, specified group is removed by removing it and following ones and then reloading the latter.

global awgdata;
if strcmp(grp, 'all')
    awgcntrl('stop');
    for a=1:length(awgdata)
        fprintf(awgdata(a).awg, 'SEQ:LENG 0'); % Clears all sequences on AWG.
        awgdata(a).pulsegroups = [];
        awgdata(a).seqpulses = [];
    end
    awgsavedata;
else
    grp = awggrpind(grp);
    grp = grp(1); 
    if isnan(grp)
        warning('First group is not currently loaded')
        return;
    end 
    awgcntrl('stop');
    if exist('ctrl','var') && isopt(ctrl, 'after')
        for a=1:length(awgdata)
            fprintf(awgdata(a).awg, 'SEQ:LENG %d', awgdata(a).pulsegroups(grp).seqind-1 + sum(awgdata(a).pulsegroups(grp).nline)); % Set nseq to the index of final sequence plus # lines
            awgdata(a).seqpulses(awgdata(a).pulsegroups(grp).seqind + sum(awgdata(a).pulsegroups(grp).npulse):end) = [];
            awgdata(a).pulsegroups(grp+1:end) = [];
        end
    else
        for a=1:length(awgdata)
            fprintf(awgdata(a).awg, 'SEQ:LENG %d', awgdata(a).pulsegroups(grp).seqind-1);
            awgdata(a).seqpulses(awgdata(a).pulsegroups(grp).seqind:end) = [];
            groups = {awgdata(a).pulsegroups(grp+1:end).name};
            awgdata(a).pulsegroups(grp:end) = [];
            awgadd(groups);
        end
    end
end
end