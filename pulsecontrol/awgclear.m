function awgclear(groups,opts)
% remove a set of waveforms from the AWG
% awgclear(groups,opts)
% groups: pack, all, unused
%   'pack' removes all groups, adds back groups loaded in sequences
%   'all' removes all groups
% opts: 
%   paranoid removes all waveforms, including those not known to be loaded.
%   unused: removes all waveforms not in sequence list. 
% can also be a list of groups to remove. 

global awgdata; 
if ~exist('opts','var'), opts=''; end
if strcmpi(groups, 'pack')
   grps={awgdata(1).pulsegroups.name};
   awgclear('all',opts);   
   awgadd(grps);
elseif strcmpi(groups, 'all') % Just delete all the waveforms and remove all groups
    for a=1:length(awgdata)
      fprintf(awgdata(a).awg,'WLIS:WAV:DEL ALL\n');
      awgdata(a).zeropls=[];
    end  
    logentry('Cleared all pulses.');
    awgrm('all');
else
    if strcmpi(groups,'unused')
        opts = 'unused'; 
        waveList=awgwaveforms;
        if isempty(waveList), waveList={}; end
        if ~isempty(awgdata(1).pulsegroups)
            groupList={awgdata(1).pulsegroups.name};
        else
            groupList={};
        end
        groups=setdiff(waveList,groupList); % Set of waveforms on AWG that aren't in sequences on the AWG.
    end
    if ischar(groups), groups = {groups}; end
    for a=1:length(awgdata) % Remove all the waves listed (not sequences) from the AWG.
        if isreal(groups)
            groups = sort(groups, 'descend');
            for i = groups
                waveName = query(awgdata(a).awg, sprintf('WLIS:NAME? %d', i)); %#ok<*PFCEL> returns name of waveform of index i
                if ~query(awgdata(a).awg, sprintf('WLIS:WAV:PRED? %s', waveName), '%s\n', '%i') % Check if predefined (i.e. comes with inst).
                    fprintf(awgdata(a).awg, 'WLIS:WAV:DEL %s', waveName); % Delete the waveform.
                end
            end
            awgcntrl('wait'); % Wait until finished.
            return;
        end
    end
    for k = 1:length(groups)
        for a=1:length(awgdata)
            wfms=awgwaveforms(groups{k},a,'delete');
            for i=1:length(wfms)
                fprintf(awgdata(a).awg, sprintf('WLIS:WAV:DEL "%s"', wfms{i})); % Delete waveforms on list 
            end
        end
        logentry('Cleared group %s.', groups{k});
        if ~isopt(opts,'unused')
            awgrm(groups{k});
        end
    end
end
end