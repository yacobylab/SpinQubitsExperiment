function awgadd(groups)
% Add groups to end of sequence on AWG, set up looping and jumps. Store group info in awgdata.pulsegroups
% awgadd(groups)
% For sequence combined groups (only), pulseind is a 2D array.  Each row corresponds to a groups, each column to a varpar in this group.
% ie, pulseind = [ 1 1 1 1 ; 1 2 3 4] will use pls 1 from group 1, pulse 1-4 on group 2.

global awgdata;
awgcntrl('clr stop');
if ~iscell(groups), groups = {groups}; end
doSave = false; % keeps track of whether awgdata changed.
for k = 1:length(groups) % we will always load, even if not changed. check later.
    grpdef=plsmakegrp(groups{k},'upload'); % Make and save, upload group's waveforms.
    awgcntrl('wait'); % Wait until upload complete.
    firstCtrl = strtok(grpdef.ctrl);
    if strcmp(firstCtrl,'grp') && isopt(grpdef.ctrl, 'seq') % combine groups at sequence level.
        chan = {grpdef.chan};
        seqmerge = true;
    else
        if ~isfield(grpdef, 'pulseind')
            if iscell(grpdef.zerolen)
                grpdef.pulseind = 1:size(grpdef.zerolen{1}, 1);
            else
                grpdef.pulseind = 1:size(grpdef.zerolen, 1);
            end
        end
        seqmerge = false;
    end
    if ~isfield(grpdef, 'nrep'), grpdef.nrep = 1; end
    for a=1:length(awgdata)
        npls = size(grpdef.pulseind, 2);
        nchans = length(awgdata(a).chans); % alternatively use awgdata or data size
        useTrig = grpdef.nrep(1) ~= Inf && ~isopt(grpdef.ctrl, 'notrig');
        if isempty(awgdata(a).pulsegroups) % Figure out sequence index.
            startLine = 1; 
            grpInd = [];
        else
            grpInd = find(strcmp(grpdef.name, {awgdata(a).pulsegroups.name}));
            if ~isempty(grpInd)
                startLine = awgdata(a).pulsegroups(grpInd(1)).seqind;
                if npls + useTrig ~= sum(awgdata(a).pulsegroups(grpInd(1)).npulse)
                    warning('Number of pulses changed in group %s. Use awgrm first!', grpdef.name);
                end
            else
                startLine = awgdata(a).pulsegroups(end).seqind + sum(awgdata(a).pulsegroups(end).nline);
            end
        end        
        zerolen = grpdef.zerolen; %Length of full sequence, for setting up zero pulses on other channels. 
        if ~iscell(zerolen), zerolen={zerolen}; end
        totPulses = [npls useTrig];
        if isempty(grpInd) % group not loaded yet, add on to awgroups.
            grpInd = length(awgdata(a).pulsegroups)+1;
            awgdata(a).pulsegroups(grpInd).seqind = startLine;            
            if strfind(grpdef.ctrl,'pack') % we end up with 1 line for all pulses, set npls to 1. 
                zlmult=npls;
                npls=1;
                awgdata(a).pulsegroups(grpInd).nline = 1+useTrig;                
            else
                zlmult=1;
                awgdata(a).pulsegroups(grpInd).nline = npls+useTrig;
            end
            fprintf(awgdata(a).awg, sprintf('SEQ:LENG %d', startLine + awgdata(a).pulsegroups(grpInd).nline-1)); % Extend sequence length on AWG
            doSave = 1;
        else
            if strfind(grpdef.ctrl,'pack') % if pack, all pulses in one sequence, and zerolen is npls times length of single pulse. 
                zlmult = npls; % 
                npls=1;
            else
                zlmult=1;
            end            
            if any(awgdata(a).pulsegroups(grpInd).nrep ~= grpdef.nrep) || (isfield(grpdef,'readout') && isfield(grpdef.readout,'readout')...
                    && any(any(any(awgdata(a).pulsegroups(grpInd).readout.readout ~= grpdef.readout.readout)))) || any(any(awgdata(a).pulsegroups(grpInd).zerolen ~= grpdef.zerolen{a})) %#ok<*COLND> % nrep or similar changed
                doSave = 1;
            end
        end
        if isfield(grpdef,'pulses') && isfield(grpdef.pulses,'data')
            for j = 1:length(grpdef.pulses)
                if isfield(grpdef.pulses(j),'data') && isfield(grpdef.pulses(j).data,'marker')
                    grpdef.pulses(j).data = rmfield(grpdef.pulses(j).data,'marker');
                end
                if isfield(grpdef.pulses(j),'data') && isfield(grpdef.pulses(j).data,'wf')
                    grpdef.pulses(j).data = rmfield(grpdef.pulses(j).data,'wf');
                end
            end
            
        end
        flds = fieldnames(grpdef);        
        for j = 1:length(flds)
            awgdata(a).pulsegroups(grpInd).(flds{j}) = grpdef.(flds{j});
        end
        
        awgdata(a).pulsegroups(grpInd).npulse = totPulses; 
        awgdata(a).pulsegroups(grpInd).changed = false;        
        awgdata(a).pulsegroups(grpInd).zerolen = grpdef.zerolen{a};
        if ~isfield(grpdef, 'jump')
            if isopt(grpdef.ctrl, 'loop')
                grpdef.jump = [npls; 1];
            else
                grpdef.jump = [];
            end
        end           
        if useTrig % Set triggers on the 1st and 3rd channels.
            for j = 1:nchans
                if mod(j,2) == 1
                    fprintf(awgdata(a).awg, sprintf('SEQ:ELEM%d:WAV%d "trig_%08d"', startLine, j, awgdata(a).triglen)); % Sets waveform for sequence startline
                else
                    fprintf(awgdata(a).awg, sprintf('SEQ:ELEM%d:WAV%d "zero_%08d_%d"', startLine, j, awgdata(a).triglen, awgdata(a).zerochan(j)));
                end
            end
            if isfield(awgdata(a),'slave') && ~isempty(awgdata(a).slave) && awgdata(a).slave % set wait trigger state on.
                fprintf(awgdata(a).awg, sprintf('SEQ:ELEM%d:TWAIT 1\n', startLine));
            end
        end
        for i = 1:npls
            ind = i-1 + startLine + useTrig;
            if ~seqmerge % pulses combined here.
                for j = 1:nchans
                    if isfield(awgdata(a),'virtualChans') && ~isempty(awgdata(a).virtualChans)
                        ch=find(awgdata(a).virtualChans(j)==grpdef.chan);
                    else
                        ch=find(awgdata(a).chans(j)==grpdef.chan);
                    end
                    if ~isempty(ch) &&  grpdef.zerolen{a}(grpdef.pulseind(i), ch) < 0  % % Set sequence to correct wave if channel in group and not zero
                        fprintf(awgdata(a).awg, sprintf('SEQ:ELEM%d:WAV%d "%s_%05d_%d"', ind, j, grpdef.name, grpdef.pulseind(i), ch));
                    else
                        fprintf(awgdata(a).awg, sprintf('SEQ:ELEM%d:WAV%d "zero_%08d_%d"', ind, j, zlmult*abs(zerolen{a}(grpdef.pulseind(i), 1)),awgdata(a).zerochan(j)));
                    end
                end
            else
                for m = 1:length(grpdef.pulses.groups)
                    for j = 1:length(chan{m}) % channels of component groups
                        ch = find(awgdata(a).chans == chan{m}(j));
                        if ~isempty(ch) %zero replacement not implemented
                            fprintf(awgdata(a).awg, sprintf('SEQ:ELEM%d:WAV%d "%s_%05d_%d"', ind, ch,grpdef.pulses.groups{m}, grpdef.pulseind(m, i), ch));
                        else
                            error('This won''t work');
                        end
                    end
                end
            end
            if isopt(grpdef.ctrl,'single') %for a single, we want the group to run just once, then go to pulseline 1. this will set it to run once
                if isopt(grpdef.ctrl,'loop')
                    error('trying to use ''loop'' and ''single'' for same group')
                end
                fprintf(awgdata(a).awg, 'SEQ:ELEM%d:LOOP:INF 0', ind); % turn off infinite looping
                fprintf(awgdata(a).awg, sprintf('SEQ:ELEM%d:LOOP:COUN 1', ind)); % set number of reps to 1.
            elseif grpdef.nrep(min(i, end)) == Inf  || grpdef.nrep(min(i, end)) == 0 || (i == npls && ~isopt(grpdef.ctrl, 'loop') && (isempty(grpdef.jump) || all(grpdef.jump(1, :) ~= i))) %
                fprintf(awgdata(a).awg, sprintf('SEQ:ELEM%d:LOOP:INF 1', ind)); % Turn on infinite looping
            else
                fprintf(awgdata(a).awg, 'SEQ:ELEM%d:LOOP:INF 0', ind); % default, turn off infinite looping
                fprintf(awgdata(a).awg, sprintf('SEQ:ELEM%d:LOOP:COUN %d', ind, grpdef.nrep(min(i, end)))); % Set number of reps to nreps.
            end
            fprintf(awgdata(a).awg, sprintf('SEQ:ELEM%d:GOTO:STAT 0', ind)); % Turn off goto part of sequencer.
            if grpdef.nrep(min(i, end)) == Inf && isreal(grpdef.pulses) && (length(awgdata(a).seqpulses) < ind || awgdata(a).seqpulses(ind) ~= grpdef.pulses(grpdef.pulseind(i)))
                doSave = 1;
                awgdata(a).seqpulses(ind) = grpdef.pulses(grpdef.pulseind(i));
            end
        end
        if isopt(grpdef.ctrl,'single') % Make single go back to all off pulse group when finished.
            offGrps=find(~cellfun('isempty',regexp({awgdata(1).pulsegroups.name},'^all_off')));
            if ~isempty(offGrps)
                offLine=awgseqind(-offGrps(1))+1; %+1 because we want to skip the trigger line.
            else
                offLine=1;
            end
            fprintf(awgdata(a).awg, sprintf('SEQ:ELEM%d:GOTO:IND %d', ind,offLine));
            fprintf(awgdata(a).awg, sprintf('SEQ:ELEM%d:GOTO:STAT 1', ind));
        end
        for j = 1:size(grpdef.jump, 2) % First row of jump is the pulse number, Second row is where to go to. 
            fprintf(awgdata(a).awg, sprintf('SEQ:ELEM%d:GOTO:IND %d', startLine+useTrig-1 + grpdef.jump(:, j))); %Jump some number of pulses at end.
            fprintf(awgdata(a).awg, sprintf('SEQ:ELEM%d:GOTO:STAT 1', startLine+useTrig-1 + grpdef.jump(1, j))); % Turn on GOTO state.
        end
    end
    awgcntrl('wait');
    nerr=0;
    for a=1:length(awgdata(a))
        err=query(awgdata(a).awg, 'SYST:ERR?');
        if ~contains(err, 'No error'), nerr=nerr+1; end
    end
    if nerr == 0, logentry('Added group %s on index %i.', grpdef.name, grpInd); end
end
if doSave, awgsavedata; end
end