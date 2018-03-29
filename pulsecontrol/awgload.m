function zerolen = awgload(grp, ind)
% load pulses from group to AWG. 
% zerolen = awgload(grp,ind)

global awgdata;
awgcntrl('stop'); %changing pulses while running is slow.
awgsyncwaveforms; % make sure the waveform list is up-to-date.
doSave = false;
for a=1:length(awgdata)   
    clockInd=find([grp.pulses(1).data.clk] == awgdata(a).clk);           
    [~,offsetchan, awgdata(a).zerochan] = unique(awgdata(a).offset./awgdata(a).scale);
    offsets=awgdata(a).offset(offsetchan);    
    % create trig pulse (and corresponding 0) if waveform list empty.
    if query(awgdata(a).awg, 'WLIS:SIZE?', '%s\n', '%i') == 25 % nothing loaded (except predefined) (wlis size = size of waveform list)
        trigData=zeros(1,awgdata(a).triglen); % triglen: length of trigger in ns
        trigMarker=ones(1,awgdata(a).triglen); % channel 1&3 marker 1
        awgloadwfm(a,trigData,trigMarker,sprintf('trig_%08d',awgdata(a).triglen),1,1);
        trigMarker=zeros(1,awgdata(a).triglen); 
        for i=1:length(offsets) % Make trig pulses with correct offsets. 
            awgloadwfm(a,trigData,trigMarker,sprintf('zero_%08d_%d',awgdata(a).triglen,i),offsetchan(i),1);
        end        
        awgdata(a).zeropls = awgdata(a).triglen;
        doSave = 1;
    end
    nzpls=0;    
    for i = 1:length(grp.pulses)
        npts = size(grp.pulses(i).data(clockInd).wf, 2); % number of points (ns) in waveform. 
        if ~any(awgdata(a).zeropls == npts) % create zero if not existing yet
            zeroData=zeros(1,npts);
            for l=1:length(offsets)
                zeroName=sprintf('zero_%08d_%d', npts, l);
                awgloadwfm(a,zeroData,zeroData,zeroName,offsetchan(l),1);
            end
            awgdata(a).zeropls(end+1) = npts;
            doSave = 1;
        end
        for j = 1:size(grp.pulses(i).data(clockInd).wf, 1) %
            if isfield(awgdata(a),'virtualChans') && ~isempty(awgdata(a).virtualChans)    %see if we are mapping virtual channels onto physical %channels
                ch=find(grp.chan(j)==awgdata(a).virtualChans);
            else
                ch=find(grp.chan(j)==awgdata(a).chans);
            end
            zerolen{a}(ind(i), j) = npts;  %#ok<*AGROW>
            if isempty(ch), continue; end                        
            if any(abs(grp.pulses(i).data(clockInd).wf(j, :)) > awgdata(a).scale(ch)/2^awgdata(a).bits) || any(grp.pulses(i).data(clockInd).marker(j,:) ~= 0) % if anything nonzero           
                name = sprintf('%s_%05d_%d', grp.name, ind(i), j);                
                if any(strcmp(name,awgdata(a).waveforms))
                    fprintf(awgdata(a).awg, 'WLIS:WAV:DEL "%s"', name);
                else
                    awgdata(a).waveforms{end+1}=name;
                end
                fprintf(awgdata(a).awg, sprintf('WLIS:WAV:NEW "%s",%d,INT', name, npts)); % creates new empty waveform
                
                err = query(awgdata(a).awg, 'SYST:ERR?');
                if ~isempty(strfind(err,'E11113'))
                    fprintf(err(1:end-1));
                    error('Error loading waveform; AWG is out of memory. Try awgclear(''all''); ');
                end
                awgloadwfm(a,grp.pulses(i).data(clockInd).wf(j,:), uint16(grp.pulses(i).data(clockInd).marker(j,:)), name, ch, 0);
                zerolen{a}(ind(i), j) = -npts;
                nzpls=1;
            end
        end
    end               
    if nzpls == 0  % If no non-zero pulses were loaded, make a dummy waveform so awgclear knows this group was in memory.
        name=sprintf('%s_1_1',grp.name);
        npts=256;
        if ~any(strcmp(name,awgdata(a).waveforms))
            fprintf(awgdata(a).awg, sprintf('WLIS:WAV:NEW "%s",%d,INT', name, npts));
            awgdata(a).waveforms{end+1}=name;
        end
        awgloadwfm(a,zeros(1,npts), zeros(1,npts), name, 1, 0);
    end
end
grpInd=awggrpind(grp.name); % if the pulse group is added, update it's load time.
if ~isnan(grpInd)
    for i=1:length(awgdata)
      awgdata(i).pulsegroups(grpInd).changed = 0; 
    end
end
if doSave, awgsavedata; end
end

function awgloadwfm(a, data, marker, name, chan,define) 
% Send waveform 'data,marker' to awg a with name 'name' intended for channel chan.
% data is scaled and offset by awgdata.scale and awgdata.offset *before* sending.
% a is the awg index.
global awgdata;
if exist('define','var') && define
    fprintf(awgdata(a).awg, sprintf('WLIS:WAV:NEW "%s",%d,INT', name, length(data)));
    awgdata(a).waveforms{end+1}=name;
end
chunksize=65536; % 2^16
if(size(data,1) > size(data,2)), data=data'; end
data=(awgdata(a).offset(min(chan,end)) + data)./awgdata(a).scale(chan) + 1;
tb=find(data > 2);
tl=find(data < 0);
if ~isempty(tb) || ~isempty(tl)
    fprintf('Pulse "%s", channel %i, exceeds allowed range on awg %d: %g - %g (0-2 allowed)\n',name,chan,a,min(data),max(data));
    data(tb) = 2;
    data(tl) = 0;
end % 14 bit data offset is hard-coded in the AWG.
data = uint16(min(data*(2^(14-1) - 1), 2^(14)-1)) + uint16(marker) * 2^14;
npts = length(data);
for os=0:chunksize:npts
    if os + chunksize >= npts
        fwrite(awgdata(a).awg, [sprintf('WLIS:WAV:DATA "%s",%d,%d,#7%07d', name, os, npts-os,2 * (npts-os)),typecast(data((os+1):end), 'uint8')]); % transfer waveform data 
    else
        fwrite(awgdata(a).awg, [sprintf('WLIS:WAV:DATA "%s",%d,%d,#7%07d', name, os, chunksize,2 * chunksize),typecast(data((os+1):(os+chunksize)), 'uint8')]);
    end
    fprintf(awgdata(a).awg,'');
end
end