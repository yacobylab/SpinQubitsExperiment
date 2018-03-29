function awgsyncwaveforms
% Make sure the list of pulses is awgdata is consistent with the awg.
% We assume if the number of pulses is right, everything is.

global awgdata;
awgcntrl('clr');
for a=1:length(awgdata)
    npls=str2double(query(awgdata(a).awg,'WLIS:SIZE?')); % Number of waves on AWG
    if isfield(awgdata(a),'waveforms') && length(awgdata(a).waveforms) == npls % If number correct, quit
        return;
    end
    fprintf('AWG waveform list out of date.  Syncing.');
    awgdata(a).waveforms=cell(npls,1);
    awgdata(a).waveforms={};
    for i=1:npls
        r=query(awgdata(a).awg,sprintf('WLIS:NAME? %d',i-1));
        awgdata(a).waveforms{i}=r(2:end-2);
    end
    fprintf('...  Done.\n');
end
end