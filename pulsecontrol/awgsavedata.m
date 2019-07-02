function awgsavedata
% save awgdata in plsdata.grpdir, with name generated from date and time.

global awgdata; global plsdata;

data = rmfield(awgdata, 'awg'); %#ok<*NASGU>
data = rmfield(data, 'pulsegroups'); %#ok<*NASGU>
time = clock;
save(sprintf('%sawgdata_%02d%02d%02d_%02d%02d', plsdata.grpdir, mod(time(1), 100), time(2:5)), 'data','-v6');
end