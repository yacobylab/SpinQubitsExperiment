function awggroups(ind)
% Print all awggroups if ind is absent, or the names of groups specified by ind otherwise.
% function awggroups(ind)

global awgdata;
if ~exist('ind','var') || isempty(ind)
    ind = 1:length(awgdata(1).pulsegroups);
end
for i = ind
    zl=awgdata(1).pulsegroups(i).zerolen;
    fprintf('%2i:  %-15s  (%3i pulses, %5.2f us, %d lines)\n', i, awgdata(1).pulsegroups(i).name, awgdata(1).pulsegroups(i).npulse(1), abs(zl(1)*1e6/awgdata(1).clk),awgdata(1).pulsegroups(i).nline);
end
end