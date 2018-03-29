function zerolen = awgzero(grp, ind, zerolen)
% determine if pulse is zero (helper function)
% zerolen = awgzero(grp, ind, zerolen,awg)
% If anything out of scale, zerolen negative.  
% Ind gives group number. 

global awgdata;
for awg=1:length(awgdata)
    for i = 1:length(grp.pulses)
        clockInd = find([grp.pulses(i).data.clk] == awgdata(awg).clk);
        npts = size(grp.pulses(i).data(clockInd).wf, 2);
        for j=1:size(grp.pulses(i).data(clockInd).wf,1) % Loops over channels. FIXME; channel mappings not honored here.
            if any(abs(grp.pulses(i).data(clockInd).wf(j,:) > awgdata(awg).scale(min(j,end))/(2^awgdata(awg).bits)))
                zerolen{awg}(ind(i),j) = -npts;
            else
                zerolen{awg}(ind(i),j) = npts;
            end
        end
    end
end
