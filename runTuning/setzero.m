function scan=setzero(scan)
% Set the gates in the scan to zero at the end. Use as configfn.
%function scan=setzero(scan)
% to set up, use scan.configfn(end+1).fn = @setzero; 
%                scan.configfn(end).args = {};
scan.cleanupfn(1).fn = @smaconfigwrap;
scan.cleanupfn(1).args{1} = @smset;
scan.cleanupfn(1).args{2} = [scan.loops.setchan];
scan.cleanupfn(1).args{3} = 0;
end