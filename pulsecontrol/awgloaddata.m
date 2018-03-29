function awgloaddata
% load latest awgdata file saved by awgsavedata.

global awgdata; global plsdata;
d = dir(sprintf('%sawgdata_*', plsdata.grpdir));
[~, mi] = max([d.datenum]);
load([plsdata.grpdir, d(mi).name]);
if exist('awgdata','var') && isfield(awgdata,'awg')
    for a=1:length(awgdata)
        data(a).awg = awgdata(a).awg;
    end
end
awgdata=data;
end
