function rawData=anaRawUnpack(scan, data)
% Unpack raw data in a scan into a more usable format.
%data=anaRawUnpack(scan, data)
% Output data is in data{channel}(group, pulse, rep)
if nargin == 0
    s = load(uigetfile('sm*.mat'));
    scan = s.scan; data = s.data;
end

nchans = length(scan.loops(1).getchan);
if iscell(scan.data.conf.datachan)
    dataChan = length(scan.data.conf.datachan);
else
    dataChan=1;
end
rawChans = nchans+(1:dataChan);
if ismatrix(data{rawChans(1)})
    for i=1:dataChan
        data{rawChans(i)} = permute(data{rawChans(i)},[1 3 2]);
    end
end
rawData=cell(dataChan,1);
npulse=scan.data.pulsegroups(1).npulse(1);
for i=1:dataChan
    rawData{i}=permute(data{rawChans(i)},[3,1,2]);
    ppg=size(rawData{i},1)*size(rawData{i},2); % points-per-group
    rawData{i}=reshape(rawData{i},npulse,ppg/npulse,size(rawData{i},3));
    rawData{i}=permute(rawData{i},[3 1 2]);
end
end