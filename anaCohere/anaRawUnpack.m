function rawData=anaRawUnpack(scan, data)
% Unpack raw data in a scan into a more usable format.
% rawData=anaRawUnpack(scan, data)
% Output data has format data{channel}(group, pulse, rep)
% Provide scan and data as args or will ask you to load data file. 
% Raw data is not averaged over "nloop," so this file expands over that
% dimension. 

if ~exist('scan','var')
    s = load(uigetfile('sm*.mat'));
    scan = s.scan; data = s.data;
end

nchans = length(scan.loops(1).getchan);
if iscell(scan.data.conf.datachan) 
    dataChan = length(scan.data.conf.datachan);
else
    dataChan=1;
end

rawChans = nchans+(1:dataChan); % Assume first nchans not raw data, raw data comes last. 
% (As is configured in pulsecontrol). 
if ismatrix(data{rawChans(1)}) % Permute data in raw data cells. 
    for i=1:dataChan
        data{rawChans(i)} = permute(data{rawChans(i)},[1 3 2]);
    end
end
rawData=cell(dataChan,1);
npulse=scan.data.pulsegroups(1).npulse(1);
for i=1:dataChan % Permute and reshape data, put in new cell. 
    rawData{i}=permute(data{rawChans(i)},[3,1,2]);
    pointsPerGrp=size(rawData{i},1)*size(rawData{i},2); 
    rawData{i}=reshape(rawData{i},npulse,pointsPerGrp/npulse,size(rawData{i},3));
    rawData{i}=permute(rawData{i},[3 1 2]);
end
end