function [data, scalefuncs, meanVals,params,histVolt,histData,fidelity]=anaHistScale(scan, data, t1s, grps,opts)
% Rescale raw voltage data to <s>=0, <t>=1.
% [data, scalefuncs, meanVals,params,histVolt,histData]=anaHistScale(scan, data,t1s,grps)
% t1 is ratio of measurement time to t1, tm/t1
% scalefuncs: funcs for rescaling data
% data: rescaled data
% meanVals: mean value for singlet and triplet states. 
% Rescale histogrammed data. 
% ASSUMES NO CROSSTALK

if ~exist('opts','var'), opts = ''; end
if length(size(data{end})) == 2 % only 1 group
    data{end}=permute(data{end},[1 3 2]);
end
if ~exist('grps','var') || isempty(grps) %second dim of histogram is # groups.
    grps=1:size(data{end},2);
end
nDataSets=floor(length(data)/2);
for i=1:length(t1s)
    procInd = [];
    for j = 1:length(scan.loops(1).procfn) %find the procfn that does histogramming
        if length(scan.loops(1).procfn(j).fn)==1 && strcmpi(func2str(scan.loops(1).procfn(j).fn.fn),'histc')
            procInd  = [procInd, j];
        end
    end
    % What these fields are: (configured in fConfSeq) 
    % scan.loops(1).procfn(3).fn=histc, 
    % scan.loops(1).procfn(3).args=set of histogram vals, from fbdata and fConfSeq
    % scan.loops(1).procfn(3).dim=500 (number of values in histogram)
    histVolt=scan.loops(1).procfn(procInd(1)).fn.args{1};
    histVolt=(histVolt(1:end-1)+histVolt(2:end))/2;    % HistC gives edges, not centers.
    % Any nans in histogrammed set to 0, ndatasets+i+1 is histogram associated w/ data set i.
    data{nDataSets+i+1}(isnan(data{end})) = 0; 
    if all(data{nDataSets+1+1}==0) % It was all populated with nans, means histogramming didin't work
        error('Histogram data was all 0 or NaN. anaHistScale wont work');
    end
    distfn = @(a, x) exp(-a(1)) * exp(-(x-1).^2./(2* a(2)^2))/sqrt(2*pi)./a(2) + a(1)/2 * exp(a(1)/2 * (a(1) * a(2)^2 - 2 * x)) ...
        .* (erf((1 + a(1) * a(2)^2 - x)/sqrt(2) ./ a(2)) + erf((-a(1) * a(2)^2 + x)/sqrt(2) ./ a(2)));
    fitfn = @(p, x) p(3) * distfn(abs(p([5, 7])), .5-(x-p(1)).*p(2)) + p(4) * distfn(abs(p([6, 7])), .5+(x-p(1)).*p(2));
    % Averages over all reps and groups to find the number of elements at each voltage.
    histData=squeeze(sum(sum(data{nDataSets+i+1}(:,grps,:),2),1)); 
    histData(end)=[];
    histData=histData/mean(histData); % Make fitwrap happy.
        
    aveV = sum(histData'.*histVolt)/sum(histData);
    sdV = sqrt(sum(histData'.*histVolt.^2)/sum(histData)-aveV^2);
    %1: mean V, 2: 1/peak spacing, 3: left peak mag 4: right peak mag 5: t/T1S 6: t/T1T, 7: noise/peak spacing
    beta0=[aveV, 1/2/sdV, 0.6*max(histData), .4*max(histData), 1e-4, t1s(i), 0.25];
    if ~isopt(opts,'noplot')
        figure(400+i); clf; hold on;
        params=fitwrap('plfit plinit samefig fine',histVolt,histData',beta0,fitfn,[1 1 1 1 0 0 1]);
        ax=axis; axis([params(1)-3/params(2), params(1)+3/params(2) ax(3) ax(4)]); % scale the x axis nicely
    else
        params=fitwrap('noplot fine',histVolt,histData',beta0,fitfn,[1 1 1 1 0 0 1]);
    end
    
    
    meanVals = ((1-exp(-abs(params(5:6))))./abs(params(5:6))-.5).*[-1 1]./params(2) + params(1); % Mean voltages for singlet, triplet
    % Makes inverse peak spacing and offset w.r.t inverse peak spacing to rescale from 0 to 1.
    scalefuncs{i}=makescalefunc(1/diff(meanVals),-meanVals(1)/diff(meanVals)); 
    data{i} = scalefuncs{i}(data{i});
    
    % Calculate fidelity
    sampNum = length(histData);
    Sfit=params; Sfit(3)=1; Sfit(4)=0; % Only keep singlet peak.
    fitSing=(params(3)+params(4))*fitfn(Sfit,histVolt)/sampNum; % Fitted hist of just sing peak
    Sfid=cumsum(fitSing); Sfid=[0,Sfid]; %Sum cumulative prob of capturing all singlets as a function of voltage
    
    Tfit=params; Tfit(3)=0; Tfit(4)=1; % Only keep triplet peak
    fitTrip=(params(3)+params(4))*fitfn(Tfit,histVolt)/sampNum; % Fitted hist of just trip peak.
    TfidRev=cumsum(fitTrip(end:-1:1)); TfidRev=[0 TfidRev];
    Tfid=TfidRev(end:-1:1); %Sum cumulative prob of capturing all singlets as a function of voltage
        
    fidelity = max(Sfid+Tfid)/2; % Fid = (correctly identified/total);
end
end

function f=makescalefunc(scale,off)
f=@(x) x*scale+off;
end