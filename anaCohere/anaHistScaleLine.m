function [data, scalefuncs, meanVals,params,volts,histData]=anaHistScaleLine(scan, data,t1s)
% Rescale raw voltage data to range 0 -> 1 using histograms from each line
% (for unstable data). 
%[data,scalefuncs, meanVals,fp,v,n]=anaHistScale(scan, data,t1s,grps)
% data is the rescaled data (from raw voltages -> 0 to 1 range.)
% scalefuncs: funcs for rescaling data 
% Rescale histogrammed data.  t1s is a vector of t1 time estimtes.
% ASSUMES NO CROSSTALK
% currently only works w/ only one channel of data.

if length(size(data{end})) == 2, data{end}=permute(data{end},[1 3 2]); end % only 1 group
nDataSets=floor(length(data)/2); 
for i=1:length(t1s)
    [fitfn, initfn] = getfn(t1s(i)); %fit function for histograms, in terms of peak spacing, noise, t1.
    procInd = [];
    for j = 1:length(scan.loops(1).procfn) %find the procfn that does histogramming
        if length(scan.loops(1).procfn(j).fn)==1 && strcmpi(func2str(scan.loops(1).procfn(j).fn.fn),'histc')
            procInd  = [procInd, j];
        end
    end
    volts=scan.loops(1).procfn(procInd(1)).fn.args{1}; %scan.loops(1).procfn(3).fn=histc, %scan.loops(1).procfn(3).args=set of histogram vals, from fbdata and fConfSeq2, scan.loops(1).procfn(3).dim=500
    volts=(volts(1:end-1)+volts(2:end))/2;    % HistC gives edges, not centers.   
    data{nDataSets+i+1}(isnan(data{end})) = 0; %any nans in histogrammed set to 0, nds+i+1 is histogram associated w/ data set i.
    if all(data{nDataSets+2}==0), error('Histogram data was all 0 or NaN. anaHistScale wont work'); end
    dataCurr = data{nDataSets+i+1}; 
    for j = 1:size(dataCurr,1)
        for k=1:size(dataCurr,2)
            histData=squeeze(dataCurr(j,k,:)); %averages over all reps and groups to find the number of elements at each voltage.
            histData(end)=[];
            histData=histData/mean(histData); % Make fitwrap happy.
            if all(isnan(histData)), continue, end
            params=fitwrap('fine',volts,histData',initfn,fitfn,[1 1 1 1 0 0 1]);
            meanVals = ((1-exp(-abs(params(5:6))))./abs(params(5:6))-.5).*[-1 1]./params(2) + params(1); % Mean voltages for singlet, triplet
            scalefuncs{i,j,k}=makescalefunc(1/diff(meanVals),-meanVals(1)/diff(meanVals)); %takes inverse peak spacing and offset w.r.t inverse peak spacing to rescale from 0 to 1.
            data{i}(j,k,:) = scalefuncs{i,j,k}(data{i}(j,k,:));
        end
    end
    % uncomment below to fit t1.
    %fp=fitwrap('plfit samefig resid',v,n',params,fitfn,[1 1 1 1 1 1 1]);        
end
end

function f=makescalefunc(scale,off)
f=@(x) x*scale+off;
end

function [fitfn, initfn] = getfn(t1)
if isnan(t1), t1=1e-5; end
distfn = @(a, x) exp(-a(1)) * exp(-(x-1).^2./(2* a(2)^2))./(sqrt(2 * pi) * a(2)) + ...
    a(1)/2 * exp(a(1)/2 * (a(1) * a(2)^2 - 2 * x)) .* ...
    (erf((1 + a(1) * a(2)^2 - x)./(sqrt(2) * a(2))) + erf((-a(1) * a(2)^2 + x)./(sqrt(2) * a(2))));
% parameters: [t_meas/T1, rms amp noise/peak spacing]

fitfn = @(a, x) a(3) * distfn(abs(a([5 7])), .5-(x-a(1)).*a(2)) + a(4) * distfn(abs(a([6 7])), .5+(x-a(1)).*a(2));
%parameters: [ center between peaks, 1/spacing, coeff left peak, coeff right peak, t_m/T1,left, t_m/T1,right, rms amp noise/peak spacing]
% sum(coefficients) = 1 corresponds to a PDF for unity peak spacing.
% If fitting raw histograms, # samples = sum(fp(:, 3:4), 2) ./(fp(:, 2) * diff(d.x(1:2)));

initfn.fn = @(x, y)[sum(x.*y)/sum(y),  1/sqrt(sum(y.*((x-sum(x.*y)/sum(y)).^2)/sum(y))), max(y), max(y), 1e-5, t1, .4];
initfn.args = {};
end 