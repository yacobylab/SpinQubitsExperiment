function [tlloc,tlfit] = tlscan(gates,scanvals,side,ind,opts)
% Takes a tl scan w/ the offset on gates given in scanvals, then fits to find tl pt.
% function [tlloc,tlfit] = tlscan(gates,scanvals,side,ind,opts)
% Returns location in mV
global tuneData;
if ~exist('opts','var'), opts = ''; end 
if isopt(opts,'amp')
    for i = 1:length(tuneData.xyChan) 
        gateInd = strcmpi(tuneData.xyChan{i},gates);
        if any(gateInd) 
            scanvals(gateInd) = scanvals(gateInd)+tuneData.measPt(i); 
        else 
            gates{end+1} = tuneData.xyChan{i}; 
            scanvals(end+1) = tuneData.measPt(i); 
        end
    end
end
scan = fConfSeq(tuneData.tl.plsgrp,{'nloop',tuneData.tl.nLoop,'nrep',tuneData.tl.nRep,'datachan',tuneData.dataChan,'opts','ampok'});
smset(gates,scanvals);
data = smrun(scan); 
smset(gates,0);
if any(isnan(data{1}(:))); tlloc=nan; tlfit = 0; return; end
% Purely empirical fit form
eps = scan.data.pulsegroups.varpar';                      
data=1e3*mean(data{1},1);
if strcmp(side,'left') || strcmp(side,'A')
    figure(11);
else
    figure(13);
end
subplot(2,3,ind);

[~,maxTlInd]=max(data);
epsMax=eps(maxTlInd); % TL is centered around where data largest.
%      1: offset, 2: upward slope, 3: upward center, 4: width both, 5: downward slope, 6: downward center
beta0 = [min(data), range(data), epsMax-range(eps)/6, 150, range(data)/2, epsMax+range(eps)/6,-2e-6];

params=fitwrap('noplot',eps,data,beta0, tuneData.tl.fitFunc);
[params,~,~,~,~,err]=fitwrap('plfit samefig',eps,data,params,tuneData.tl.fitFunc);

fitFunc = str2func(tuneData.tl.fitFunc);
[~,tlInd] = max(fitFunc(params,eps)); 
tlPt = eps(tlInd);
tlerr = err(1,4,1);
% if tlerr > range(eps)/6 
%     tlfitfunc_Err = @(p,x) p(1)+p(2)*(tanh((x-p(3))/p(4))+1)/2 - p(5)*(tanh((x-p(6))/p(4))+1)/2+p(7)*(tanh((x-p(8))/p(9))+1)-p(10)*(tanh((x-p(11))/p(12))+1);
%     beta0 = [min(data), range(data), epsMax-range(eps)/6, range(eps)/12, range(data)/2, epsMax+range(eps)/6 ...
%         ,range(data)/5,eps(end)-range(eps)/5,range(eps)/8,range(data)/5,eps(end)-range(eps)/10,range(eps)/8,];
%     [params,~,~,~,~,err]=fitwrap('plfit samefig',eps,data,beta0,tlfitfunc_Err);
%     tlerr = err(1,4,1);
% end
a = gca; a.YTickLabelRotation=-30;
a.XLim = [min(eps),max(eps)];
tlfit = tlerr/1e3 < 3;
% if params(6) > params(3)+ 2500
%     tlPt = params(3)+200;
%     fprintf('Broad peak. Using left edge + 200.\n')
%     tlfit = 1;
% else
%     tlPt = (params(3)+params(6))/2;
% end
plot(tlPt,fitFunc(params,tlPt),'kx','MarkerSize',8);
plot(params(3),fitFunc(params,params(3)),'x','MarkerSize',8);
plot(params(6),fitFunc(params,params(6)),'x','MarkerSize',8);
plot(mean(eps),fitFunc(params,mean(eps)),'o'); % circle on center point

title(sprintf('%3.f, wid %3.f',tlPt,params(6)-params(3)));
tlloc=tlPt/1000;
if ~tlfit, fprintf('TL didn''t fit. Try retuning and try again \n'); end
end