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
for i=1:length(gates)
    scan.consts(end+1).setchan=gates{i};
    scan.consts(end).val=scanvals(i);
end
data = smrun(scan); smset(gates,0); 
if any(isnan(data{1}(:))); tlloc=nan; return; end
% Purely empirical fit form
eps = scan.data.pulsegroups.varpar';                      
data=mean(data{1},1);
[~,ci]=max(data);
epsMax=eps(ci);
if strcmp(side,'left') || strcmp(side,'A')
    figure(11);
else
    figure(13);
end
subplot(2,4,ind);
beta0 = [min(data), range(data), epsMax-range(eps)/6, range(eps)/12, range(data)/2, epsMax+range(eps)/6];
params=fitwrap('plinit plfit samefig',eps,data,beta0, tuneData.tl.fitFunc);
[params,~,~,~,~,err]=fitwrap('plinit plfit samefig',eps,data,params,tuneData.tl.fitFunc);
tlerr = err(1,4,1);
if tlerr > range(eps)/6 
    tlfitfunc_Err = @(p,x) p(1)+p(2)*(tanh((x-p(3))/p(4))+1)/2 - p(5)*(tanh((x-p(6))/p(4))+1)/2+p(7)*(tanh((x-p(8))/p(9))+1)-p(10)*(tanh((x-p(11))/p(12))+1);
    beta0 = [min(data), range(data), epsMax-range(eps)/6, range(eps)/12, range(data)/2, epsMax+range(eps)/6 ...
        ,range(data)/5,eps(end)-range(eps)/5,range(eps)/8,range(data)/5,eps(end)-range(eps)/10,range(eps)/8,];
    [params,~,~,~,~,err]=fitwrap('plinit plfit samefig',eps,data,beta0,tlfitfunc_Err);
    tlerr = err(1,4,1);
end
tlfit = tlerr/1e3 < 3;
if params(6) > params(3)+ 800
    tl_de = params(3)+200;
    fprintf('Broad peak. Using left edge + 200.\n')
    tlfit = 1;
else
    tl_de = (params(3)+params(6))/2;
end
title(sprintf('%3.f',tl_de));
tlloc=tl_de/1000;
end