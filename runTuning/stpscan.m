function [stploc,stpfit] = stpscan(gates,scanvals,side,ind,opts)
%function stploc = stpscan(gates,scanvals,side,ind) 
%takes an stp scan w/ the offset on gates given in scanvals, then fits to find stp pt. 
% Returns stp location in mV. 
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
scan = fConfSeq(tuneData.stp.plsGrp,{'nloop',tuneData.stp.nLoop,'nrep',tuneData.stp.nRep, 'datachan',tuneData.dataChan,'opts','ampok'});
for i=1:length(gates)
    scan.consts(end+1).setchan=gates{i};
    scan.consts(end).val=scanvals(i);
end
data = smrun(scan);
if any(isnan(data{1}(:))); stploc=nan; return; end
smset(gates,0);
% Purely empirical fit form
data = (mean(data{1},1));
eps = scan.data.pulsegroups.varpar'*1e3;
if strcmp(side,'left') || strcmp(side,'A')
    figure(10);
else
    figure(12);
end
subplot(2,4,ind);
pf=polyfit(eps,data,1); % Remove offset, linear slope
dataLin=smooth(data-pf(1)*eps - pf(2));

ign=5; % points at start and end to ignore.
[maxSTP,maxSTPind]=max(dataLin(ign:end-ign));
% offset, STP height, STP loc, STP width, linslope
beta0 = [pf(2),maxSTP,eps(maxSTPind+ign),range(eps)/8,pf(1)];
mask = [0 1 1 1 0];

params=fitwrap('plfit plinit samefig',eps,data,beta0,tuneData.stp.fitFn,mask);
[params,~,~,~,~,err]=fitwrap('plfit plinit samefig',eps,data,params,tuneData.stp.fitFn);
stperr = err(1,4,1);
stpfit = stperr<5;%4;
stploc=params(3);
stpw=params(4);
title(sprintf('%3.f; wd %3.f',stploc,stpw));
stploc = stploc/1e3; 
end