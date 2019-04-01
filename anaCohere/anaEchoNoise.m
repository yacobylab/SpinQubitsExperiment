function [fileSet,out]=anaEchoNoise(fileSet)
% Perform noise analysis on echo data.
% function [fileSet,out]=anaEchoNoise(fileSet)
% Grab files if not given as argument. 
% Perform fit of all excho sets, using "mfitEchoFish" 
% Return the values for beta, 

figStart = 10;
if ~exist('fileSet','var') || isempty(fileSet), fileSet = getFiles; end
for i = 1:length(fileSet)
    out{i}=mfitEchoFish(fileSet{i},struct('mfitopts','none','grps',[4 Inf],'fignum',figStart,'xrng',[-29 Inf]));
    figStart = figStart+2;
    alpha(i) =out{i}.alpha;
    t2(i) = out{i}.t2;
    j(i) = nanmean(out{i}.params(1,1:10))/2/pi;
    eps(i) = out{i}.s.scan.data.pulsegroups(1).params(2);
    tMax(i) = out{i}.tMax;
end
djdeps = diff(j)./diff(eps);% jmean = (j(1:end-1)+j(2:end))/2;
djdepsMean(1) = djdeps(1);
% Fixme: should this be done through fitting? 
for i = 2:length(out)-1
    djdepsMean(i) = 1/2*(djdeps(i-1)+djdeps(i));
end
djdepsMean(end+1) = djdeps(end);
if any(isinf(djdepsMean))
    warning('Some djdepsMean are infinite') 
    djdepsMean(isinf(djdepsMean)) = nan; 
end
%% Plot useful data. 
figure(1111); clf;
ha = tightSubplot([2,2]);

axes(ha(1));
plot(djdepsMean,t2,'.-');
xlabel('dJ/d\epsilon'); ylabel('T_2 (\mu s)');

axes(ha(2));
plot(j*1e3,t2,'.-');
xlabel('J (MHz)'); ylabel('T_2 (\mu s)');

axes(ha(3));
plot(djdepsMean,alpha-1,'.-');
xlabel('dJ/d\epsilon'); ylabel('\beta');

axes(ha(4));
plot(t2./tMax,'.-');
ylabel('T_2/tMax'); xlabel('Group')
formatFig(1111,'exch full',2,2); 
end