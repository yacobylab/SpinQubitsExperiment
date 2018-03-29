function [fileSet,out]=anaEchoNoise(fileSet) 
figStart = 10; 
if ~exist('fileSet','var') || isempty(fileSet)
    fileSet = get_files; 
end
for i = 1:length(fileSet)
    out{i}=mfitEchoFish2(fileSet{i},struct('mfitopts','none','grps',[4 Inf],'fignum',figStart,'xrng',[-29 Inf]));
    figStart = figStart+1; 
    alpha(i) =out{i}.alpha; 
    t2(i) = out{i}.t2;     
    amp(i) = out{i}.amp; 
    j(i) = mean(out{i}.params(1,1:10))/2/pi; 
    eps(i) = out{i}.s.scan.data.pulsegroups(1).params(2);   
    tMax(i) = out{i}.tMax; 
end
djdeps = diff(j)./diff(eps);% jmean = (j(1:end-1)+j(2:end))/2; 
djdepsMean(1) = djdeps(1); 
for i = 2:length(out)-1 
    djdepsMean(i) = 1/2*(djdeps(i-1)+djdeps(i)); 
end
djdepsMean(end+1) = djdeps(end); 
figure(44); clf; 
ha = tight_subplot(2,2); 
axes(ha(1)); 
plot(djdepsMean,t2,'.'); 
xlabel('dJ/d\epsilon'); ylabel('T_2'); 
axes(ha(2)); 
plot(j,t2,'.'); 
xlabel('J'); ylabel('T_2'); 
axes(ha(3)); 
plot(djdepsMean,alpha,'.'); 
xlabel('dJ/d\epsilon'); ylabel('\alpha'); 
axes(ha(4)); 
plot(t2./tMax,'.'); 
ylabel('T_2/tMax'); 
end