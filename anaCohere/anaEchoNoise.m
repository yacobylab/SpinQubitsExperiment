function [fileSet,out,params]=anaEchoNoise(fileSet,opts)
% Perform noise analysis on echo data.
% function [fileSet,out]=anaEchoNoise(fileSet)
% Grab files if not given as argument. 
% Perform fit of all excho sets, using "mfitEchoFish" 
% Return the values for beta, 
if ~exist('opts','var'), opts = ''; end
figStart = 10;
if ~exist('fileSet','var') || isempty(fileSet), fileSet = getFiles; end
if isopt(opts,'noplot') 
    fitopts = 'nocenter'; 
else
    fitopts = 'residuals colorplot nocenter plotfit';
end
for i = 1:length(fileSet)
    if ~isopt(opts, 'white')
       out(i)=mfitEchoFish(fileSet{i},struct('mfitopts','none','grps',[4 Inf],'fignum',figStart,'xrng',[-29 Inf],'opts',fitopts));
    else
       out(i)=mfitEchoFishWhite(fileSet{i},struct('mfitopts','none','grps',[4 Inf],'fignum',figStart,'xrng',[-29 Inf],'opts',fitopts));        
    end
    figStart = figStart+2;        
    jList(i,:) = out(i).J;
    jErr(i,:) = out(i).jErr;
    eps(i) = out(i).s.scan.data.pulsegroups(1).params(2);        
    dBz(i) = out(i).s.scan.data.setpt(1); 
end
%badJ = jList-nanmean(jList) > 2*std(jList); 
%jList(badJ) = nan; 
%j = nanmean(jList,2); 
alpha =[out(i).alpha];
tMax = [out.tMax];
t2s = [out.t2s];
chisq=[out.chisq]; 
fishAmp = [out.fishAmp];    
amp = [out.amp];     
t2 = [out.t2];
j=mleParFit(jList,jErr); 
j=j';
jNorm=sqrt(j.^2-dBz.^2); %subtract of dBz value. 

beta0=[0 3500 .1];
fitFn = @(p,x) p(1)+p(2)*exp(-x/p(3));
params.jPars=fitwrap('plinit plfit',eps,jNorm,beta0,fitFn); % fit j vs. eps

djdeps=params.jPars(2)/params.jPars(3)*exp(-eps/params.jPars(3)); %djdeps

djdepsDiff = diff(j)./diff(eps);% jmean = (j(1:end-1)+j(2:end))/2;
djdepsMean(1) = djdepsDiff(1);

for i = 2:length(out)-1
    djdepsMean(i) = 1/2*(djdepsDiff(i-1)+djdepsDiff(i));
end
 djdepsMean(end+1) = djdepsDiff(end);
% if any(isinf(djdepsMean))
%     warning('Some djdepsMean are infinite, probably using eps of same value.') 
%     djdepsMean(isinf(djdepsMean)) = nan; 
% end
%% Plot useful data. 
figure(1111); clf;
ha = tightSubplot([3,3]);

axes(ha(1));
plot(eps,djdeps,'.-'); hold on; 
plot(eps,-djdepsMean,'.-'); 
xlabel('\epsilon (mV)'); ylabel('dJ/d\epsilon');

axes(ha(2));
errorbar(j,t2,[out.t2Err],'.-'); 
xlabel('J (MHz)'); ylabel('T_2 (\mu s)');

axes(ha(3));
errorbar(alpha-1,[out.alphaErr],'.-');
%xlabel('dJ/d\epsilon'); 
ylabel('\beta');

axes(ha(4));
plot(t2./tMax,'.-');
ylabel('T_2/tMax'); xlabel('Group')

axes(ha(5)); 
plot(j,amp,'.-'); hold on; 
plot(j,fishAmp,'.-'); 
xlabel('J (MHz');
ylabel('Amplitude');

axes(ha(6)); 
errorbar(t2s,[out.t2sErr],'.-'); 
%xlabel('J (MHz'); 
ylabel('T_2^* (ns)');

axes(ha(7)); 
plot(chisq,'.-'); 
ylabel('\chi^2'); 

axes(ha(8)); 
semilogy(eps,j,'.-'); hold on; 
semilogy(eps,jNorm,'.-'); 
xlabel('\epsilon (mV)'); ylabel('J (MHz)'); 

axes(ha(9)); 
Q = t2 .* j; 
params.Q = Q; 
plot(j,Q,'.-'); 
formatFig(1111,'exch full',3,3); 

%% Fit the noise 
% t2 in us, djdeps in MHz / uV 
figure(1112); clf;
ga = tightSubplot([2,2]); 
scaledT2=(1./t2).^mean(alpha);
scaleFactor = 1e-3;
djdepsmu = djdeps*scaleFactor; 
scaleT2Err=sqrt(([out.alphaErr] .* ((1./t2).^alpha) .* log(1./t2)).^2 + ([out.t2Err] .*alpha .* ((1./t2).^(alpha+1))).^2);
fitfn = @(p,x) p(1)*x+p(2); 
axes(ga(1)); 
fpp = fitwrap('plfit',djdepsmu.^2,scaledT2,[1 1], fitfn,[1 1]);
fpp = fitwrap('plfit samefig',djdepsmu.^2,scaledT2,fpp, fitfn,[1 1]);
% we could do an mfitwrap 

data.x = djdepsmu.^2; 
data.y = scaledT2; 
data.vary = scaleT2Err; 
model.fn = fitfn; 
[params.noiseFit,chisq,cov]=mfitwrap(data,model,[1,1]);
 
errorbar(djdepsmu.^2,scaledT2,scaleT2Err,'LineStyle','None'); hold on;
xlabel('(dJ/d\epsilon)^2, MHz/\muV');
ylabel('T_2^{(-\alpha)}');
title(sprintf('%g x + %g',fpp))

axes(ga(2)); 
plot(djdepsmu,sqrt(scaledT2),'.-'); hold on;
beta = nanmean(alpha-1);
gval=gamma(-1-beta)*sin(pi*beta/2)*2^-beta*(-2+2^beta);
params.Se = 1e-18*1e6^(beta)*noiseFit(1)*(2*pi)^(-beta)/abs(gval); %1e12 converts from MHz/uV to Hz/V . %The first part deals with units
% In limit as beta->0, this is 1/(4*pi^2)      

params.Seps = @(f) params.Se/f^beta;
title(sprintf('%f nV/sqrt{Hz} at 1 MHz', 1e9*sqrt(params.Seps(1e6)))); 

end