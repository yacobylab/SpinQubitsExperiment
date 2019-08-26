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
    fitopts = 'residuals colorplot nocenter plotfit fishless';
end
f = [10,11]; 
pptControl('start')

for i = 1:length(fileSet)
    if ~isopt(opts, 'white')
       out(i)=mfitEchoFish(fileSet{i},struct('mfitopts','none','grps',[2 Inf],'fignum',figStart,'opts',fitopts));
    else
        out(i)=mfitEchoFishWhite(fileSet{i},struct('mfitopts','none','grps',[4 Inf],'fignum',figStart,'xrng',[-27 Inf],'opts',fitopts));
    end
    %figStart = figStart+2;
    if isempty(out) || ~isfield(out,'s'), continue; end    
    jList(i,:) = out(i).J;    
    jErr(i,:) = out(i).jErr;
    j(i)=mlePar(jList(i,:),jErr(i,:)); 
    prettyFile = removePath(out(i).file); 
    str = sprintf('File %s. \n alpha = %3.3g. T_2 = %3.3g us. Amp = %3.3g. Fish Amp = %3.3g. T_2* = %3.3g ns', prettyFile,out(i).alpha, out(i).t2, out(i).amp,out(i).fishAmp,out(i).t2s);
    str2 = sprintf('J = %3.3g MHz, T1 = %3.3g us. Peak separation %3.3g mV. Fidelity %3.3f%%.',j(i),out(i).s.t1*1e6,1e3*diff(out(i).s.meanvals),out(i).s.fidelity*100);
    slideInfo.body = str;
    slideInfo.body2 = str2;
    slideInfo.comments = ''; slideInfo.title = '';    
    save2pptauto(slideInfo,f)
    eps(i) = out(i).s.scan.data.pulsegroups(1).params(2);        
    dBz(i) = out(i).s.scan.data.setpt(1);     
end
params.eps = eps;
params.J = j; 
%badJ = jList-nanmean(jList) > 2*std(jList); 
%jList(badJ) = nan; 
%j = nanmean(jList,2); 
alpha =[out.alpha];
tMax = [out.tMax];
t2s = [out.t2s];
chisq=[out.chisq]; 
fishAmp = [out.fishAmp];    
amp = [out.amp];     
t2 = [out.t2];

jNorm=sqrt(j.^2-dBz.^2); %subtract of dBz value. 

beta0=[35 3500 .2];
fitFn = @(p,x) p(1)+p(2)*exp(-x/p(3));
params.jPars=fitwrap('plinit plfit',eps,j,beta0,fitFn); % fit j vs. eps

djdeps=params.jPars(2)/params.jPars(3)*exp(-eps/params.jPars(3)); %djdeps

djdepsDiff = diff(j)./diff(eps);% jmean = (j(1:end-1)+j(2:end))/2;
djdepsMean(1) = djdepsDiff(1);

for i = 2:length(out)-1
    djdepsMean(i) = 1/2*(djdepsDiff(i-1)+djdepsDiff(i));
end
 djdepsMean(end+1) = djdepsDiff(end);
 djdepsMean = abs(djdepsMean); 
% if any(isinf(djdepsMean))
%     warning('Some djdepsMean are infinite, probably using eps of same value.') 
%     djdepsMean(isinf(djdepsMean)) = nan; 
% end
%% Plot useful data. 
figure(1111); clf;
ha = tightSubplot([3,3]);

axes(ha(1));
plot(eps,djdeps,'.-'); hold on; 
plot(eps,djdepsMean,'.-'); 
xlabel('\epsilon (mV)'); ylabel('dJ/d\epsilon');

axes(ha(2));
errorbar(j,t2,[out.t2Err],'.-'); 
xlabel('J (MHz)'); ylabel('T_2 (\mus)');

axes(ha(3));
errorbar(j,alpha-1,[out.alphaErr],'.-');
xlabel('J (MHz)');
%xlabel('dJ/d\epsilon'); 
ylabel('\beta');

axes(ha(4));
plot(t2./tMax,'.-');
ylabel('T_2/tMax'); xlabel('Group')

axes(ha(5)); 
plot(j,abs(amp),'.-'); hold on; 
plot(j,abs(fishAmp),'.-'); 
xlabel('J (MHz');
ylabel('Amplitude');

axes(ha(6)); 
%errorbar(t2s,[out.t2sErr],'.-'); 
errorbar(j,t2s,[out.t2sErr],'.-'); 
xlabel('J (MHz'); 
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
xlabel('J (MHz)'); ylabel('Q'); 
formatFig(1111,'exch full',3,3); 
%% Fit the noise 
% t2 in us, djdeps in MHz / uV 
inds = 1:length(t2s); 
figure(1112); clf; 
ga = tightSubplot([2,2],'title'); 
%ga = tightSubplot([2,2]); 
alphaErr = [out.alphaErr]; 
alphaMean = mlePar(alpha(inds),alphaErr(inds)); 
beta = alphaMean-1; 
gval=gamma(-1-beta)*sin(pi*beta/2)*2^-beta*(-2+2^beta);

scaledT2=(1./t2).^alphaMean;
scaleFactor = 1e-3;
djdepsSc = djdeps*scaleFactor; 
% Below is just standard error propagation 
scaleT2Err=sqrt(([out.alphaErr] .* ((1./t2).^alpha) .* log(1./t2)).^2 + ([out.t2Err] .*alpha .* ((1./t2).^(alpha+1))).^2);
fitfn = @(p,x) p(1)*x+p(2);

% Start by trying nlinfit
axes(ga(1)); 
fpp = fitwrap('plfit samefig',djdepsSc(inds).^2,scaledT2(inds),[1 1], fitfn,[1 1]);
xlabel('(dJ/d\epsilon)^2, MHz/\muV');
ylabel('T_2^{(-\alpha)}');
% 1e12 converts from MHz/uV to Hz/V . %The first part deals with units
params.Se0 = 1e-18*1e6^(beta)*fpp(1)*(2*pi)^(-beta)/abs(gval); 
params.Seps0 = @(f) params.Se0/f^beta;
title(sprintf('%1.2f nV/sq{Hz} at 1 MHz. beta = %1.2f', 1e9*sqrt(params.Seps0(1e6)),beta)); 

% Then do mlefit (uses error bars on params, should be more accurate)
data.x = djdepsSc(inds).^2; 
data.y = scaledT2(inds); 
data.vary = scaleT2Err(inds); 
model.fn = fitfn; 
[params.noiseFit,chisq,cov]=mfitwrap(data,model,[1,1],'none');
axes(ga(2)); 

% 1e12 converts from MHz/uV to Hz/V . %The first part deals with units
params.Se = 1e-18*1e6^(beta)*params.noiseFit(1)*(2*pi)^(-beta)/abs(gval); 
params.Seps = @(f) params.Se/f^beta;

errorbar(djdepsSc(inds).^2,scaledT2(inds),scaleT2Err(inds),'.-'); hold on;
plot(djdepsSc(inds).^2,params.noiseFit(1)*djdepsSc(inds).^2+params.noiseFit(2))
xlabel('(dJ/d\epsilon)^2, MHz/\muV');
ylabel('T_2^{(-\alpha)}');
title(sprintf('%1.2f nV/sq{Hz} at 1 MHz. beta = %1.2f', 1e9*sqrt(params.Seps(1e6)),beta)); 
%% Now fit low frequency noise using T2* data
axes(ga(3)); hold on;
inds2 = 1:length(t2s); 
fitfn = @(p,x) sqrt(p(2) * x.^2 + p(1));
beta0 = [10 150];
[p,~,~,~,mseFit] = fitwrap('plfit samefig',djdeps(inds2)/1e3,1/sqrt(2)/pi./(1e-3*t2s(inds2)),beta0,fitfn);
title(sprintf('\\epsilon_{RMS} =%3.3g (\\muV), \\sigma_t = %3.3g MHz',sqrt(p(2)),sqrt(p(1))));
xlabel('(dJ/d\epsilon) (MHz/\muV)'); ylabel('\sigma (MHz)');
save2pptauto(slideInfo,[1111, 1112])
%pptControl('save'); 
end