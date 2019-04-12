function [fidelity,tAdj,STdiff,t1] = anaT1Meas(opts,scan,data,dropTime,bins)
% finds T1, plots measurement histograms, finds optimal meas time.
% function [fidelity,tAdj,STdiff,t1] = anaT1Meas(opts,scan,data,dropTime,bins)
% bins is set to 512 if not given.
% dropTime set to 2.4 us if not given. This is 2 us for manipulation. 0.09 waiting at meas before marker sent, then 310 ns while tank circuit rings up.
% opts: full : fit all timesteps (slow).
%       pl : plot histogram for meas time
%       fpl : produce histogram plots for 4 different times
%       nograd : assume dBz gradient not locked for scan. 
%       fid: plot fidelity as a function of measurement time. 
% Computes T1 and overall, singlet, and triplet fidelities.
% Scan loads singlet and measures (oversampled) and then triplet and measures (oversampled).
% format of data is [S T S T....] (100 of each per row).
% Uses fitting equation in Barthel paper, to find curve showing dist. voltages for singlets 
% and triplets, fidelity found by how much of curve past Vt.
% We may want to update to use the meas dict element to keep offsets up to date.

if ~exist('scan','var') || isempty(scan)
    d=loadFiles('t1'); data = d.data{1}; scan = d.scan;
end
if ~exist('opts','var'), opts = ''; end
pulseLength = abs(scan.data.pulsegroups.zerolen(1));
dt = 1/scan.configfn(1).args{3}(2); %integration time bin
timesteps=1e-9*2*pulseLength/dt; % # time steps in an [S T]. (pulselength given in ns, dt in s).
npoints=size(data,1)*size(data,2);
dataCol = reshape(data',timesteps,npoints/timesteps); % from [S' S'..x num rows;T' T'..;S'...]we get array with each column=[S' T'] (so row# modulo 425 gives time step)
sampNum = 2*size(dataCol,2); %number of S/T runs
if ~exist ('dropTime','var'), dropTime=2.4e-6; end  %mask for getting rid of manipulation time
mask = 1:ceil(dropTime/dt);

singData = dataCol(1:size(dataCol,1)/2,:);
tripData = dataCol(1+size(dataCol,1)/2:size(dataCol,1),:);
singData(mask,:)= [];
tripData(mask,:)= [];
T = dt*(1:size(singData,1)); %array of times of measurements

% Calc T1
diffSig = mean(tripData,2)- mean(singData,2);  %difference between singlet and triplet average voltage as a function of time
fitFcn = @(p,x)p(1)+p(2).*exp(-1*x./p(3)); % p(1) is offset between singlet/triplet scans, p(2) the spacing between peaks, p(3) t1.
beta0 = [0, range(diffSig), 1/abs(range(diffSig))*dt];
params = fitwrap('',T, smooth(diffSig)', beta0, fitFcn, [0 1 1]);
t1=params(3);

%Plot T1 data, skipping last few points.
endStop = 14; T=T(1:end-endStop);
diffSig = diffSig(1:end-endStop);
figure(80); clf; hold on;
set(gcf,'Name','T1 Histograms');
subplot(3,2,3); hold on;
plot(T, diffSig',T, fitFcn(params, dt*(1:length(diffSig))));
a = gca; a.YLim = [min(diffSig),max(diffSig)];
title(sprintf('T_{1} = %.2f \\mus', 1e6*params(3)));

%this makes V(t) the average V from 0 to t for a specific run
%Note: Histogram is made with data averaged from 0 to t
singAve = cumsum(singData,1)./repmat((1:size(singData,1))',1,size(singData,2));
tripAve = cumsum(tripData,1)./repmat((1:size(tripData,1))',1,size(tripData,2));
allData = [singAve tripAve];
cen = mean(allData(:));
rng = 5*std(allData(:));

if ~exist ('bins','var') || isempty(bins), bins=512; end
Vt=linspace(cen-rng,cen+rng,bins); %512 voltage bins with a range of n standard deviations
pfunc=@log;
singHist = histc(singAve', Vt);
singHist = singHist(:,1:end-endStop);
tripHist = histc(tripAve', Vt);
tripHist = tripHist(:,1:end-endStop);

subplot(3,2,1); hold on;
imagesc(1e6*T, Vt,pfunc(singHist)); title('Singlet Histogram');
a = gca; a.XLim = [min(1e6*T),max(1e6*T)];
a.YLim = [min(Vt),max(Vt)];
xlabel('T_{meas}(\mus)'); ylabel('Voltage');

subplot(3,2,2); hold on;
imagesc(1e6*T, Vt,pfunc(tripHist)); title('Triplet Histogram');
a = gca; a.XLim = [min(1e6*T),max(1e6*T)];
a.YLim = [min(Vt),max(Vt)];
xlabel('T_{meas}(\mus)'); ylabel('Voltage');

%Now, find V_thresh and T_meas by fitting.
hist=singHist+tripHist;
if isopt(opts,'nograd'), hist = singHist; end % Only take one histogram.
hist = hist';
if isopt(opts,'full') % Fit each data step
    short=1;
else % Fit every 4th data step.
    short=4;
end
timestep=floor(length(T)/short);

% params: 1: t/T1, 2: noise/peak spacing x: voltage]
% See Barthel 2010 single shot read out paper for equation. This gives function for the gaussian + decay function.
% Gaussian has form % exp[(-v-vt)^2/2sig^2].
% Decaying part finds probability of average voltage v (based on time of decay), constructs gaussian centered at v, and integrates over v.
% erf(x) = 2/sqrt(pi) int_0^x exp(-t^2) dt,  % 1/(sqrt(2pi)/sig e^(-t/t1) e^(-v-1)^2/(2 sig^2)+t/2t1 e^(t/2t1 *(t/t1* vNoise/STDiff-2 v)
distfn = @(a, x) exp(-a(1)) * exp(-(x-1).^2./(2* a(2)^2))/sqrt(2*pi)./a(2) + a(1)/2 * exp(a(1)/2 * (a(1) * a(2)^2 - 2 * x)) ...
    .* (erf((1 + a(1) * a(2)^2 - x)/sqrt(2) ./ a(2)) + erf((-a(1) * a(2)^2 + x)/sqrt(2) ./ a(2)));
fitfn = @(p, x) p(3) * distfn(abs(p([5, 7])), .5-(x-p(1)).*p(2)) + p(4) * distfn(abs(p([6, 7])), .5+(x-p(1)).*p(2));
for i=1:timestep
    t=short*i;
    histCurr = hist(t,:);
    aveV = sum(histCurr.*Vt)/sum(histCurr);
    sdV = sqrt(sum(histCurr.*Vt.^2)/sum(histCurr)-aveV^2);
    %1: mean V, 2: 1/peak spacing, 3: left peak mag 4: right peak mag 5: t/T1S 6: t/T1T, 7: noise/peak spacing
    beta0=[aveV, 1/2/sdV, 0.6*max(histCurr), .4*max(histCurr), t*dt*1e-4, t*dt/t1, 2.5/sqrt(t)];
    fitpar(i,:)=fitwrap('',Vt,histCurr,beta0,fitfn,[1 1 1 1 0 0 1]); %#ok<*AGROW> % Don't fit the t1s.
    fitparCurr=fitpar(i,:);
    
    Sfit=fitparCurr; Sfit(3)=1; Sfit(4)=0; % Only keep singlet peak.
    fitSing=(fitparCurr(3)+fitparCurr(4))*fitfn(Sfit,Vt)/sampNum; % Fitted hist of just sing peak
    SfidCurr=cumsum(fitSing); SfidCurr=[0,SfidCurr]; %Sum cumulative prob of capturing all singlets as a function of voltage
    Sfid(:,i)=SfidCurr;
    Tfit=fitparCurr; Tfit(3)=0; Tfit(4)=1; % Only keep triplet peak
    fitTrip=(fitparCurr(3)+fitparCurr(4))*fitfn(Tfit,Vt)/sampNum; % Fitted hist of just trip peak.
    TfidRev=cumsum(fitTrip(end:-1:1)); TfidRev=[0 TfidRev];
    Tfid(:,i)=TfidRev(end:-1:1); %Sum cumulative prob of capturing all singlets as a function of voltage
end
STdiffVec=1./fitpar(:,2); STdiff=median(STdiffVec);
fprintf('Peak Separation: %2.3g mV. ',STdiff*1e3);

fidArr=(Sfid+Tfid)/2; % Fid = (correctly identified/total);
[maxFidVec,VtInd] = max(fidArr); %VtFit is optimal Vthreshold ind
[fidelity,TmeasInd] = max(maxFidVec);   % Tmeasfit gives the index of the time at which there is max fidelity. use this also to find the vT at that time. Note that since we don't fit all data points, cannot just apply to T.
tMeasStep=TmeasInd*short; tMeas=dt*tMeasStep;
tAdj=tMeas+dropTime-2e-6+0.15e-6; % Recommended measurement time is opt fid, ignoring manip time,adding on final off time.

subplot(3,2,4);
plotter(hist(tMeasStep,:),Vt, sampNum,tMeas,fitfn,fitpar(TmeasInd,:));
subplot(3,2,5);
plot(1e6*dt*short*(1:timestep),maxFidVec);
xlabel('T (\mus)'); title('Best Fidelity');

if isopt(opts,'fpl') % Plot histograms for 4 times.
    figure(81); clf;
    tList = [16,52,160,300];
    for i=1:4
        subplot(2,2,i); hold on;
        plotter(hist(tList(i),:),Vt, sampNum,dt*tList(i),fitfn,fitpar(tList(i)/short,:));
    end
end
fprintf('Fidelity = %2.3f.  Tmeas = %.2f usec. Vthreshold = %2.3f \n', 100*fidelity, tAdj*1e6, 1e3*Vt(VtInd(TmeasInd)));
fprintf('T1 = %2.3g us, F_s = %2.2f, F_t = %2.2f \n',1e6*t1,100*Sfid(VtInd(TmeasInd),TmeasInd),100*Tfid(VtInd(TmeasInd),TmeasInd));

if isopt(opts,'fid')
    figure(91); clf; hold on;
    inds=sub2ind(size(fidArr),VtInd,1:length(VtInd));
    plot(1e6*tAdj,1-Sfid(inds),'.-')
    plot(1e6*tAdj,1-Tfid(inds),'.-')
    plot(1e6*tAdj,1-maxFidVec,'.-')
    xlabel('Time (\mu s)'); ylabel('1-Fidelity')
    title(sprintf('T_1 %g \mu s, V_{sep} %g mV',1e6*t1,1e3*STdiff))
end
end

function plotter(hist,Vt , sampNum,tm,fitfn,fitpar)
Sfit=fitpar; Sfit(3)=1; Sfit(4)=0;
singFit=(fitpar(3)+fitpar(4))*fitfn(Sfit,Vt)/sampNum;
Tfit=fitpar; Tfit(3)=0; Tfit(4)=1;
tripFit=(fitpar(3)+fitpar(4))*fitfn(Tfit,Vt)/sampNum;

npoints = 1;%sum(hist(t,:));
plot(Vt,sampNum*fitpar(3)/(fitpar(3)+fitpar(4))*singFit/npoints,'LineWidth',2); hold on;
plot(Vt,sampNum*fitpar(4)/(fitpar(3)+fitpar(4))*tripFit/npoints,'LineWidth',2);
plot(Vt,smooth(hist)/npoints,'.');
plot(Vt,fitfn(fitpar,Vt)/npoints,'LineWidth',1);
a = gca; a.XLim = [min(Vt),max(Vt)];
xlabel('Voltage (V)'); ylabel('Prob');
title(sprintf('t=%.2f us', tm*1e6))
end