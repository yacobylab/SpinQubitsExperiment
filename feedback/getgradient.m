function [grad,config]=getgradient(grpname,config)
% Measure magnetic field gradient using fitting or fft.
% To speed things up, we run the scan fcns manually, not the scan itself. 
% function [grad, opts]=getgradient(grpname,opts)
% grpname; name/number of group to use.
%          empty/does not exist: use 2nd group on AWG whose name starts with dBz_
% config fields: 
%   datachan: default to DAQ2.
%   opts: nopol, nodisp, reget
%       nopol:
%       nodisp:
%       reget: assumes the most recent previous scan to run was sm_getgradient, and skips re-initializing the DAQ card.
% 	ampthresh: 0.1; If visiblity is less than this, assume fit is bad and make grad_dev huge.
% 	chisqthresh: 5; If chi^2 is more than this, assume fit is bad and make grad, grad_dev huge.
% Sign of the gradient is only trustworthy if the gradient is locked.
% Negative means triplet-side.
%    grad_dev: estimate of error bar on gradient measurement
% Some things in this function a little awkward to make things go faster.

%% fill in default options, make the scan
global fbdata; global awgdata; global tuneData;
if ~exist('config','var'), config=struct();  end
config=def(config,'chisqthresh',5);
config=def(config,'ampthresh',.2);
config=def(config,'fitwrapopts','');

if nargin == 0, config.figure = 1035; end
if ~exist('grpname','var') || ~isstruct(grpname)
    switch tuneData.activeSetName %guess the data channel
        case 'right'
            config=def(config,'datachan','DAQ2');
            ind = 2;
        case 'left'
            config=def(config,'datachan','DAQ1');
            ind =1;        
    end
    config=def(config,'nloop',fbdata.params(ind).nloopfb);
    config=def(config,'opts','nopol');
    % Have to create scan
    if ~exist('grpname','var') || isempty(grpname)
        switch tuneData.activeSetName
            case 'right'
                dbzgrps=find(~cellfun('isempty',regexp({awgdata(1).pulsegroups.name},'^dBz_'))&~cellfun('isempty',regexp({awgdata(1).pulsegroups.name},'_R$')));
            case 'left'
                dbzgrps=find(~cellfun('isempty',regexp({awgdata(1).pulsegroups.name},'^dBz_'))&~cellfun('isempty',regexp({awgdata(1).pulsegroups.name},'_L$')));
            otherwise
                dbzgrps=find(~cellfun('isempty',regexp({awgdata(1).pulsegroups.name},'^dBz_')));
        end
        grpname=dbzgrps(1);% Default to first dBz group.
    end
    %if ~iscell(config.datachan), config.datachan={config.datachan}; end
    scan = makeFeedbackScan(grpname,config.nloop,config.datachan); 
    % add pol back? 
    scan.datachan=config.datachan; % Protect cell arrays.
    if isopt(config.opts,'nodisp'), scan.disp=[]; end    
else
    scan = grpname;
end
%% Run scan (manually)
if ~isopt(config.opts,'reget') % First time, configure DAQ.
    scan = scan.configfn(1).fn(scan,scan.configfn(1).args{:});
end
smatrigfn(1,smchaninst(scan.loops(1).getchan{1}),4); % Arm DAQ. 
for i=2:length(scan.loops(1).prefn) 
    f = funcify(scan.loops(1).prefn(i).fn);
    f(1,scan.loops(1).prefn(i).args{:});    
end
d=smget(scan.loops(1).getchan{1});
smset('PulseLine',1,[],'quiet');
d{3}=d{1};
d{4}=histc(d{1},scan.loops(1).procfn(4).fn.args{1})';
d{1}=mean(reshape(d{1},scan.loops(1).procfn(1).fn(3).args{1}),2)';

if any(isnan(d{1})), grad=nan; return; end
xv=scan.data.pulsegroups(1).varpar';
ind = str2double(config.datachan(end)); % determine side
switch fbdata.params(ind).fitType
    case 'fit'
        data_std=std(reshape(d{3},length(d{1}),config.nloop),0,2)/sqrt(config.nloop); % Work out x,y, sigma_y
        params=fioscill(xv,squeeze(d{1}),1); % Initial guess for fit
        params(5)=0; params(6) = 0.001;
        cosfn2 = '@(y, x)y(1)+(y(2)*cos(y(4)*x) + y(3) * sin(y(4)*x)).*exp(-(x-y(5)).^2 * y(6).^2)';
        [params,chisq,cov]=mfitwrapfb(struct('x',xv,'y',d{1},'vy',(data_std.^2)'),struct('fn',str2func(cosfn2)),params,'',[1 1 1 1 0 1]); % Actual fit
        grad = 1e3*params(4)/(2*pi);
        gradDev = 1e3*sqrt(cov(4,4))/(2*pi);
        
        vc=scan.loops(1).procfn(end).fn.args{1}; % Histogram voltages
        if size(vc,1)>1, vc = vc'; end
        vc=vc+(vc(2)-vc(1))/2; % Find edges
        n=sum(d{end}); % Total data points
        xbar = sum(d{4} .* vc)/n;
        xxbar = sum(d{4} .* vc .* vc)/n;
        sigstd=sqrt(xxbar-xbar.*xbar); %
        normSig=sqrt(params(2)^2+params(3)^2)/sigstd;
        x = xv; y = d{1};
        if ~isnan(config.ampthresh) && normSig < config.ampthresh, gradDev=1000; end
        if ~isnan(config.chisqthresh) && chisq > config.chisqthresh, gradDev=1000; end
    case 'fft'
        L = length(xv);
        fdata = d{1};
        T_samp = abs(diff(xv(1:2))); F_samp = 1/T_samp;
        NFFT = 2^nextpow2(25*L); % Padded with zeros
        nyq = abs(F_samp/2); % Nyquist frequency
        nAlias = 0; % Number of times signal is aliased
        freqs = (nAlias+.5)*nyq+(-1)^nAlias*(F_samp/2)*linspace(-.5,.5,NFFT/2+1);
        fftData = fft(fdata(1,:)-mean(fdata(1,:)),NFFT)/L;
        fftData = 2*abs(fftData(1:NFFT/2+1));
        startInd=3;
        %[val, peakInd] = max(fftData(startInd:goodInds)); % Largest fft peak
        peakInds=peakfinder(fftData);
        peakFreqs = freqs(peakInds);
        peakMags = fftData(peakInds);
        grad = 1e3*peakFreqs(1);
        if length(peakFreqs)>1
            gradDev= 1e3*std(peakFreqs);
        else
            [~,ind]=min(abs(fftData(startInd:end)-fftData(peakInds)/2));  % Next closest peak
            gradDev =1e3*freqs(ind);
        end
        %grad=1e3*freqs(peakInd+startInd-1); % +1 because looking at 2:end
        %gradDev = 1e3*(freqs(ind)-freqs(peakInd));
        if peakMags(1) < 4e-4
            grad = nan;
        end
        x = 1e3*freqs;
        y = fftData;
end
config.gradDev = gradDev;
if isfield(config,'figure') && config.figure
    if ~isfield(config,'graddata')
        figure(config.figure); subplot(2,1,1); cla;
        config.graddata = plot(x,y);
    else
        config.graddata.YData = y;
    end
end
end