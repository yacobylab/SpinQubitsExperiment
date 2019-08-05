function [grad,config]=getgradientMulti(grpname,config)
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
global fbdata; global tuneData;
if ~exist('config','var'), config=struct();  end
config=def(config,'chisqthresh',5);
config=def(config,'ampthresh',.2);
config=def(config,'fitwrapopts','');
config = def(config,'opts','both'); 

if nargin == 0, config.figure = 1036; end
if ~exist('grpname','var') || ~isstruct(grpname)
    if isopt(config.opts,'both')
        config = def(config,'datachan','both');    
        ind = 3;
    else
        switch tuneData.activeSetName %guess the data channel
            case 'right'
                config=def(config,'datachan','DAQ2');
                ind = 2;
            case 'left'
                config=def(config,'datachan','DAQ1');
                ind =1;
        end
    end
    config=def(config,'nloop',fbdata.params(ind).nloopfb);
    config=def(config,'opts','nopol');
    % Have to create scan
    if ~exist('grpname','var') || isempty(grpname)
        fbGroup = fbdata.params(ind).fbInit;
    end
    %if ~iscell(config.datachan), config.datachan={config.datachan}; end
    scan = makeFeedbackScan(fbGroup,config.nloop,config.datachan); 
    % add pol back? 
    scan.datachan=config.datachan; % Protect cell arrays.
    if isopt(config.opts,'nodisp'), scan.disp=[]; end    
else
    scan = grpname;
end
nQub = length(scan.loops(1).getchan)-1; 
if nQub > 1, ind = 3; end
%% Run scan (manually)
if ~isopt(config.opts,'reget') % First time, configure DAQ.
    scan = scan.configfn(1).fn(scan,scan.configfn(1).args{:});
end
smatrigfn(1,smchaninst(scan.loops(1).getchan{1}),4); % Arm DAQ. 
for i = 2:length(scan.loops(1).prefn) 
    f = funcify(scan.loops(1).prefn(i).fn);
    f(1,scan.loops(1).prefn(i).args{:});    
end

d=smget(scan.loops(1).getchan(1:nQub));
smset('PulseLine',1,[],'quiet');
rawInd = nQub; histInd = 2*nQub; 
d((1:nQub)+rawInd)=d(1:nQub); % Raw data
for i = 1:nQub        
    d{i+histInd}=histc(d{i},scan.loops(1).procfn(end-nQub+i).fn.args{1})'; % Histogram of data
    d{i}=mean(reshape(d{i},scan.loops(1).procfn(1).fn(3).args{1}),2)'; % Averaged data. 
end
if any(isnan(d{1})), grad=nan; return; end
if ~exist('ind','var')
    if strcmp(config.datachan,'both')
        ind = 3;
    else
        ind = str2double(config.datachan(end)); % determine side
    end
end
for i = 1:nQub
    if isfield(scan.data.pulsegroups(1),'varpar') && ~isempty(scan.data.pulsegroups(1).varpar)
        xv=scan.data.pulsegroups(1).varpar';
    else
        xv=scan.data.pulsegroups(1).pgList{i}.varpar';
    end
    switch fbdata.params(ind).fitType
        case 'fit'
            data_std=std(reshape(d{rawInd+i},length(d{i}),config.nloop),0,2)/sqrt(config.nloop); % Work out x,y, sigma_y
            params=fioscill(xv,squeeze(d{i}),1); % Initial guess for fit
            params(5)=0; params(6) = 0.001;
            cosfn2 = '@(y, x)y(1)+(y(2)*cos(y(4)*x) + y(3) * sin(y(4)*x)).*exp(-(x-y(5)).^2 * y(6).^2)';
            [params,chisq,cov]=mfitwrapfb(struct('x',xv,'y',d{i},'vy',(data_std.^2)'),struct('fn',str2func(cosfn2)),params,'',[1 1 1 1 0 1]); % Actual fit
            grad(i) = 1e3*params(4)/(2*pi);
            gradDev(i) = 1e3*sqrt(cov(4,4))/(2*pi);
            
            vc=scan.loops(1).procfn(end).fn.args{1}; % Histogram voltages
            if size(vc,1)>1, vc = vc'; end
            vc=vc+(vc(2)-vc(1))/2; % Find edges
            n=sum(d{histInd+i}); % Total data points
            xbar = sum(d{histInd+i} .* vc)/n;
            xxbar = sum(d{histInd+i} .* vc .* vc)/n;
            sigstd=sqrt(xxbar-xbar.*xbar); %
            normSig=sqrt(params(2)^2+params(3)^2)/sigstd;
            x = xv; y = d{i};
            if ~isnan(config.ampthresh) && normSig < config.ampthresh, gradDev(i)=1000; end
            if ~isnan(config.chisqthresh) && chisq > config.chisqthresh, gradDev(i)=1000; end
        case 'fft'
            L = length(xv);
            fdata = d{i};
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
            [~,indP] = max(peakMags);
            grad(i) = 1e3*peakFreqs(indP);
            if length(peakFreqs) > 1 % is this the best way to do this?
                gradDev(i)= 1e3*std(peakFreqs);
            else
                [~,indS]=min(abs(fftData(startInd:end)-fftData(peakInds)/2));  % Next closest peak
                gradDev(i) =1e3*freqs(indS);
            end
            %grad=1e3*freqs(peakInd+startInd-1); % +1 because looking at 2:end
            %gradDev = 1e3*(freqs(ind)-freqs(peakInd));
            if peakMags(indP) < 4e-4
                grad(i) = nan;
            end
            x = 1e3*freqs;
            y = fftData;
    end
    if isfield(config,'figure') && config.figure
        if ~isfield(config,'graddata') || ~any(isgraphics(config.graddata)) || length(config.graddata)<i
            figure(config.figure); subplot(2,2,i); cla;
            config.graddata(i) = plot(x,y);
        else
            config.graddata(i).YData = y;
        end
    end
end
config.gradDev = gradDev;
end