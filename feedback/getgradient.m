function [grad,opts]=getgradient(grpname,opts)
% function [grad out]=sm_getgradient(grpname,opts)
%  grpname; name/number of group to use.
%           empty/does not exist: use 2nd group whose name starts with dBz_
%  opts.datachan: default to DAQ2.
%  opts.opts: nopol, nodisp, reget
%  opts.ampthresh: 0.1; If visiblity is less than this, assume fit is bad and make grad_dev huge.
%  opts.chisqthresh: 10; If chi^2 is more than this, assume fit is bad and make grad, grad_dev huge.
%     reget assumes the most recent previous scan to run was sm_getgradient, and skips re-initializing the DAQ card.
% Return the current magnetic field gradient.
% Sign of the gradient is only trustworthy if the gradient is locked.
% Negative means triplet-side.
%    grad_dev: estimate of error bar on gradient measurement

% fill in default options, make the scan
global fbdata; global awgdata; global tuneData;
if ~exist('opts','var'),       opts=struct();  end
opts=def(opts,'chisqthresh',5);
opts=def(opts,'ampthresh',.2);
opts=def(opts,'fitwrapopts','');
opts=def(opts,'nloop',fbdata.nloopfb);
if nargin == 0 
    opts.figure = 1035; 
end
if ~exist('grpname','var') || ~isstruct(grpname)
    switch tuneData.activeSetName %guess the data channel
        case 'right'
            opts=def(opts,'datachan',{'DAQ2'});
        case 'left'
            opts=def(opts,'datachan',{'DAQ1'});
        otherwise
            opts=def(opts,'datachan',{'DAQ2'});
    end
    opts=def(opts,'opts','nopol');    
    if ~exist('grpname','var') || isempty(grpname)
        switch tuneData.activeSetName
            case 'right'
                dbzgrps=find(~cellfun('isempty',regexp({awgdata(1).pulsegroups.name},'^dBz_'))&~cellfun('isempty',regexp({awgdata(1).pulsegroups.name},'_R$')));
            case 'left'
                dbzgrps=find(~cellfun('isempty',regexp({awgdata(1).pulsegroups.name},'^dBz_'))&~cellfun('isempty',regexp({awgdata(1).pulsegroups.name},'_L$')));
            otherwise
                dbzgrps=find(~cellfun('isempty',regexp({awgdata(1).pulsegroups.name},'^dBz_')));
        end
        grpname=dbzgrps(1);% Default to first DBZ group.
    end
    scanfunc = @fConfSeq;
    if ~iscell(opts.datachan), opts.datachan={opts.datachan}; end
    if isstruct(grpname)
        scan=grpname;
    else
        if isopt(opts.opts,'nopol')
            scan=scanfunc(grpname,struct('nloop',opts.nloop,'nrep',1,'opts','raw ampok','datachan',opts.datachan));
        else
            scan=scanfunc(grpname,struct('nloop',opts.nloop,'nrep',1,'opts','pol raw ampok','datachan',opts.datachan));
        end
        scan.datachan=opts.datachan; % Protect cell arrays.
    end
    if isopt(opts.opts,'nodisp')
        scan.disp=[];
    else
        for i=1:length(scan.datachan)
            scan.disp=struct('loop',1,'channel',i,'dim',1);
        end
    end
    scan.configch=[];
else 
    scan = grpname; 
end % Have to create scan 
if ~isopt(opts.opts,'reget') % Take the data
    scan = scan.configfn(1).fn(scan,scan.configfn(1).args{:});
end
smatrigfn(1,smchaninst(scan.loops(1).getchan{1}),4);
for i=2:length(scan.loops(1).prefn)
    if ischar(scan.loops(1).prefn(i).fn)
        f=str2func(scan.loops(1).prefn(i).fn);
        f(1,scan.loops(1).prefn(i).args{:});
    else
        scan.loops(1).prefn(i).fn(1,scan.loops(1).prefn(i).args{:});
    end
end
d=smget(scan.loops(1).getchan{1});
smset('PulseLine',1,[],'quiet');
d{3}=d{1}; d{4}=d{1};
d{1}=mean(reshape(d{1},scan.loops(1).procfn(1).fn(3).args{1}),2)';
d{4}=histc(d{4},scan.loops(1).procfn(4).fn.args{1})';

if any(isnan(d{1})), grad=nan; return; end
xv=scan.data.pulsegroups(1).varpar';
switch fbdata.fitType
    case 'fit'
        data_std=std(reshape(d{3},length(d{1}),opts.nloop),0,2)/sqrt(opts.nloop); % Work out x,y, sigma_y        
        params=fioscill(xv,squeeze(d{1}),1); % Initial guess for fit
        params(5)=0; params(6) = 0.001;
        cosfn2 = '@(y, x)y(1)+(y(2)*cos(y(4)*x) + y(3) * sin(y(4)*x)).*exp(-(x-y(5)).^2 * y(6).^2)';
        [params,chisq,cov]=mfitwrapfb(struct('x',xv,'y',d{1},'vy',(data_std.^2)'),struct('fn',str2func(cosfn2)),params,'',[1 1 1 1 0 1]); % Actual fit
        grad = 1e3*params(4)/(2*pi);
        gradDev = 1e3*sqrt(cov(4,4))/(2*pi);
        
        vc=scan.loops(1).procfn(end).fn.args{1};
        vc=vc+(vc(2)-vc(1))/2;
        n=sum(d{end});
        xbar = sum(d{end} .* vc)/n;
        xxbar = sum(d{end} .* vc .* vc)/n;
        sigstd=sqrt(xxbar-xbar*xbar);
        normSig=sqrt(params(2)^2+params(3)^2)/sigstd;
        x = xv; y = d{1}; 
        if ~isnan(opts.ampthresh) && normSig < opts.ampthresh, gradDev=1000; end
        if ~isnan(opts.chisqthresh) && chisq > opts.chisqthresh, gradDev=1000; end
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
        if peakMags(1)<5e-4  
            grad = nan;
        end
        x = 1e3*freqs;
        y = fftData;
end
opts.gradDev = gradDev;
if isfield(opts,'figure') && opts.figure 
    if ~isfield(opts,'graddata') 
        figure(opts.figure); subplot(2,1,1); cla; 
        opts.graddata = plot(x,y);
    else
        opts.graddata.YData = y; 
    end
end
end