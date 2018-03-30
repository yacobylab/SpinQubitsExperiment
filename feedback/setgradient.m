function [good,pumpTime,pumpHist]=setgradient(~,fbScan,tgt,config)
%  Lock the gradient to tgt (default:setpt given in fbdata.params.setpt)
% function [good,pumpTime,pumpHist=setgradient(~,fbScan,tgt,opts)
%  config:
%      grpname   Group name or number for dBz measurement
%                   default: first group that starts with dBz_
%     pumptime   Time to pump when trying to lock.  defined in fbdata.params.
%     attempts   How many times to try.  Default 100
%       updown   fbdata.buttonpls indices for up, down pump (singlet, triplet).  Default [4 3]
%       figure   Figure to use for status display.  Default 1035
%   opts: 
% The heart of this code is a Kalman filter.  See wikipedia page for similar notation.

global fbdata; global tuneData;
params = fbdata.params; 
flipcount=0; good = 0;
if ~exist('config','var'), config=struct(); end
config=def(config,'datachan',tuneData.dataChan);
config=def(config,'grpname','');
config=def(config,'pumptime',[]);
config=def(config,'attempts',50);
config=def(config,'opts','');
config=def(config,'gopts','nopol nodisp');
config=def(config,'figure',1035);

ind = str2double(config.datachan(end));
config=def(config,'tol',fbdata.params(ind).tolInit); 
config=def(config,'dir',fbdata.params(ind).dir); 
switch tuneData.activeSetName
    case 'right'
        config=def(config,'updown',[4 3]);  % triplet then singlet pump
    case 'left'
        config=def(config,'updown',[2 1]);  % triplet then singlet pump
    otherwise
        config=def(config,'updown',[4 3]);  % triplet then singlet pump
end

if ~exist('tgt','var') || isempty(tgt), tgt=params(ind).setpt; end
if isempty(config.pumptime), config.pumptime=fbdata.pumptime(ind); end
gradOpts=config; gradOpts.opts=gradOpts.gopts; gradOpts = rmfield(gradOpts,'gopts'); 
xHist=[]; obsHist=[]; pumpHist=[];
if ~exist('fbScan','var') || isempty(fbScan)
    fbGroup = params(ind).fbInit;
    if isopt(config.opts,'fine')
        fbGroup = params(ind).fbGroup;
    end
    fbScan=fConfSeq(fbGroup,struct('nloop',fbdata.nloopfb,'nrep',1,'opts','raw','datachan',tuneData.dataChan)); %make the feedback scan: this will run as a prefunction
    fbScan.configch=[];
    fbScan.figure=1111;
    fbScan.loops(1).setchan='count2';
    fbScan.disp=fbScan.disp([fbScan.disp.dim]==1);
    fbScan.consts(1) = []; % remove clock setting.
    fbScan.xv=fbScan.data.pulsegroups(1).varpar(:,1)';
end
[grad,gradOpts] = getgradient(fbScan,gradOpts); % measure gradient.
gradOpts.opts = [gradOpts.opts 'reget']; % after first time, faster to not run whole scan.
x=[grad; params(ind).prate]; % state of filter: gradient, singlet rate, triplet rate. 
if isfield(params,'pMatrix') && ~isempty(params(ind).pMatrix)
    P = params(ind).pMatrix; 
else
    P=diag([gradOpts.gradDev^2, 1,1]);  % estimated covariance matrix.  We start with no knowledge of pump rates.
end
% The gradient fluctuates a lot, but the pump rate should be pretty stable.
Q=diag([5 1 1]); % Process noise, half-assed guess.
j=1;
err(j)=grad-params(ind).setpt;
H = [1 0 0]; % Perform a correction step.
maxAccel=20;
pumpTime = [];
if isopt(config.opts,'long'), config.attempts=500; end
nullTime = 2.5e-3; nullOffTime = 2.5e-3;
while (abs(err(j))>config.tol || isnan(err(j))) && j < config.attempts 
    j = j+1;
    if strcmp(config.dir,'sing') && grad < tgt % this tries to pump gradient on singlet side. To change, change this, change F, pumpHIst. 
        pulseLine=params(ind).singletPulse;
        F=[1 1 0; 0 1 0 ; 0 0 1]; % Matrix for change in x. 
        pumpHist=[pumpHist 1]; %#ok<*AGROW>
        prate=x(2);
    elseif strcmp(config.dir,'sing') && grad > tgt
        pulseLine = params(ind).tripletPulse;
        F=[ 1 0 -1 ; 0 1 0 ; 0 0 1];
        pumpHist=[pumpHist -1];
        prate=x(3);
    elseif strcmp(config.dir,'trip') && grad < tgt % Pump triplet. 
        pulseLine=params(ind).tripletPulse;
        F=[1 0 1; 0 1 0 ; 0 0 1]; % Matrix for change in x. 
        pumpHist=[pumpHist 1]; %#ok<*AGROW>
        prate=x(3);
    else
        pulseLine = params(ind).singletPulse;
        F=[ 1 -1 0 ; 0 1 0 ; 0 0 1];
        pumpHist=[pumpHist -1];
        prate=x(2);
    end
    % F represents change matrix. So we expect pump rates to be the same, gradient to change re pump time. 
    idealTime=abs((grad-tgt)/prate)*config.pumptime; % prate is given in terms of each pumptime cycle
    idealTime=max(config.pumptime,min(maxAccel*config.pumptime,idealTime)); % minimum pump time is config.pumptime
    pumpTime = [pumpTime,idealTime]; % Store
    tic;
    smset('PulseLine',pulseLine,[],'quiet') % Turn on pumping 
    mpause(idealTime-nullTime); % wait 
    smset('PulseLine',params(ind).offPulse,[],'quiet') % Turn off pumping
    totTime=toc; 
    totTime = totTime -nullOffTime;
    b=totTime/config.pumptime; % Number of "cycles" pumped for. 
    F(1,2:3)=F(1,2:3)*b; % Correct process matrix for actual pump time.
    
    x_priori = F * x ; % Predict the new state
    P_priori = F * P * (F') + Q; % Predict new process matrix. 
    
    [grad,gradOpts] = getgradient(fbScan,gradOpts); % Measure the new state.    
    gradDev = gradOpts.gradDev;
    grad = abs(grad)*sign(x_priori(1));  % Guess sign of gradient.
    err(j)=grad-tgt;
    obsHist=[obsHist , [grad ; gradDev]]; % List grad and noise 
        
    y = grad - x_priori(1);   % Innovation
    S = P_priori(1,1) + gradDev^2; % Innovation covariance
    K = P_priori * (H') / S; % Optimal gain 
    x = x_priori + K * y; % a-posteriori value of x
    P = (eye(size(P)) - K * H) * P_priori; % New process matrix  
    xHist=[xHist, x]; % List of xvals (grads and rates)
    if isnan(x(1))
        fprintf('Got Lost \n'); 
        return
    elseif y > params(ind).MaxError
        fprintf('Error too big \n'); 
        return
    elseif (xHist(end)-xHist(end-1))>params(ind).MaxChange
        fprintf('Changed too much \n');
        return
    end
    % Check for sign error; if the apparent pump rate is negative, we probably have the wrong sign on the gradient.
    if flipcount >= 6 && ((x(2) < -0.5*sqrt(P(2,2))) || (x(3) < -0.5*sqrt(P(3,3))))  % We appear to have misidentified gradient
        x(1)=-x(1);
        x(2:3)=abs(x(2:3));
        flipcount = 0;
    else
        flipcount=flipcount+1;
    end
    if x(2) > 20 || x(3) > 20
        x(2:3)=[1 1];
        P(2,2)=1e4; P(3,3)=1e4;
    end
    if config.figure % Update the plots        
        xvals=1:size(xHist,2);
        pumpon=find(pumpHist > 0);
        if ~isfield(config,'filterhistory')
            subplot(2,2,3); cla;
            config.filterhistory=plot(xHist(1,:),'k'); hold on; % 
            config.pumppoints=plot(xvals(pumpon),xHist(1,pumpon),'r.');            
            xlabel('Iteration'); ylabel('Gradient');
            subplot(2,2,4); cla;
            config.uprate=plot(xHist(2,:),'g.-'); hold on;
            config.downrate=plot(xHist(3,:),'b.-');
            xlabel('Iteration'); ylabel('S(g)/T(b) rate MHz/cycle');
        else
            set(config.filterhistory,'XData',xvals,'YData',xHist(1,:));
            set(config.pumppoints,'XData',xvals(pumpon),'YData',xHist(1,pumpon));            
            set(config.uprate,'XData',xvals,'YData',xHist(2,:));
            set(config.downrate,'XData',xvals,'YData',xHist(3,:));
        end
        if get(config.figure, 'CurrentCharacter') == char(27)
            set(config.figure, 'CurrentCharacter', char(0));
            sleep;
            return
        end
    end
end
fbdata.params(ind).prate = x(2:3); % Store info learned in fbdata. 
fbdata.params(ind).pMatrix = P; 
clearMask(config.datachan); % when finishing, delete the mask. 
if j < config.attempts 
    fbdata.set(end+1) = j;
    fbdata.pumpHist{end+1} = pumpTime; 
    fbdata.gradHist{end+1} = xHist; 
    good = true; 
else
    fbdata.set(end+1) = false;
end
end