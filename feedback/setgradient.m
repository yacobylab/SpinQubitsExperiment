function [good,pumpTime,pumpHist]=setgradient(~,fbScan,tgt,config)
% Lock the gradient to tgt (default:setpt given in fbdata.params.setpt)
% function [good,pumpTime,pumpHist=setgradient(~,fbScan,tgt,config)
% Code used in two ways: from command line, to check if you can set
% gradient, and from within scan, to keep gradient locked during other
% measurements. First case can be called directly, for scanning, use
% swfbConfig. 
% config:
%     grpname   Group name or number for dBz measurement
%                   default: first group that starts with dBz_
%     pumptime   Minimum time to pump when trying to lock.  defined in fbdata.params.
%     attempts   How many times to try before giving up.  Default 50
%     figure   Figure to use for status display.  Default 1035
%     datachan: which channel to take data  on. 
%     tol: Setpoint tolerance. Default from fbdata 
%     dir: What side to stabilize gradient on (sing/trip). default from fbdata 
%     gopts: nopol, nodisplay
%     opts: 
%          long: change number of attempts to 500 instead of 50.
%          fine: Use fbGroup with wider spacing to set gradient (more
%          precise). 
%          init: initialize new covariance matrix
% The heart of this code is a Kalman filter.  See wikipedia page for similar notation.
% Null first argument due to sm prefn convention.

%% Configure scan, set up Kalman filter
global fbdata; global tuneData; 
flipcount=0; good = 0;
if ~exist('config','var'), config=struct(); end
config=def(config,'datachan',tuneData.dataChan);
config=def(config,'grpname','');
config=def(config,'attempts',50);
config=def(config,'opts','');
if isopt(config.opts,'long'), config.attempts=500; end
config=def(config,'gopts','nopol nodisp');
config=def(config,'figure',1036);
ind = str2double(config.datachan(end)); % determine side
config=def(config,'tol',fbdata.params(ind).tolInit);
config=def(config,'dir',fbdata.params(ind).dir);
config=def(config,'pumptime',fbdata.params(ind).pumptime);
config=def(config,'nloop',fbdata.params(ind).nloopfb); 
if ~exist('tgt','var') || isempty(tgt), tgt=fbdata.params(ind).setpt; end

% gradOpts configures getting gradient. 
gradOpts=config; gradOpts.opts=gradOpts.gopts; gradOpts = rmfield(gradOpts,'gopts');
xHist=[]; obsHist=[]; pumpHist=[];
% Default scan, usually used when calling setgradient from cmd line. 
if ~exist('fbScan','var') || isempty(fbScan)    
    if ~isopt(config.opts,'fine')
        fbGroup = fbdata.params(ind).fbInit; % Init group 1 ns spacing, less aliasing. 
    else
        fbGroup = fbdata.params(ind).fbGroup; 
    end
    fbScan = makeFeedbackScan(fbGroup,config.nloop,tuneData.dataChan);     
end
%% Measure gradient, initalize filter. 
[grad,gradOpts] = getgradient(fbScan,gradOpts); % Measure gradient.
gradOpts.opts = [gradOpts.opts 'reget']; % After first time, faster to not run whole scan.
if isnan(grad), [grad,gradOpts] = getgradient(fbScan,gradOpts); end
if isnan(grad), [grad,gradOpts] = getgradient(fbScan,gradOpts); end
pauseTime = 0.5; 
if isnan(grad)
    fprintf('Gradient nan, trying to reset \n'); 
    if ~isopt(config.opts,'two')
        pause(pauseTime);
        setGrad = setgradient([],[],[],struct('opts','two'));
        [grad,gradOpts] = getgradient(fbScan,gradOpts); % Measure gradient.
        if ~setGrad || isnan(grad)
            pause(pauseTime);
            fprintf('Couldn''t set at start \n');
            return
        end
    else
        return;
    end    
end

if isopt(config.opts, 'init') || any(isnan(fbdata.params(ind).prate))
    fbdata.params(ind).prate = 250 *config.pumptime *[1;1]; % 0.25 Mhz / ms. 
end
x=[grad; fbdata.params(ind).prate]; % State of filter: gradient, singlet rate, triplet rate.

% Estimated covariance matrix.  We start with no knowledge of pump rates.
if isfield(fbdata.params,'pMatrix') && ~isempty(fbdata.params(ind).pMatrix) && ~isopt(config.opts,'init')
    P = fbdata.params(ind).pMatrix;
else
    P=diag([gradOpts.gradDev^2, 5,5]);
end

% Initialize Kalman
% The gradient fluctuates a lot, but the pump rate should be pretty stable.
Q=diag([5 1 1]); % Process noise, half-assed guess.
j=1;
err(j)=grad-tgt; % Initialize error vector. 
if abs(err) > fbdata.params(ind).MaxError
    if ~isopt(config.opts,'two')
        pause(pauseTime);
        fbdata.params(ind).MaxError = 150; 
        setGrad = setgradient([],[],[],struct('opts','two'));
        fbdata.params(ind).MaxError = 50; 
        [grad,gradOpts] = getgradient(fbScan,gradOpts); % Measure gradient.
        if ~setGrad || isnan(grad)
            pause(pauseTime);
            fprintf('Couldn''t set at start \n');
            return
        end
    else
        return;
    end
end

H = [1 0 0]; % Perform a correction step.
maxAccel=25; % Maximum number of 'pumptime' increments allowed; too large and will lose
% track of gradient, too small and things will take forever.
pumpTime = []; 

nullTime = 2.5e-3; nullOffTime = 2.5e-3; % Time it takes to call functions (in which we're not pumping)
% This was measured on computer, may want to check this periodically. 
%% Run Kalman filter
% Kalman process runs until tolerance satisfied or max attempts reached.
while (abs(err(j)) > config.tol || isnan(err(j))) && j <= config.attempts
    j = j+1;
    % Configure type of pumping depending on gradient value, desired gradient side.
    if strcmp(config.dir,'sing') && grad < tgt % Tries to pump gradient on singlet side. To change, change this, change F, pumpHist.
        pulseLine=fbdata.params(ind).singletPulse;
        F=[1 1 0; 0 1 0 ; 0 0 1]; % Matrix for change in x.
        pumpDir = 1; 
        pumpHist=[pumpHist pumpDir]; %#ok<*AGROW>
        prate=x(2); % pump rate
    elseif strcmp(config.dir,'sing') && grad > tgt
        pulseLine = fbdata.params(ind).tripletPulse;
        F=[ 1 0 -1 ; 0 1 0 ; 0 0 1];
        pumpDir = 1;
        pumpHist=[pumpHist pumpDir];
        prate=x(3);
    elseif strcmp(config.dir,'trip') && grad < tgt % Pump triplet.
        pulseLine=fbdata.params(ind).tripletPulse;
        F=[1 0 1; 0 1 0 ; 0 0 1]; % Matrix for change in x.
        pumpDir = 1; 
        pumpHist=[pumpHist pumpDir]; %#ok<*AGROW>        
        prate=x(3);
    else
        pulseLine = fbdata.params(ind).singletPulse;
        F=[ 1 -1 0 ; 0 1 0 ; 0 0 1];
        pumpDir = -1; 
        pumpHist=[pumpHist pumpDir];
        prate=x(2);
    end
    % F represents change matrix. So we expect pump rates to be the same, gradient to change re pump time.
    idealTime=abs((grad-tgt)/prate)*config.pumptime; % prate is given in units of config.pumptime
    idealTime=max(config.pumptime,min(maxAccel*config.pumptime,idealTime)); % minimum pump time is config.pumptime
    
    tic;
    % Quiet speeds things up by not updating graphics.
    smset('PulseLine',pulseLine,[],'quiet') % Turn on pumping.
    mpause(idealTime-nullTime); % wait
    smset('PulseLine',fbdata.params(ind).offPulse,[],'quiet') % Turn off pumping
    totTime=toc;
    totTime = totTime - nullOffTime;
    pumpTime = [pumpTime,totTime]; % Store
    b=totTime/config.pumptime; % Number of "cycles" pumped for.
    F(1,2:3)=F(1,2:3)*b; % Correct process matrix for actual pump time.
    
    x_priori = F * x ; % Predict the new state
    P_priori = F * P * (F') + Q; % Predict new process matrix.
    
    [grad,gradOpts] = getgradient(fbScan,gradOpts); % Measure the new state.
    % Try to remeasure twice if it fails. 
    if isnan(grad), [grad,gradOpts] = getgradient(fbScan,gradOpts); end
    if isnan(grad), [grad,gradOpts] = getgradient(fbScan,gradOpts); end
    if isnan(grad), fprintf('Reading gradient nan \n'), break; end    
    gradDev = gradOpts.gradDev;
    grad = abs(grad)*sign(x_priori(1));  % Guess sign of gradient.
    
    err(j)=grad-tgt;    
    obsHist=[obsHist , [grad ; gradDev]]; % List grad and noise
    
    y = grad - x_priori(1); % Innovation
    S = P_priori(1,1) + gradDev^2; % Innovation covariance
    K = P_priori * (H') / S; % Optimal gain
    x = x_priori + K * y; % a-posteriori value of x
    P = (eye(size(P)) - K * H) * P_priori; % New process matrix
    xHist=[xHist, x]; % List of xvals (grads and rates)
    if isnan(x(1))
        fprintf('Got Lost \n');
        break
    elseif y > fbdata.params(ind).MaxError
        fprintf('Error %2.1f MHz, too big \n',y);
        break
    elseif (xHist(end)-xHist(end-1)) > fbdata.params(ind).MaxChange
        fprintf('Changed too much \n');
        break
    end
    % Check for sign error; if the apparent pump rate is negative, we probably have the wrong sign on the gradient.
    if flipcount >= 6 && ((x(2) < -0.05*sqrt(P(2,2))) || (x(3) < -0.05*sqrt(P(3,3)))) % We appear to have misidentified gradient
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
            drawnow
        end
        if get(config.figure, 'CurrentCharacter') == char(27)
            set(config.figure, 'CurrentCharacter', char(0));
            sleep('fast');
            return
        end
    end
end

clearMask(config.datachan); % When done, delete the mask.
if abs(err(end)) < config.tol
    fbdata.set(end+1) = j;
    good = true;
    
    fbdata.params(ind).prate = x(2:3); % Store info learned in fbdata.
    fbdata.params(ind).pMatrix = P;
else
    fbdata.set(end+1) = NaN;
end
fbdata.pumpHist{end+1} = pumpTime;
fbdata.gradHist{end+1} = xHist;
fbdata.err{end+1}= err;   
end