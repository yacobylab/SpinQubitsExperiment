function [xHist,obsHist]=characPump(fbScan,config)
% Check how feedback is working, measuring pump rate in both directions. 
% function [xHist,obsHist]=characPump(fbScan,config)
%  options:
%      grpname   Group name or number for dBz measurement
%                   default: first group that starts with dBz_
%     pumptime   Time to pump. Default 0.005
%     attempts   How many times to try. Default 50
%       figure   Figure to use for status display.  Default 1035
% The heart of this code is a Kalman filter.  See wikipedia page for similar notation.
%   opts: tl, do tl pumping (default is stp). 

global fbdata; global tuneData;
params = fbdata.params; 
if ~exist('config','var')
    config = struct;
elseif ischar(config)
    config.opts = config;
elseif iscell(config)
    config=struct(config{:});
end
awgcntrl('on start wait err');

config=def(config,'datachan',tuneData.dataChan);
config=def(config,'grpname','');
config=def(config,'pumptime',0.005); 
config=def(config,'attempts',50);
config=def(config,'opts','');
config=def(config,'gopts','nopol nodisp');
config=def(config,'figure',1035);
ind = str2double(config.datachan(end)); % side
config = def(config,'nloop',2*fbdata.params(ind).nloopfb); 
gradOpts=config; gradOpts.opts=gradOpts.gopts; gradOpts = rmfield(gradOpts,'gopts'); 
flipcount=0; xHist=[]; obsHist=[]; pumpHist=[];
if ~exist('fbScan','var') || isempty(fbScan) % Make the scan
    fbGroup = fbdata.params(ind).fbInit; 
    fbdata.params(ind).fitType = 'fit'; 
    fbScan = makeFeedbackScan(fbGroup,config.nloop,tuneData.dataChan);     
end
[grad,gradOpts] = getgradient(fbScan,gradOpts); % measure gradient.
gradOpts.opts = [gradOpts.opts 'reget']; % After first time, can just rerun prefns. 
x = [grad; params(ind).prate]; % state of filter: gradient, singlet rate, triplet rate. 

if isfield(params,'pMatrix') && ~isempty(params(ind).pMatrix)
    P = fbdata.params(ind).pMatrix; 
else
    % Estimated covariance matrix.  We start with no knowledge of pump rates.
    P = diag([gradOpts.gradDev^2, 1,1]); 
end
% The gradient fluctuates a lot, but the pump rate should be pretty stable.
Q=diag([5 1 1]); % Process noise, half-assed guess.
H = [1 0 0]; % Perform a correction step.
nullTime = 2.5e-3; nullOffTime = 2.5e-3; % Approx time for computer to turn on or off.
if isopt(config.opts,'long'), config.attempts=500; end
for j = 2:config.attempts    
    if ~isopt(config.opts,'tl') % Pump singlet direction
        pulseLine=params(ind).singletPulse;
        F=[1 1 0; 0 1 0 ; 0 0 1];
        pumpHist=[pumpHist 1]; %#ok<*AGROW>
        pumpInd = 2; 
    else % Pump triplet direction. 
        pulseLine = params(ind).tripletPulse;
        F=[ 1 0 -1 ; 0 1 0 ; 0 0 1];
        pumpHist=[pumpHist 1];
        pumpInd = 3;
    end    
    tic;
    smset('PulseLine',pulseLine,[],'quiet') % Turn on pumping
    mpause(config.pumptime-nullTime); % Wait
    smset('PulseLine',params(ind).offPulse,[],'quiet') % Turn off pumping
    totTime = toc; % Check correct time 
    totTime = totTime -nullOffTime;
    [grad,gradOpts] = getgradient(fbScan,gradOpts); % Measure the new state.  
    
    b=totTime/config.pumptime;
    F(1,2:3)=F(1,2:3)*b; % Correct process matrix for actual pump time.
    
    x_priori = F * x ; % Predict the new state
    P_priori = F * P * (F') + Q; % Predict new process matrix. 
        
    gradDev = gradOpts.gradDev;
    grad = abs(grad)*sign(x_priori(1));  % Guess sign of gradient.
    obsHist=[obsHist , [grad ; gradDev]]; % List grad and noise 
        
    y = grad - x_priori(1);   % Innovation
    S = P_priori(1,1) + gradDev^2; % Innovation covariance
    K = P_priori * (H') / S; % Optimal gain 
    x = x_priori + K * y; % a-posteriori value
    P = (eye(size(P)) - K * H) * P_priori; % New process matrix  
    xHist=[xHist x]; % List of xvals (grads and rates)
    
    % Check for sign error; if the apparent pump rate is negative, we probably have the wrong sign on the gradient.
    if flipcount >= 3 && (x(2) < -0.5*sqrt(P(2,2)) || x(3) < -0.5*sqrt(P(3,3)))  % We appear to have misidentified gradient
        x(1)=-x(1);
        x(2:3)=abs(x(2:3));
        flipcount = 0;
    else
        flipcount=flipcount+1;
    end
    if x(2) > 20 || x(3) > 20 % Huge pump rates
        x(2:3)=[1 1];
        P(2,2)=1e4; P(3,3)=1e4;
    end
    if config.figure % Update the plots        
        xvals=1:size(xHist,2);
        if ~isfield(config,'filterhistory')
            subplot(2,2,3); cla;
            config.filterhistory=plot(xHist(1,:)); hold on; % 
            xlabel('Iteration'); ylabel('Gradient (MHz)');
            subplot(2,2,4); cla;
            config.uprate=plot(xHist(pumpInd,:),'.-'); hold on;
            xlabel('Iteration'); ylabel('Pump rate MHz/cycle');
        else
            set(config.filterhistory,'XData',xvals,'YData',xHist(1,:));
            set(config.uprate,'XData',xvals,'YData',xHist(2,:));
        end
        drawnow; 
        if get(config.figure, 'CurrentCharacter') == char(27)
            set(config.figure, 'CurrentCharacter', char(0));
            sleep;
            clearMask(config.datachan);
            fbdata.params(ind).fitType = 'fft';
            return
        end
    end
end
clearMask(config.datachan);
fbdata.params(ind).fitType = 'fft'; 
fbdata.params(ind).pMatrix=P; 
fbdata.params(ind).prate = x(2:3); 
end