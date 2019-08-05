function [good,pumpTime,pumpHist]=setgradientMulti(~,fbScan,tgt,config)
% Lock the gradient to tgt (default:setpt given in fbdata.params.setpt). Perform on multiple qubits. 
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

config=def(config,'grpname','');
config=def(config,'attempts',50);
config=def(config,'opts','both');
if isopt(config.opts,'both'), config.datachan = 'both'; end
if isopt(config.opts,'long'), config.attempts=500; end
config=def(config,'gopts','nopol nodisp');
config=def(config,'figure',1035);

if isopt(config.opts,'both')
    config = def(config,'datachan','both');
    ind = 3; indL = [1,2]; 
else
    config=def(config,'datachan',tuneData.dataChan);    
    ind = str2double(tuneData.dataChan(end));
    indL = ind; 
end
nQub = length(indL); 
config=def(config,'tol',fbdata.params(ind).tolInit);
config=def(config,'nloop',fbdata.params(ind).nloopfb); 
config=def(config,'pumptime',fbdata.params(ind).pumptime);
if ~exist('tgt','var') || isempty(tgt)
    for i = 1:nQub
        tgt(i)=fbdata.params(indL(i)).setpt;
    end
end
% gradOpts configures getting gradient. 
gradOpts=config; gradOpts.opts=gradOpts.gopts; gradOpts = rmfield(gradOpts,'gopts');
xHist=[]; obsHist=[]; pumpHist=[];
% Default scan for measuring gradient, usually used when calling setgradient from cmd line. 
if ~exist('fbScan','var') || isempty(fbScan)    
    if ~isopt(config.opts,'fine')
        fbGroup = fbdata.params(ind).fbInit; % Init group 1 ns spacing, less aliasing. 
    else
        fbGroup = fbdata.params(ind).fbGroup; 
    end
    fbScan = makeFeedbackScan(fbGroup,config.nloop,config.datachan);     
end
%% Measure gradient, initalize filter. 
[grad,gradOpts] = getgradientMulti(fbScan,gradOpts); % Measure gradient.
if any(isnan(grad)), [grad,gradOpts] = getgradientMulti(fbScan,gradOpts); end
if any(isnan(grad)), [grad,gradOpts] = getgradientMulti(fbScan,gradOpts); end

if any(isnan(grad))
    if ~isopt(config.opts,'two')
        pause(6);
        setGrad = setgradientMulti([],[],[],struct('opts',[config.opts ' two']));
        [grad,gradOpts] = getgradientMulti(fbScan,gradOpts); % Measure gradient.
        %if ~setGrad || isnan(grad)
        %    pause(6);
        %    fprintf('Couldn''t set at start \n');
        %    return
        %end
    %else
    %    return;
    %end
    end
end
gradOpts.opts = [gradOpts.opts 'reget']; % After first time, faster to not run whole scan.

for i = 1:nQub
    if isopt(config.opts, 'init') || any(isnan(fbdata.params(ind).prate))
        fbdata.params(indL(i)).prate = 250 * config.pumptime *[1;1]; % 0.25 Mhz / ms.
    end
    kalm(i).x=[grad(i); fbdata.params(indL(i)).prate]; % State of filter: gradient, singlet rate, triplet rate.
    
    % Estimated covariance matrix.  We start with no knowledge of pump rates.
    if isfield(fbdata.params,'pMatrix') && ~isempty(fbdata.params(indL(i)).pMatrix) && ~isopt(config.opts,'init')
        kalm(i).P = fbdata.params(indL(i)).pMatrix;
    else
        kalm(i).P=diag([gradOpts.gradDev(i)^2, 5,5]);
    end
    kalm(i).y = NaN; 
    % Initialize Kalman    
    j=1;
    err(i,j) = grad(i)-tgt(i);  % Initialize error vector.        
end
if any(err(i,j) > fbdata.params(ind).MaxError)
    if ~isopt(config.opts,'two')
        pause(6);
        fbdata.params(ind).MaxError = 150;
        setGrad = setgradientMulti([],[],[],struct('opts',[config.opts ' two']));
        fbdata.params(ind).MaxError = 50;
        [grad,gradOpts] = getgradientMulti(fbScan,gradOpts); % Measure gradient.
        if ~setGrad || isnan(grad)
            pause(6);
            fprintf('Couldn''t set at start \n');
            return
        end
    else
        return;
    end
end
% The gradient fluctuates a lot, but the pump rate should be pretty stable.

maxAccel=25; % Maximum number of 'pumptime' increments allowed; too large and will lose
% track of gradient, too small and things will take forever.
pumpTime = []; 

nullTime = 2.5e-3; nullOffTime = 2.5e-3; % Time it takes to call functions (in which we're not pumping)
% This was measured on computer, may want to check this periodically. 
%% Run Kalman filter
% Kalman process runs until tolerance satisfied or max attempts reached.
singInd = 2; tripInd = 3; 
while (any(abs(err(:,j)) > config.tol) || any(isnan(err(:,j)))) && j <= config.attempts      
    for i = 1:nQub % this is the set of qubits 
        % Configure type of pumping depending on gradient value, desired gradient side.
        kalm(i).F = eye(3);
        if abs(err(i,j)) < config.tol
            pump{i} = 'null';
            pumpDir = 0; pumpInd=[];
        else
            if strcmp(fbdata.params(indL(i)).dir,'sing') && grad(i) < tgt(i) % Tries to pump gradient on singlet side. To change, change this, change F, pumpHist.
                pump{i} = 'sing'; pumpInd = singInd;            
                pumpDir = 1;
            elseif strcmp(fbdata.params(indL(i)).dir,'sing') && grad(i) > tgt(i)                
                pump{i} = 'trip'; pumpInd = tripInd;
                pumpDir = -1;
            elseif strcmp(fbdata.params(indL(i)).dir,'trip') && grad(i) < tgt(i) % Pump triplet.                
                pump{i} = 'trip'; pumpInd = tripInd;
                pumpDir = 1;
            else                
                pump{i} = 'sing'; pumpInd = singInd;
                pumpDir = -1;
            end
        end
        kalm(i).F(1,pumpInd) = pumpDir;
        
        pumpHist(i,j)=pumpDir; 
        % F represents change matrix. So we expect pump rates to be the same, gradient to change re pump time.
        if ~isempty(pumpInd)
            prate=kalm(i).x(pumpInd); % pump rate
            idealTime(i)=abs((grad(i)-tgt(i))/prate)*config.pumptime; % prate is given in units of config.pumptime
            idealTime(i)=max(config.pumptime,min(maxAccel*config.pumptime,idealTime(i))); % minimum pump time is config.pumptime
        else
            idealTime(i) = Inf;
        end
    end
    % Now, go through and figure out which pulse to use. 
    if strcmp(pump{1},'sing')
        if strcmp(pump{2},'sing')
            pulseLine = fbdata.multi.SS; 
        elseif strcmp(pump{2},'trip')
            pulseLine = fbdata.multi.ST; 
        else
            pulseLine = fbdata.params(indL(1)).singletPulse; 
        end
    elseif strcmp(pump{1},'trip')
        if strcmp(pump{2},'sing')
            pulseLine = fbdata.multi.TS; 
        elseif strcmp(pump{2},'trip')
            pulseLine = fbdata.multi.TT; 
        else
            pulseLine = fbdata.params(indL(1)).tripletPulse; 
        end
    else % No need to pump first qubit.
        if strcmp(pump{2},'sing')
            pulseLine = fbdata.params(indL(2)).singletPulse; 
        elseif strcmp(pump{2},'trip')
            pulseLine = fbdata.params(indL(2)).tripletPulse; 
        else
           fprintf('No pumping needing. Why are you here? \n'); 
           break; 
        end
    end
    idealTimeAll = min(idealTime);
    tic;
    % Quiet speeds things up by not updating graphics.
    smset('PulseLine',pulseLine,[],'quiet') % Turn on pumping.
    mpause(idealTimeAll-nullTime); % wait
    smset('PulseLine',fbdata.params(ind).offPulse,[],'quiet') % Turn off pumping
    totTime=toc;
    totTime = totTime - nullOffTime;
    pumpTime(j) = totTime; % Store
    b=totTime/config.pumptime; % Number of "cycles" pumped for.    
    
    [grad,gradOpts] = getgradientMulti(fbScan,gradOpts); % Measure the new state.
    % Try to remeasure twice if it fails.
    if any(isnan(grad)), [grad,gradOpts] = getgradientMulti(fbScan,gradOpts); end
    if any(isnan(grad)), [grad,gradOpts] = getgradientMulti(fbScan,gradOpts); end
    if any(isnan(grad)), fprintf('Reading gradient nan \n'), break; end
    gradDev = gradOpts.gradDev;
    
    j = j+1;    
    err(:,j)=grad-tgt;    
    obsHist=[obsHist , [grad ; gradDev]]; % List grad and noise
    
    for i = 1:nQub
        kalm(i).F(1,2:3)=kalm(i).F(1,2:3)*b; % Correct process matrix for actual pump time.
        kalm(i) = updateKalman(kalm(i),grad(i),gradDev(i));
        xHist{i}(1:3,j-1)=kalm(i).x; % List of xvals (grads and rates)
        if isnan(kalm(i).x(1))
            fprintf('Got Lost \n');
            break
        elseif kalm(i).y > fbdata.params(indL(i)).MaxError
            fprintf('Error %2.1f MHz, too big \n',kalm(i).y);
            break
        elseif j>2 && (xHist{i}(end,1)-xHist{i}(end-1,1)) > fbdata.params(indL(i)).MaxChange
            fprintf('Changed too much \n');
            break
        end
        % Check for sign error; if the apparent pump rate is negative, we probably have the wrong sign on the gradient.
        % this may need to be improved.  % Not sure why error has been so big lately, there may be an issue.
        if flipcount >= 6 && ((kalm(i).x(2) < -0.04*sqrt(kalm(i).P(2,2))) || (kalm(i).x(3) < -0.04*sqrt(kalm(i).P(3,3)))) % We appear to have misidentified gradient
            kalm(i).x(1)=-kalm(i).x(1);
            kalm(i).x(2:3)=abs(kalm(i).x(2:3));
            flipcount = 0;
        else
            flipcount=flipcount+1;
        end
        if kalm(i).x(2) > 20 || kalm(i).x(3) > 20
            kalm(i).x(2:3)=[1 1];
            kalm(i).P(2,2)=1e4; kalm(i).P(3,3)=1e4;
        end
        if config.figure % Update the plots
            xvals=1:size(xHist{i},1);
            pump1=find(pumpHist(i,:) > 0);            
            pump2=find(pumpHist(i,:) < 0);            
            if ~isfield(config,'filterhistory') % Create plot first time
                a1=subplot(2,2,3); cla; 
                xlabel('Iteration'); ylabel('Gradient (MHz)');
                 legend(a1,'show')
                a2 = subplot(2,2,4); cla;
                xlabel('Iteration'); ylabel('Pump rate (MHz/cycle)');
                legend(a2,'show')
            end
            if j==2 % First round
                config.filterhistory(i)=plot(a1,xHist{i}(1,:),'DisplayName',num2str(i)); 
                hold on; 
                %config.pumppoints(i)=plot(a1,xvals(pumpOn),xHist{i}(pumpOn,1),'r.');
                config.uprate(i)=plot(a2,xHist{i}(2,:),'.-','DisplayName',sprintf('Sing %d',i),'Color',a1.ColorOrder(i,:)); 
                hold on;
                config.downrate(i)=plot(a2,xHist{i}(3,:),'.-','DisplayName',sprintf('Trip %d',i),'Color',a1.ColorOrder(2+i,:));                
            else % After first time, update plot
                set(config.filterhistory(i),'YData',xHist{i}(1,:));
                %set(config.pumppoints(i),'XData',xvals(pumpOn),'YData',xHist{i}(pumpOn,1));
                set(config.uprate(i),'XData',pump1,'YData',xHist{i}(2,pump1));
                set(config.downrate(i),'XData',pump2,'YData',xHist{i}(3,pump2));
                drawnow
            end
             % Exit program if esc key pressed on figure.
            if get(config.figure, 'CurrentCharacter') == char(27)
                set(config.figure, 'CurrentCharacter', char(0));
                sleep('fast');
                return
            end
        end
    end           
end

clearMask('DAQ1'); % When done, delete the mask.
if all(abs(err(:,end)) < config.tol)
    fbdata.set(end+1) = j;
    good = true;
    for i = 1:nQub
        fbdata.params(indL(i)).prate = kalm(i).x(2:3); % Store info learned in fbdata.
        fbdata.params(indL(i)).pMatrix = kalm(i).P;
    end
else
    fbdata.set(end+1) = NaN;
end
fbdata.pumpHist{end+1} = pumpTime;
fbdata.gradHist{end+1} = xHist;
fbdata.err{end+1}= err;   
end