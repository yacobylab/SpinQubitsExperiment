function biasTimer(opts)
% Set up and run several different types of timers. 
% opts: init bias, init dmm, init lock, start, stop 
% init: create a new timer, clearing global variables. 
% by default, this also starts the timer. To just create it, include stat
% in opt. 
% start: start a timer that's already been initialized. 
% stop: stop timer and save data. 
%dmm: read and plot dmm channel every 100 s 
% lock: read and plot dmm channel every 10 s
% bias: read and plot resistance for bias cooling every 10 s, assuming 18
% gates are connected through 100K resistor to dmm. 
% Edit channel names as necessary.
global logdata
if ~exist('opts','var'), error('Please provide argument. See help'); end

if isopt(opts,'bias')
    logdata.timer=timer('BusyMode','queue','ExecutionMode','fixedRate','Name','biasCurr','Period',10,'TimerFcn',@biascurr);
    ylab='Resistance (MOhm)';
elseif isopt(opts,'lock')
    ylab='Current';
    logdata.timer=timer('BusyMode','queue','ExecutionMode','fixedRate','Name','lockin','Period',10,'TimerFcn',@lockRead);
elseif isopt(opts,'dmm')
    ylab='Resistance';
    logdata.timer=timer('BusyMode','queue','ExecutionMode','fixedRate','Name','dmm','Period',100,'TimerFcn',@dmmRead);
end

if isopt(opts,'init')
    clear global rbias; clear global tm;
    logdata.timer.ErrorFcn = @errfunction; % Try restarting the timer with a 5 min start delay
    logdata.timer.StopFcn = @stopperfn; %print out a stopping msg
    global rbias; global tm;
    tm = []; rbias = [];
    figure(5); clf; hold on; xlabel('Time'); ylabel(ylab);
    if ~isopt(opts,'stat'), biasTimer('start'); end
end

if isopt(opts,'start'), start(logdata.timer); end
if isopt(opts,'stop')
    global tm; global rbias;
    stop(logdata.timer);
    save(sprintf('Z:/Shannon/data/biasCool-%s',datestr(now,1)),'tm','rbias');
end
end

function biascurr(~,~)
% use a 100k resistor. For now, assume 10 gates. 
r = 1e5; ngates = 18; 
global rbias; global tm;
tm(end+1) = now;
vres=cell2mat(smget('DMM'));
vin = diff(cell2mat(smget({'SD1mid','SD1top'})));
vrat = vres / vin; 
rbias(end+1) =ngates *  r *  (1 - vrat) / vrat; 

figure(5);
plot(datetime(datestr(tm)),rbias/1e6,'b.');
end

function lockRead(~,~)
global rbias; global tm;
tm(end+1) = now;
rbias(end+1)=cell2mat(smget('LockinA'));
figure(5);
plot(datetime(datestr(tm)),rbias,'b.');
end

function dmmRead(~,~)
global rbias; global tm;
tm(end+1) = now;
rbias(end+1)=cell2mat(smget('DMM'));
figure(5);
plot(datetime(datestr(tm)),rbias,'b.');
end

function errfunction(~, ~)
% try to restart the timer with a 50 sec start delay
global logdata;
stop(logdata.timer);
t = get(logdata.timer);
t.StartDelay = 50;
start(t);
end

function stopperfn(timerhandle, ~)
Tinfo = get(timerhandle);
fprintf('Stopping Timer %s \n', Tinfo.Name);
end