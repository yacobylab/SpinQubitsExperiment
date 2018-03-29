function biasTimer(cntl)
% cntl: initBias, initDMM, start, stop 
global logdata

switch cntl
    case 'initBias'
        clear global rbias; 
        clear global tm; 
      
        logdata.timer=timer('BusyMode','queue','ExecutionMode','fixedRate','Name','biascurr','Period',10,'TimerFcn',@biascurr);
        logdata.timer.ErrorFcn = @errfunction; %try restarting the timer with a 5 min start delay
        logdata.timer.StopFcn = @stopperfn; %print out a stopping msg
        global rbias; global tm; 
        tm = []; rbias = [];
        figure(5); clf; hold on; xlabel('Time'); ylabel('Resistance (MOhm)'); 
    case 'initDMM'
        global rbias; global tm; 
        tm = []; rbias = [];
        figure(5); clf; hold on; xlabel('Time'); ylabel('Resistance'); 
        logdata.timer=timer('BusyMode','queue','ExecutionMode','fixedRate','Name','biascurr','Period',100,'TimerFcn',@dmmRead);
        logdata.timer.ErrorFcn = @errfunction; %try restarting the timer with a 5 min start delay
        logdata.timer.StopFcn = @stopperfn; %print out a stopping msg
    case 'initLock'
        global rbias; global tm; 
        tm = []; rbias = [];
        figure(5); clf; hold on; xlabel('Time'); ylabel('Current'); 
        
        logdata.timer=timer('BusyMode','queue','ExecutionMode','fixedRate','Name','biascurr','Period',10,'TimerFcn',@lockRead);
        logdata.timer.ErrorFcn = @errfunction; %try restarting the timer with a 5 min start delay
        logdata.timer.StopFcn = @stopperfn; %print out a stopping msg
    case 'start'
        start(logdata.timer);
    case 'stop'
        global tm 
        global rbias
        stop(logdata.timer);        
        save(sprintf('Z:\\Shannon\\data\\biasCool-%s',datestr(now,1)),'tm','rbias');
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
% use a 100k resistor. For now, assume 10 gates. 
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