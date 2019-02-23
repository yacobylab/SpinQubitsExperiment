function measureQPCauto(dots,opts,div, badGates) 
% function measureQPCauto(dots,opts,div badGates) 
% measures all QPC gates on two qubit system. 
% Assumes that you are using LockinA, with all labels consistent. If not, you will have to dig through code.
% function measureQPCauto(dots,opts) 
% dot names are : SL, SR, QL, QR
% if only some dots bonded up, give cell or str of list to test for arg dots. Default is
% to do all. 
% opts : 
%   nobond: don't need to turn off other set of ohmics. (assumes ohmics
%   across dot are pinned together). 
%   quiet: don't tell you to change cables 
%   restart: refind closeVal. 
%   res: device with resonator; 1a and 3b should be pinned together on 4k. 
% div: resistance values in divider (for finding QOC) (default is 47e3 and
% 47 ohms)
% badGates: cell of gate names that you know don't work (won't test those
% ones). Note that problems will crop up 
% if you are running this on the fridge, will save the QOC values to
% scandata to set up hyst scan. 

global smdata; 
if ~exist('badGates','var'), badGates = {}; end
if ~exist('opts','var'), opts = ''; end
persistent firstTime closeVal % remember the closeVal so that you don't have to test each time you run. % fix me 
if ~exist('firstTime','var') || isopt(opts,'restart') % Retest for closeval 
    firstTime = 1; closeVal =[];
end
fulldotList = {'SL','SR','QL','QR'}; 
if ~exist('dots','var') || isempty(dots) 
    dots = fulldotList; 
elseif ischar(dots)
    dots = {dots};         
end

if isopt(opts,'nobond') % Do we need to cut off the other set of ohmics?
    zeroOL = 0;
elseif (~any(strcmpi(dots,'QL')) && any(strcmpi(dots,'QR'))) || ...
        (~any(strcmpi(dots,'QR')) && any(strcmpi(dots,'QL'))) && ~isopt(opts,'quiet')
    % If only testing one dot, check that the other one is bonded up. 
    action = questdlg('Is the other set of ohmics bonded up?','','Yes','No','No');
    if strcmp(action,'No')
        zeroOL = 0; 
    else 
        zeroOL = 1; 
    end
else
    zeroOL = 1; 
end
    
if ~exist('res','var') || isempty(res)
    divA = 47; divB = 47e3; 
else 
    div = sort(div); divA = div(1); divB = div(2); 
end

if isfield('smdata','folder') && isfield('smdata','qpcFolder') 
    cd([smdata.folder, smdata.qpcFolder]); 
elseif ~strcmp('Z:\qDots\data\data_2015_11_05\qpc_2016_05_16',pwd)
    cd Z:\qDots\data\data_2015_11_05\qpc_2016_05_16
end
is4k = contains(smdata.name,'4k'); % Use this to configure lock-in due to increased noise at 4k station. 

vOut = cell2mat(smget('LockExcA')); 
vOut = vOut * divA/(divA + divB); 
r0 = 80; % expected resistance per square of 2DEG 
gqoc = 7.748e-5; % one quantum of conductance. 
%% Set up scan
channels = {'N12','3b','2a','T34','4b','SD4mid','SD4top','SD4bot','4a','3a','N34','1b','2b','T12','SD1top','SD1mid','SD1bot'};%put back in 1a
ohmicsList = {'OS1top, OS1bot','OS4top, OS4bot', 'OL13, OL24','OL13, OL24'}; 
smset(channels,0);
logstr = '';

if is4k 
    scan.loops(1).ramptime = .2; 
    channels(end+1:end+2) = {'RF3','RF4'}; 
    scan.loops.rng = [0 -1.9];
else
    scan.loops(1).ramptime = .06; 
    smset('LockinTauA',0.03); 
    scan.loops(1).rng = [0 -2]; 
end

scan.loops(1).npoints= 100; scan.saveloop= [2 1]; 
lockin = 'LockinA'; 
scan.disp(1).channel = 1; scan.disp(1).dim = 1; scan.disp(1).loop=1;  
scan.loops(1).getchan=lockin;
scan.loops(1).testfn.fn = @atZero; scan.loops(1).testfn.args = {1e-10,10};
% to create a scan to use without this program, copy above 5 lines, set ramptime, npoints, setchan. 
scan.disp(2).channel = 2; scan.disp(2).dim = 1; scan.disp(2).loop=1;  
scan.loops(1).procfn(2).fn.fn = @(i,v,r0,gq) log(abs(1./(v./i-r0)/gq));
scan.loops(1).procfn(2).fn.inchan = 1; scan.loops(1).procfn(2).fn.outchan = 2;
scan.loops(1).procfn(2).dim = 1;

for j = 1:length(dots)
    name = dots{j};    
    try
        dotInd(j) = find(strcmp(fulldotList,name));     
    catch 
        error('Not a correct Dot Name!') 
    end % make sure that you do need to change the ohmics... 
    if isopt(opts,'quiet') || (j > 1 && (dotInd(j-1) == find(strcmp(fulldotList,'QL')) || dotInd(j-1) == find(strcmp(fulldotList,'QR'))) && (dotInd(j) == find(strcmp(fulldotList,'QL')) || dotInd(j) == find(strcmp(fulldotList,'QR'))))
    else
        action = questdlg(sprintf('Please change the cables to %s',ohmicsList{dotInd(j)}),'','Done','Quit','Done');
        if strcmp(action,'Quit')
            break
        end
    end
    switch name
        case 'SL' %Sensor L,
            title_str = 'Sensor_L';
            gate2={'SD1top','SD1mid','SD1bot'}; gate1 = {'1b';'1b';'1b'};
            fn = @(a) find(strcmp(gate2,a)); badGateNum=cell2mat(cellfun(fn,badGates,'UniformOutput',false));
            gate2(badGateNum)=[]; gate1(badGateNum)=[]; 
        case 'QL' %Qubit L, 
            title_str = 'Qubit_L';             
            if is4k
                gate2={'N12' '1a' '2a' '1b','2b','RF4','RF3'};
            elseif is4k && isopt(opts,'res')
                gate2 = {'N12','3b','2a','1b','2b','RF4','RF3'};
            else
                gate2={'N12' '1a' '2a' '1b','2b'};            
            end                        
            gate1val={'T12'}; gate1 = repmat(gate1val,5,1); gate1{end+1} = '2b'; gate1{end+1} = '1b';   
            fn = @(a) find(strcmp(gate2,a)); badGateNum=cell2mat(cellfun(fn,badGates,'UniformOutput',false));
            gate2(badGateNum)=[]; gate1(badGateNum)=[]; 
            if firstTime && zeroOL
                chans = {'T34','4b','3b','1b','2b','T12'};
                zeroval = 1e-10;
                closeVal = closeParallel(chans, lockin, zeroval);
                firstTime = 0;
            end
            closeChans ={'4b','3b','T34'};
            smset(channels,0); 
            if zeroOL              
            smset(closeChans,closeVal);
            end
        case 'SR' %Sensor R
            title_str = 'Sensor_R';
            gate2={'SD4top' 'SD4mid','SD4bot'}; gate1 = {'4b';'4a';'4a'};
            fn = @(a) find(strcmp(gate2,a)); badGateNum=cell2mat(cellfun(fn,badGates,'UniformOutput',false));
            gate2(badGateNum)=[]; gate1(badGateNum)=[]; 
        case 'QR' %Qubit R 
            title_str = 'Qubit_R';             
            if is4k
                gate2={'3a','3b','4a','4b','N34','RF4','RF3'};                    
            else
                gate2={'N34' '3a' '4a' '3b' '4b'};            
            end
            gate1val={'T34'}; gate1 = repmat(gate1val, 5,1);  gate1{end+1} = '4b'; gate1{end+1} = '3b';
            fn = @(a) find(strcmp(gate2,a)); badGateNum=cell2mat(cellfun(fn,badGates,'UniformOutput',false));
            gate2(badGateNum)=[]; gate1(badGateNum)=[];
            if firstTime && zeroOL
                chans = {'T34','4b','3b','1b','2b','T12'};
                zeroval = 1e-10;
                closeVal = closeParallel(chans, lockin, zeroval);
                firstTime = 0;
            end
            if zeroOL
                closeChans = {'1b','2b','T12'};
                smset(closeChans,closeVal);
            end
    end
    
    logstr = [logstr, sprintf('%s\n',title_str)];
    i = cell2mat(smget(lockin)); rfilt = abs(vOut) / i - r0;
    scan.loops(1).procfn(2).fn.args = {vOut,rfilt,gqoc};
    qocVal=[];
    for i=1:length(gate2)
        fprintf('Starting scan %i of %i, Gate %s\n',i,length(gate2),gate2{i});
        scan.loops(1).setchan={gate2{i},gate1{i,:}};
        if is4k
            scanName=smnext(sprintf('qpc_4k_%s_%s',gate1{i,1},gate2{i}));
        else
            scanName=smnext(sprintf('qpc_%s_%s',gate1{i,1},gate2{i}));
        end
        d=smrun(scan,scanName);
        dStore{i,1} = d{1}; dStore{i,2}=d{2};  %#ok<*SAGROW>
        f = gcf; 
        if f.CurrentCharacter == char(27)
            action = lower(questdlg('Go to next scan?','Question','yes','no','yes'));
            switch action
                case 'yes'
                    fprintf('Next scan...\n');
                otherwise                
                    fprintf('Exiting loop. Not setting gates to 0\n');
                    break
            end
        end
        if any(d{2}<0)
            [~,mi]= min(abs(d{2}));
            gvals = linspace(scan.loops(1).rng(1),scan.loops(1).rng(2),scan.loops(1).npoints);
            fprintf('1e^2/h at gate voltage %.3fV\n', gvals(mi));
            logstr = [logstr,sprintf('%s and %s: 1e^2/h at gate voltage %.3fV\n', gate2{i},gate1{i,1},gvals(mi))];
            qocVal(i)=gvals(mi);
            figure(1000); hold on;
            plot(scan.loops(1).rng,sign(nanmean(d{1}))*0*[1 1],'r');            
        else
            fprintf('Never reached 1 conductance quantum\n')
            logstr = [logstr,sprintf('Gates %s and %s never reached 1 Rq \n',gate2{i},gate1{i,1})];
        end
        smset(scan.loops(1).setchan,0);
    end
    if ~is4k && any(strcmp(name,{'QL','QR'}))
        global scandata
        scandata.config.hystVal = min(qocVal)-.15;
    end
    figure(10+j); clf;
    for i =1:length(dStore)
        subplot(2,1,1);
        plot(scanRng(scan),dStore{i,1},'DisplayName',gate2{i}); hold on;        
        subplot(2,1,2);
        plot(scanRng(scan),dStore{i,2},'DisplayName',gate2{i}); hold on;
    end
    xlabel('Voltage'); ylabel('Conductance');
    legend('show');
    subplot(2,1,1);    
    legend('show');
    xlabel('Voltage'); ylabel('Current');
    dStore = [];
    smset(channels,0);
end
fprintf('\n\n');  disp(logstr);
if ~is4k && zeroOL
    global scandata;
    scandata.config.closeVal = closeVal; 
end
cd ..
end