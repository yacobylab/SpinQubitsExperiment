classdef LoadPos < autotune.Op
    %LoadPos < autotune.Op % represents the load position scan
    %   see help for individual properties and methods.
    %   e.g. >> help autotune.LoadPos.plsGrp
    %   e.g. >> help autotune.LoadPos.Run
    
    properties
        plsGrp = 'loadPos_1_L'; %plulse group used in scan
        fitFn = '@(p,x) p(7) + p(1)*(tanh((x-p(3))/p(2))+1)/2+ p(4)*(tanh((x-p(6))/p(5))+1)/2'; %fit function 
        nRep = 20; %for scan
        nLoop = 500; %for scan
        subPlot = 5; %for plot in tuneData.figHandle
        pulseScan; 
        target; 
    end
    
    properties (SetAccess= {?autotune.Data, ?autotune.Op})
        width; %characteristic width of load region, (double size Nrun)
    end
    
    methods        
        function this = LoadPos
            global tuneData
            if strcmp(tuneData.activeSetName,'right')
                this.plsGrp ={'loadPos_1_R'};
            end           
        end
        
        function out = getData(this,runNumber)
            %function out = getData(this,runNumber)
            % returns a struct with the width from run runNumber
            out = struct();
            out.width = this.time(runNumber);
        end
        
        function makeNewRun(this,runNumber)
            %function makeNewRun(this,runNumber)
            % make a new run- populate this.width(end+1) = NaN;
            if runNumber ~= length(this.width)+1
                warning('runNumber not consistent with know chrg runs');
            end
            this.width(end+1) = nan;
        end
        
        function run(this)
            %function run(this,runNumber)
            %run the load position scan
            global tuneData;
            if ~awgcntrl('ison')
                awgcntrl('on start wait err');
            end
            file = sprintf('%s/sm_loadPos%s_%04i',tuneData.dir,upper(tuneData.activeSetName(1)), tuneData.runNumber);
            scan=fConfSeq(this.plsGrp,struct('nloop',this.nLoop,'nrep',this.nRep,'datachan',tuneData.dataChan,'opts','ampok'));%,'hwsampler',100e6));
            scan.consts(1).setchan=tuneData.xyChan(1);
            scan.consts(2).setchan=tuneData.xyChan(2);
            scan.consts(1).val = tuneData.measPt(1);
            scan.consts(2).val = tuneData.measPt(2);
            scan.loops(1).stream = 1; 
            
            data = smrun(scan, file);
            if any(isnan(data{1}(:))); return; end
            data = data{1};
            this.ana('',data,scan); 
        end
        
        function out=ana(this,opts,data,scan)
            % function out=ana(this,opts,data,scan)
            global tuneData
            if ~exist('opts','var'), opts = '';       end
            runNumber = tuneData.runNumber;
            if ~exist('data','var') || isempty(data) || ischar(data) || numel(data)==1 % Check if loading old scan or analyzing new data.
                if (~exist('data','var') || isempty(data)) &&~isopt(opts,'last')
                    [data,scan,~,time]=loadAna('sm_loadPos*');
                    out.time = time;
                elseif exist('data','var') && ~isempty(data) && ischar(data)
                    [data,scan,~,time]=loadAna(data);
                    out.time = time;
                else
                    side = upper(tuneData.activeSetName(1));
                    if isopt(opts,'last'), data = tuneData.runNumber; end
                    fileName = sprintf('sm_loadPos%s_%04.f.mat',side,data);
                    [data,scan,~,time]=loadAna(fileName);
                    if isempty(data)
                        out = struct;
                        return
                    else
                        out.time = time;
                    end
                end
                anaData=1;
            else
                anaData=0;
            end
            if ndims(data) == 3 % ??
                data = -diff(squeeze(mean(data)));
            else
                data = mean(data);
            end
            if ischar(this.fitFn)
                func = str2func(this.fitFn);
            else
                func = this.fitFn;
            end
            xv = scan.data.pulsegroups.varpar';
            [~,loadInd] = min(data);
            this.target = xv(loadInd);
            %        tanh amp, 1st tanh width, 1st step loc,
            beta0 = [range(data), .1*range(xv),xv((diff(data) == max(diff(data)))), range(data),-.1*range(xv),xv((diff(data) == min(diff(data)))), min(data)];
            axes(tuneData.axes(this.subPlot)); cla;
            params = fitwrap('plinit plfit samefig woff', xv, data, beta0, func, [1 1 1 1 1 1 1]);
            if ~anaData
                this.width(runNumber) = abs(params(3)-params(6));                
            end
            title(sprintf('Load wdth: %3.1f uV',1e3*this.width(runNumber)));
            a = gca; a.YTickLabelRotation=-30;
            a.XLim = [min(xv),max(xv)];
            a.YLabel.Position(1) = a.XLim(1) - range(a.XLim)/18;
            a.XLabel.Position(2) = a.YLim(1) - range(a.YLim)/7;
            eps = scan.data.pulsegroups.varpar'; % tl scan sweeps epsilon value (along TL curve)
            params = scan.data.pulsegroups.params; 
            loadCenter = params(3:4); 
            
            figure(tuneData.chrg.figHandle); hold on;
            plot(tuneData.measPt(1)+loadCenter(1)*1e-3+this.target*1e-3,tuneData.measPt(2)+loadCenter(2)*1e-3+this.target*1e-3,'kx');
            loadRng=[loadCenter+[1,1]*max(eps); loadCenter+[1,1]*min(eps)];
            plot(tuneData.measPt(1)+loadRng(:,1)*1e-3,tuneData.measPt(2)+loadRng(:,2)*1e-3,'g-');
            for i=1:2
                for j = 1:length(scan.consts)
                    if strcmp(tuneData.xyChan{i},scan.consts(j).setchan)
                        out.measPt(i) = scan.consts(j).val;
                    end
                end
            end
        end
        
        function updateGroup(this,opts,config)
            global tuneData;
            %function updateGroup(this,opts,config)
            % using the 'target' opt, uses the this.target value to update
            % the load position.
            % if config is a pulsegroup struct is uses that.
            % Parameters you might want to control are ramptime, loadtime,
            % load center, range of scan.
            % otherwise: remake loadPos group from dictionary
            if ~exist('opts','var'), opts = ''; end
            if ~exist('config','var'), config = []; end
            if isstruct(config) %regular pulse group struct
                if isfield(config,'name') && ~strcmp(config.name,this.plsGrp)
                    error('pg.name = %s, loadPos plsGrp.name = %s\n',config.name,this.plsGrp);
                end
                try
                    plsupdate(config);
                    awgadd(config.name);
                    awgcntrl('on start wait err');
                catch
                    fprintf('Problem making groups %s. quitting...\n',config.name);
                end
                return;
            elseif isempty(config) %make default group;
                if isopt(opts,'target')
                    dict = pdload(tuneData.activeSetName);
                    dict.reload.val = dict.reload.val + [1,1]*this.target;
                    pdsave(tuneData.activeSetName,dict);
                    tuneData.updateAll('nodict');
                    return
                end
                pg.chan=[str2double(char(regexp(tuneData.xyChan{1},'\d+','match'))),str2double(char(regexp(tuneData.xyChan{2},'\d+','match')))];
                rangeScale = 3;
                pg.pulses = 6;
                pg.dict=tuneData.activeSetName;
                dict=pdload(pg.dict);
                cntr = dict.reload.val;
                %params=[ramp to/from load (ns), loadTime (ns), load cntr loadpos offset(mV)]
                pg.params = [20 500 cntr 0];
                pg.varpar = (-.5:.01:.5)' * rangeScale; %[range of scan from the current reload val               
                pg.name = ['loadPos_1_' upper(tuneData.activeSetName(1))];
                pg.ctrl = 'loop pack';
                plsdefgrp(pg);
                awgadd(pg.name);
            else
                error('must past pulsegroup struct or empty to updateGroup')
            end
        end
        
        function runTwoD(this,opts)
            global tuneData;
            if ~exist('opts','var'), opts = ''; end
            if ~isopt(opts,'run')
                pg.chan=[str2double(char(regexp(tuneData.xyChan{1},'\d+','match'))),str2double(char(regexp(tuneData.xyChan{2},'\d+','match')))];
                rangeScale = 6;
                pg.pulses = 6;
                pg.dict=tuneData.activeSetName;
                dict=pdload(pg.dict);
                cntr = dict.reload.val;
                %params=[ramp to/from load (ns), loadTime (ns), load cntr loadpos offset(mV)]
                pg.params = [20 500 cntr 0];
                pg.varpar = (-.5:.01:.5)' * rangeScale; %[range of scan from the current reload val
                pg.ctrl = 'loop pack';
                centVal = linspace(-rangeScale/2,rangeScale/2,10);
                yVal = pg.params(4);
                pg.params(3) = pg.params(3)+1.5; 
                for i =1:10
                    pg.params(4) = yVal+centVal(i);
                    pg.name = sprintf('loadPos_%d_%s',i, upper(tuneData.activeSetName(1)));
                    plsdefgrp(pg);
                    loadGrp{i}=pg.name;
                end
                awgadd(loadGrp);
                awgcntrl('on start wait err')                
                this.pulseScan = loadGrp; 
            end            
            scan=fConfSeq(this.pulseScan,struct('nloop',this.nLoop,'nrep',this.nRep,'datachan',tuneData.dataChan,'opts','ampok'));%,'hwsampler',100e6));
            smrun(scan,smnext(sprintf('Load2D_%s',upper(tuneData.activeSetName(1)))));
        end
    end    
end