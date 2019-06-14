classdef LoadPos < autotune.Op
    %LoadPos < autotune.Op % represents the load position scan
    %   see help for individual properties and methods.
    %   e.g. >> help autotune.LoadPos.plsGrp
    %   e.g. >> help autotune.LoadPos.Run
    
    properties
        plsGrp = 'loadPos_1_L'; %pulse group used in scan
        fitFn = '@(p,x) p(7) + p(1)*tanh((x-p(3))/p(2))/2 + p(4)*tanh((x-p(6))/p(5))/2';
        nRep = 10; % for scan
        nLoop = 500; % for scan
        subPlot = 5; % for plot in tuneData.figHandle
        pulseScan; 
        target =0;
        location; 
        rangeScale = 3;
        dist = 2; 
        slope = -5; % this is likely to be the same as tl slope, the x/y lead slope
        fineIndex = 1; 
        filePat = 'loadPos';
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
            this.width(runNumber) = nan;
            this.fineIndex = 1;
        end
        
        function run(this)
            %function run(this,runNumber)
            %run the load position scan
            global tuneData;
            if ~awgcntrl('ison')
                awgcntrl('on start wait err');
            end
            side = upper(tuneData.activeSetName(1)); 
            file = sprintf('%s/sm_%s%s_%04i_%03i',tuneData.dir,this.filePat, side, tuneData.runNumber,this.fineIndex);
            scan=fConfSeq(this.plsGrp,struct('nloop',this.nLoop,'nrep',this.nRep,'datachan',tuneData.dataChan,'opts','ampok'));
            scan = measAmp(scan); % Add measPt voltage to scan. 
            this.fineIndex = this.fineIndex+1;
            data = smrun(scan, file);
            if any(isnan(data{1}(:))); return; end
            data = data{1};
            this.ana('',data,scan);
        end
        
        function out=ana(this,opts,data,scan)
            % function out=ana(this,opts,data,scan)
            global tuneData
            if ~exist('opts','var'), opts = ''; end
            runNumber = tuneData.runNumber;
            % Check if loading old scan or analyzing new data.
            if ~exist('data','var'), data = []; end
            [data,out] = loadTunes(data,opts,this.filePat);                 
            if isempty(data), return; end
            if ~isfield(out,'scan'), out.scan = scan; end
            data = data*1e3; 
            if ndims(data) == 3 % ?? oh myabe for two d data
                data = -diff(squeeze(mean(data)));
            else
                data = mean(data);
            end            

            xv = out.scan.data.pulsegroups.varpar';             
            loadCenter = out.scan.data.pulsegroups.params(3:4);
            xvX = xv(1,:)+loadCenter(1); % Use only x coord's for fitting. 
            [minData,loadInd] = min(data);
            
            this.target = xv(:,loadInd)'+loadCenter; % Add load center to target
            [~,posStep]=max(diff(data)); 
            [~,negStep]=min(diff(data));             
            %       tanh amp, 1st tanh width, step loc,   tanh amp, 2nd tanh width, 2nd loc offset,  
            beta0 = [range(data), 0.1*range(xvX),xvX(posStep),range(data),-0.1*range(xvX),xvX(negStep), min(data)];
            if ~any(isgraphics(tuneData.axes,'axes')), tuneData.rePlot([],'fig'); end             
            axes(tuneData.axes(this.subPlot)); cla;             
            
            params = fitwrap('plfit samefig woff', xvX, data, beta0, this.fitFn);
            plot(xvX(floor(end/2)),data(floor(end/2)),'o'); % Plot o for center of scan (often current load point). 
            plot(this.target(1),minData,'kx'); % Plot o for center of scan (often current load point). 
            width = abs(params(3)-params(6));
            if ~out.anaData
                this.width(runNumber) = width; %#ok<*PROPLC>
            end
            title(sprintf('Load wdth: %3.1f uV',1e3*width));
            a = gca; a.YTickLabelRotation=-30;
            a.XLim = [min(xvX),max(xvX)];
            
            % Plotting 
            figure(tuneData.chrg.figHandle); hold on;
            % Plot x on charge diagram with best charge point. 
            plot(tuneData.measPt(1)+this.target(1)*1e-3,tuneData.measPt(2)+this.target(2)*1e-3,'kx');
            loadRng=[loadCenter+xv(:,1)'; loadCenter+xv(:,end)'];
            % Plot line showing range of load scan, starting with white point
            plot(tuneData.measPt(1)+loadRng(:,1)*1e-3,tuneData.measPt(2)+loadRng(:,2)*1e-3,'g-');
            plot(tuneData.measPt(1)+loadRng(1,1)*1e-3,tuneData.measPt(2)+loadRng(1,2)*1e-3,'w.');
            plot(tuneData.measPt(1)+loadRng(end,1)*1e-3,tuneData.measPt(2)+loadRng(end,2)*1e-3,'k.');
            for i=1:2
                for j = 1:length(out.scan.consts)
                    if strcmp(tuneData.xyChan{i},out.scan.consts(j).setchan)
                        out.measPt(i) = out.scan.consts(j).val;
                    end
                end
            end
        end
        
        function updateGroup(this,opts,config)
            global tuneData;
            % function updateGroup(this,opts,config)
            % using the 'target' opt, uses the this.target value to update
            % the load position.
            % if config is a pulsegroup struct is uses that.
            % Parameters you might want to control are ramptime, loadtime,
            % load center, range of scan.
            % otherwise: remake loadPos group from dictionary
            if ~exist('opts','var'), opts = ''; end
            if ~exist('config','var'), config = []; end
            
            if tuneData.sepDir(1)>0  % TL meas point. Set up direction of lead.
                leadDir = -[1 this.slope];
                loadDir = [1 1./leadDir(2)]; loadDir=-sqrt(2)*loadDir/(norm(loadDir));
            else  % BR meas point
                leadDir = [1 this.slope];
                loadDir = [1 -1./leadDir(2)]; loadDir=sqrt(2)*loadDir/(norm(loadDir));
            end
            leadDir = leadDir / norm(leadDir); % lead Direction, towards the right            
                                   
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
                    dict.reload.val = this.target; % Sqrt 2?
                    pdsave(tuneData.activeSetName,dict);
                    tuneData.updateAll('nodict');
                    return
                end
                if isopt(opts,'init')
                    juncDist = tuneData.chrg.trTriple(end,:)-tuneData.chrg.blTriple(end,:);
                    loadCenter=-1e3*(1/2*juncDist+tuneData.chrg.defaultOffset) + leadDir*this.dist;
                else
                    dict = pdload(tuneData.activeSetName);
                    loadCenter = dict.reload.val;
                end
                pg.chan=[getNum(tuneData.xyChan{1}),getNum(tuneData.xyChan{2})];
                pg.pulses = 120;
                pg.dict=tuneData.activeSetName;
                
                %params=[ramp to/from load (ns), loadTime (ns), load cntr loadpos offset(mV)]
                pg.params = [1 500 loadCenter 0 0];
                pg.varpar = (-.5:.01:.5)' * this.rangeScale*loadDir-tuneData.measPt; %[range of scan from the current reload val
                pg.name = ['loadPos_1_' upper(tuneData.activeSetName(1))];
                pg.ctrl = 'loop pack';
                plsdefgrp(pg);
                awgadd(pg.name);
            else
                error('must pass pulsegroup struct or empty to updateGroup')
            end
        end
        
        function runTwoD(this,opts)
            global tuneData;
            if ~exist('opts','var'), opts = ''; end
            if ~isopt(opts,'run')
                pg.chan=[getNum(tuneData.xyChan{1}),getNum(tuneData.xyChan{2})];                
                pg.pulses = 120;
                pg.dict=tuneData.activeSetName;
                dict=pdload(pg.dict);
                cntr = dict.reload.val; % Start at current load point. 
                %params=[ramp to/from load (ns), loadTime (ns), load cntr loadpos offset(mV)]
                pg.params = [20 500 cntr 0 0];
                xVals = (-.5:.01:.5)' * this.rangeScale; %[range of scan from the current reload val
                pg.ctrl = 'loop pack';
                yVal = linspace(-this.rangeScale/2,this.rangeScale/2,10);
                for i =1:10
                    pg.varpar = [xVals, yVal(i)*ones(size(xVals))];                    
                    pg.name = sprintf('loadPos_%02d_%s',i, upper(tuneData.activeSetName(1)));
                    plsdefgrp(pg);
                    loadGrp{i}=pg.name;
                end
                awgadd(loadGrp);
                awgcntrl('on start wait err')
                this.pulseScan = loadGrp;
            end
            scan=fConfSeq(this.pulseScan,struct('nloop',this.nLoop,'nrep',this.nRep,'datachan',tuneData.dataChan,'opts','ampok'));%,'hwsampler',100e6));
            scan = measAmp(scan); 
            data=smrun(scan,smnext(sprintf('Load2D_%s',upper(tuneData.activeSetName(1)))));
            data = squeeze(nanmean(data{1}));
            figure(100); clf;
            imagesc(data);
            a = gca; a.YDir = 'normal';

        end
    end
end