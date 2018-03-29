classdef LoadTime < autotune.Op
    %LoadTime < autotune.Op % calculate the speed of the singlet load
    %   see help for individual properties and methods.
    %   e.g. >> help autotune.LoadTime.plsgrp
    %   e.g. >> help autotune.LoadTime.Run
    
    properties
        plsgrp = 'loadTime_1_L'; %pulsegroup to take data
        fitFn = '@(p, x)p(1)+p(2)*exp(-x./p(3))'; %function to fit data
        nRep = 20;  % for scan
        nLoop = 1000;% for scan
        subPlot = 4; %for plotting in tuneData.figHandle
    end
    
    properties (SetAccess= {?autotune.Data, ?autotune.Op})
        time; %characteristic load time double size Nrun)
        amp; %amplitude in mV of the load scan, double size Nrun)
    end
    
    methods        
        function this = LoadTime
            global tuneData
            if strcmp(tuneData.activeSetName,'right')
                this.plsgrp ={'loadTime_1_R'};
            end            
        end
        
        function out = getData(this,runNumber)
            out = struct();
            out.time = this.time(runNumber);
            out.amp = this.amp(runNumber);
        end
        
        function makeNewRun(this,runNumber)
            if runNumber ~= length(this.time)+1 
                warning('runNumber not consistent with know chrg runs');
            end
            this.time(end+1) = nan;
            this.amp(end+1) = nan;
        end
        
        function run(this)
            global tuneData;
            if ~awgcntrl('ison')
                awgcntrl('on start wait err')
            end
            file = sprintf('%s/sm_loadTime%s_%04i',tuneData.dir, upper(tuneData.activeSetName(1)), tuneData.runNumber);
            scan=fConfSeq(this.plsgrp,struct('nloop',this.nLoop,'nrep',this.nRep,'datachan',tuneData.dataChan,'opts','ampok'));%,'hwsampler',100e6));
            scan.consts(1).setchan=tuneData.xyChan(1);
            scan.consts(2).setchan=tuneData.xyChan(2);
            scan.consts(1).val = tuneData.measPt(1);
            scan.consts(2).val = tuneData.measPt(2);
            scan.loops(1).stream = 1; 
            data = smrun(scan, file);
            if any(isnan(data{1}(:))); return; end
            if ndims(data{1}) == 3 %
                data = -diff(squeeze(mean(data{1})));
            else
                data = mean(data{1});
            end
            this.ana('',data,scan);
        end
        
        function ana(this, opts,data,scan)
            global tuneData;
            if ~exist('opts','var')
                opts = '';
            end
            runNumber = tuneData.runNumber;
            if ~exist('data','var') || isempty(data) || ischar(data) || numel(data)==1
                if (~exist('data','var') || isempty(data)) &&~isopt(opts,'last')
                    [data,scan]=loadAna('sm_loadTime*');
                elseif exist('data','var') && ~isempty(data) && ischar(data)
                    [data,scan]=loadAna(data);
                else
                    side = upper(tuneData.activeSetName(1));
                    if isopt(opts,'last')
                        data = tuneData.runNumber;
                    end
                    fileName = sprintf('sm_loadTime%s_%04.f.mat',side,data);
                    [data,scan]=loadAna(fileName);
                    if isempty(data)
                        return
                    end
                end
                anaData=1;
            else
                anaData=0;
            end
            eps = scan.data.pulsegroups.varpar(:,1)' * 1e-3; %fix me!
            if ischar(this.fitFn)
                func = str2func(this.fitFn);
            else
                func = this.fitFn;
            end            
            axes(tuneData.axes(this.subPlot)); 
            beta0 = [min(data), range(data), .01];
            pars = fitwrap('woff plinit plfit samefig', eps,data, beta0,func);
            if ~anaData
                this.time(runNumber) = pars(3);
                this.amp(runNumber) = pars(2);
            end
            title(sprintf('Load: %g',this.time(runNumber)*1e3));        
            a = gca; a.YTickLabelRotation=-30;            
            a.XLim = [min(eps),max(eps)]; 
            a.YLabel.Position(1) = a.XLim(1) - range(a.XLim)/14;
            a.XLabel.Position(2) = a.YLim(1) - range(a.YLim)/7;
        end
        function updateGroup(this,config)
            %function updateGroup(this,config)
            % if config is a pulsegroup struct is uses that.
            % otherwise: remake loadPos group from dictionary
            global tuneData; 
            if ~exist('config','var')
                config = [];
            end
            if isstruct(config) %regular pulse group struct
                if isfield(config,'name') && ~strcmp(config.name,this.plsGrp)
                    error('pg.name = %s, loadTime plsGrp.name = %s\n',config.name,this.plsGrp);
                end
                try
                    plsupdate(config);
                    awgadd(config.name);
                    awgcntrl('on start wait err');
                catch
                    fprintf('problem making groups %s. quitting...\n',config.name);
                end
                return;
            elseif isempty(config) %make default group;
                pg.chan=[str2double(char(regexp(tuneData.xyChan{1},'\d+','match'))),str2double(char(regexp(tuneData.xyChan{2},'\d+','match')))];                dict = 'left';
                pg.pulses = 6;
                pg.dict=tuneData.activeSetName;
                
                dict=pdload(pg.dict);
                cntr = dict.reload.val;
                pg.name = ['loadTime_1_', upper(tuneData.activeSetName(1))];
                pg.params = [20 0 0 0 0];
                pg.varpar(:,2:3)= repmat(cntr', 1,51)';
                pg.varpar(:,1)=(0:5:250)';
                pg.varpar(:,4)=0;
                
                pg.ctrl = 'loop pack';
                plsdefgrp(pg); 
                awgadd(pg.name);
            else
                error('must past pulsegroup struct or empty to updateGroup')
            end
        end
        
    end
end

