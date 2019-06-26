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
        fineIndex = 1; 
        filePat = 'loadTime';
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
            this.time(runNumber) = nan;
            this.amp(runNumber) = nan;
            this.fineIndex = 1;
        end
        
        function run(this)
            global tuneData;
            if ~awgcntrl('ison')
                awgcntrl('on start wait err')
            end
            file = sprintf('%s/sm_loadTime%s_%04i_%03i',tuneData.dir, upper(tuneData.activeSetName(1)), tuneData.runNumber,this.fineIndex);
            scan=fConfSeq(this.plsgrp,struct('nloop',this.nLoop,'nrep',this.nRep,'datachan',tuneData.dataChan,'opts','ampok'));
            scan = measAmp(scan);
            this.fineIndex = this.fineIndex+1;
            data = smrun(scan, file);
            if any(isnan(data{1}(:))); return; end
            if ndims(data{1}) == 3 %
                data = -diff(squeeze(mean(data{1})));
            else
                data = data{1};
            end
            this.ana('',data,scan);
        end
        
        function ana(this, opts,data,scan)
            global tuneData;
            if ~exist('opts','var'),  opts = ''; end            
            if ~exist('data','var'), data = []; end
            [data,out] = loadTunes(data,opts,this.filePat);                 
            if isempty(data), return; end
            if ~isfield(out,'scan'), out.scan = scan; end
            
            data = 1e3*nanmean(data); 
            tms = out.scan.data.pulsegroups.varpar(:,1)' * 1e-3; %fix me!                             
            
            beta0 = [min(data), range(data), .01];
            axes(tuneData.axes(this.subPlot)); cla; 
            try                
                params = fitwrap('woff plfit samefig', tms,data, beta0,this.fitFn);            
            catch
                params = nan(size(beta0)); 
            end
            tm = params(3); ampl = params(2); 
            if tm > 10
                fprintf('Load time didn''t fit or too high too measure \n'); 
                tm = nan; 
            end
            if ~out.anaData
                this.time(tuneData.runNumber) = tm;
                this.amp(tuneData.runNumber) = ampl;
            end
            title(sprintf('Load: %g ns',tm*1e3));        
            a = gca; a.YTickLabelRotation=-30;            
            a.XLim = [min(tms),max(tms)]; 
        end
        
        function updateGroup(this,config)
            %function updateGroup(this,config)
            % if config is a pulsegroup struct is uses that.
            % otherwise: remake loadTime group from dictionary
            global tuneData; 
            if ~exist('config','var'), config = []; end
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
                pg.chan=[getNum(tuneData.xyChan{1}),getNum(tuneData.xyChan{2})];
                pg.pulses = 6;
                pg.dict=tuneData.activeSetName;
                
                dict=pdload(pg.dict);
                cntr = dict.reload.val;
                pg.name = ['loadTime_1_', upper(tuneData.activeSetName(1))];
                %       1
                %params=[ramp to/from load (ns), loadTime (ns), load cntr loadpos offset(mV)]
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

