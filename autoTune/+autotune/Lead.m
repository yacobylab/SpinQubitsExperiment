classdef Lead < autotune.Op
    %LoadTime < autotune.Op % calculate the speed of the singlet load
    %   see help for individual properties and methods.
    %   e.g. >> help autotune.LoadTime.plsgrp
    %   e.g. >> help autotune.LoadTime.Run
    % 5 mV pulse in the X and Y directions. 
    % Y lead = one with larger slope (cross it by moving in x-direction)
    % v.v. for X lead. 
    % 1: black green, 2: cyan magenta. 
    % Needs the leads to be marked in charge scan to properly run scan
    properties
        scan; 
        plsGrp = {'sqrX_L', 'sqrY_L'}; %pulsegroup to take data        
        nRep = 10;  % for scan
        subPlot = 1:2;  
        samprate = 1.1e8;
    end
    
    properties (SetAccess= {?autotune.Data, ?autotune.Op})        
        timeX;
        timeY; 
        posX;
        posY;
    end
    
    methods
        function this = Lead
            global tuneData
            if strcmp(tuneData.activeSetName,'right')
                this.plsGrp = {'sqrX_R','sqrY_R'};
            end
            this.makeNewScan;
        end
        
        function out = getData(this,runNumber)
            out = struct();
            out.time = this.time(runNumber);
            out.amp = this.amp(runNumber);
        end
        
        function makeNewRun(this,runNumber)
            if runNumber ~= length(this.timeX)+1
                warning('runNumber not consistent with know chrg runs');
            end
            this.timeX(end+1,:) = nan;
            this.posX(end+1,:) = nan;
            this.timeY(end+1,:)=nan;
            this.posY(end+1,:)=nan;
        end        
        
        function run(this, runNumber)
            global tuneData; global smdata;
            if ~exist('runNumber','var')
                runNumber = tuneData.runNumber;
            end
            sqrAmp=1.5e-3; %hard coded square wave amplitude (pk-pk). This appears to be set in the pulsedef.
            nLeads = 2;
            figure(77); clf;
            for i = 1:nLeads
                file = sprintf('%s/sm_lead%d%s_%04i',tuneData.dir, i, upper(tuneData.activeSetName(1)), runNumber);
                if i == 1 % Y lead, x square wave [larger slope, affected by X gate]. 
                    leadSlp = [1; tuneData.chrg.yLeadSlope(runNumber)]; % slope of vertical Left lead (lower left).
                    leadSlp = leadSlp/norm(leadSlp);
                    % move down lead by amount sqrAmp, then across in xdir, starting 1/2 sqr amp over.
                    leadPos1 = tuneData.chrg.blTriple(runNumber,:)' +2*sqrAmp * leadSlp; %get the right offset
                    %leadPos1 = leadPos1+[-sqrAmp/4;0];
                    leadPos2 = leadPos1 + [sqrAmp;0];
                else % x lead, y square wave [smaller slope, affected by y gate]. 
                    leadSlp = [1; tuneData.chrg.xLeadSlope(runNumber)]; % slope of horizonal left lead. (upper left)
                    %leadSlp = leadSlp / norm(leadSlp);
                    leadPos1 = tuneData.chrg.blTriple(runNumber,:)' - sqrAmp * leadSlp; %get the right offset
                    %leadPos1 = leadPos1+[0; -sqrAmp/4];
                    leadPos2 = leadPos1 + [0; sqrAmp];
                end
                this.scan.loops(2).trafofn(1).fn = @(x,y,leadPos) leadPos(x(2));
                this.scan.loops(2).trafofn(1).args = {[leadPos1(1),leadPos2(1)]};
                this.scan.loops(2).trafofn(2).fn = @(x,y,leadPos) leadPos(x(2));
                this.scan.loops(2).trafofn(2).args = {[leadPos1(2),leadPos2(2)]};
                this.scan.loops(2).prefn(2).args = {awgseqind(this.plsGrp{i})}; % set pulsegroup for scan
                inst = smchaninst(tuneData.dataChan);
                %smdata.inst(inst(1)).data.extclk = 2;
                clearMask(tuneData.dataChan)
                data = smrun(this.scan, file);
                sleep('fast'); 
                %smdata.inst(inst(1)).data.extclk = 0;
                smdata.inst(inst(1)).data.combine = [];
                if any(isnan(data{1}(:))); return; end
                this.ana('',data{1},i);
            end
        end
        
        function ana(this,opts,data,i)
            global tuneData;
            if ~exist('opts','var'), opts = ''; end
            % Load data
            if ~exist('data','var') || isempty(data) || ischar(data) || numel(data)==1
                if (~exist('data','var') || isempty(data)) &&~isopt(opts,'last') %&&~ischar(data)
                    [data,scan,file]=loadAna('sm_lead*');
                    i=str2double(file(end-10));
                elseif exist('data','var') && ~isempty(data) && ischar(data)
                    [data,scan]=loadAna(data);
                    if isempty(data)
                        return
                    end
                else
                    side = upper(tuneData.activeSetName(1));
                    if isopt(opts,'last')
                        data = tuneData.runNumber;
                    end
                    fileName = sprintf('sm_lead1%s_%04.f.mat',side,data);
                    this.ana('',fileName,1); 
                    fileName = sprintf('sm_lead2%s_%04.f.mat',side,data);
                    i = 2; 
                    [data,scan]=loadAna(fileName);
                    if isempty(data)
                        return
                    end
                end
            else
                scan = this.scan;
            end
            samprate = scan.consts(1).val;
            leadName = {'Y', 'X'}; 
            colorList = {'black and green','magenta and cyan'}; 
            nleads = 2; 
            data = squeeze(mean(data))'; 
            t = (0:length(data)-1)./samprate * 1e6; %#ok<*PROPLC> % given in us
            figure(77); subplot(nleads,1,i);                        
            plot(t,data); hold on; 
            xlabel('Time (us)'); title(sprintf('Lead %s, color %s',leadName{i},colorList{i})); 
            data = diff(data,[],2); %subtract off opposite direction.                                         figure(3); subplot(3,3,this.subPlot(i));            
            
            slopeData = range(data)/6*sign(data(round(end/4))-data(round(3*end/4))); 
            %      offset ,     slope,   width of left, right, wrap around (in case pulse timing not quite synced)
            beta0 = [mean(data), slopeData, .35, .35, .1];
                                
            axes(tuneData.axes(this.subPlot(i)));  
            params = fitwrap('plinit plfit samefig', t, data', beta0, @leadfn);            
            params = fitwrap('plinit plfit samefig', t, data', params, @leadfn);            
            runNumber = tuneData.runNumber;
            if i == 1
                this.timeX(runNumber,1:2) = params(3:4);
                title(sprintf('Lead %s',leadName{i}));
                fprintf('Lead %s: %g us, %g us\n',leadName{i},this.timeX(runNumber,1),...
                    this.timeX(runNumber,2));
            else
                this.timeY(runNumber,1:2) = params(3:4);
                title(sprintf('Lead %s',leadName{i}));
                fprintf('Lead %s: %g us, %g us\n',leadName{i},this.timeY(runNumber,1),...
                    this.timeY(runNumber,2));
            end
            a = gca; a.YTickLabelRotation=-30; a.XLim = [min(t),max(t)]; 
            a.YLabel.Position(1) = a.XLim(1) - range(a.XLim)/14;
            a.XLabel.Position(2) = a.YLim(1) - range(a.YLim)/7;
            figure(1); hold on; 
            % Plot the places where lead scans are centered on the charge diagram. 
            trafofn = this.scan.loops(2).trafofn; 
            if i == 1
                plot(trafofn(1).args{1}(1),trafofn(2).args{1}(1),'.k','MarkerSize',12)
                plot(trafofn(1).args{1}(2),trafofn(2).args{1}(2),'.g','MarkerSize',12)
            else
                plot(trafofn(1).args{1}(1),trafofn(2).args{1}(1),'.m','MarkerSize',12)
                plot(trafofn(1).args{1}(2),trafofn(2).args{1}(2),'.c','MarkerSize',12)
            end
        end
        
        function makeNewScan(this)
            global tuneData; global smdata;
            this.scan=[];
            this.scan.saveloop = [3 20];
            this.scan.consts(1).setchan = 'samprate';
            this.scan.consts(1).val =this.samprate;            
            this.scan.configfn(2).fn = @setzero; % Set voltages to 0 at end
            this.scan.configfn(2).args = {};
            
            plsLength=4e-6;
            this.scan.loops(1).ramptime = 0.08;
            npoints = this.samprate*this.scan.loops(1).ramptime; 
            npointsPulse = plsLength*this.samprate;  
            % -> downsamp will be npointsPulse here. 
            this.scan.loops(1).npoints = 1; 
            
            % Ramp in different X direction on first point, Y on second. 
            this.scan.loops(2).npoints = 2;
            this.scan.loops(2).setchan = tuneData.xyChan; 
            this.scan.loops(2).getchan = {tuneData.dataChan,'Time'}; 
            daqInst = smchaninst(tuneData.dataChan);
            this.scan.loops(2).prefn(2).fn = @(x,n) smset('PulseLine',n); 
            this.scan.loops(2).prefn(2).args = {awgseqind(this.plsGrp(1))};
            fn = sprintf('@(x,ico) %s([ico 4])',func2str(smdata.inst(daqInst(1)).cntrlfn));            
            this.scan.loops(2).prefn(1).fn = fn;
            this.scan.loops(2).prefn(1).args = {daqInst}; 
            
            % Repeat for averaging in first loop 
            this.scan.loops(3).npoints = this.nRep;
            this.scan.loops(3).setchan = 'count';
            
            fn = sprintf('@(ico,npoints,samprate,nsm,ctrl) %s([ico 5],npoints,samprate,nsm,ctrl)',func2str(smdata.inst(daqInst(1)).cntrlfn));
            this.scan.configfn(1).fn = @smaconfigwrap; 
            this.scan.configfn(1).args = {fn,daqInst,npoints/npointsPulse,1/plsLength,'mean',npointsPulse};            
            this.scan.loops(1).stream = 1; 
            this.scan.loops(2).stream = 1; 
        end
              
    end
end

function y = leadfn(p, t)
% Fit function for lead
% x = time
tmid = t(end/2+1); %mean time.
t = mod(t-p(5), 2*tmid);
y1 = cosh(tmid/2/p(3))-exp((tmid/2-t)./p(3))./sinh(tmid/2/p(3)); 
y2 = cosh(tmid/2/p(4))-exp((3*tmid/2-t)./p(4))./sinh(tmid/2/p(4)); 
y =  p(1)+.5*p(2)*(y1.*(t < tmid)-y2.*(t >= tmid)); %function to fit data
end