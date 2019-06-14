classdef Line < autotune.Op
   % Run scan to across junction to find tunnel coupling / temperature. 
    
    properties
        scan; % scan to take data
        fitFn = '@(p,x) p(1)+p(2)*(x)-p(4)*(1+tanh((x-p(3))./p(5)))/2'; %fit function for fitting
        perpDist = 0.55; % Distance parallel to jcn to scan. 0 = blTriple, 1 = trTriple. 0.5 = middle
        subPlot = 3; % for plotting in tuneData.figHandle
        nrep = 5; 
        range = 1e-3; % Range of scan 
        filePat = 'line'; 
        fineIndex = 1;
    end
    
    properties (SetAccess= {?autotune.Data, ?autotune.Op})
       width; % double, size Nrun 
    end
    
    methods         
        function this = Line            
            this.makeNewScan;             
        end
        
        function out = getData(this,runNumber)
            out = struct();
            out.width = this.width(runNumber);
        end
        
        function makeNewRun(this,runNumber)
            if runNumber ~= length(this.width)+1 
                warning('runNumber not consistent with know chrg runs');
            end
            this.width(runNumber) = nan;
            this.fineIndex = 1;
        end
        
        function run(this)
            global tuneData;
            side = upper(tuneData.activeSetName(1)); 
            file = sprintf('%s/sm_%s%s_%04i_%03i', tuneData.dir, this.filePat,side,tuneData.runNumber,this.fineIndex);
            %this.scan.consts = tuneData.chrg.scan.consts;
            %smset({this.scan.consts.setchan}, [this.scan.consts.val]);
            % Permute this away from the center a bit to avoid the mousebite.
            %       [x,y,x,y]
            chrg = [tuneData.chrg.blTriple(tuneData.runNumber,:),tuneData.chrg.trTriple(tuneData.runNumber,:)];
            cntr = [chrg(1)*this.perpDist+chrg(3)*(1-this.perpDist), chrg(2)*this.perpDist+chrg(4)*(1-this.perpDist)];
            this.scan.loops(1).trafofn(1).args = {cntr(1)};
            this.scan.loops(1).trafofn(2).args = {cntr(2)};
            this.fineIndex = this.fineIndex+1;
            data = smrun(this.scan, file); data = data{1}; 
            sleep('fast');
            this.ana('',data,this.scan);            
        end
        
        function ana(this,opts,data,scan)  
            global tuneData;   
            if ~exist('opts','var'), opts = ''; end 
            runNumber = tuneData.runNumber; 
            if ~exist('data','var'), data = []; end
            [data,out] = loadTunes(data,opts,this.filePat);                 
            if isempty(data), return; end            
            if ~isfield(out,'scan'), out.scan = scan; end
            
            eps = scanRng(out.scan,1); 
            data = 1e3*nanmean(data);
            slp=median(smooth(diff(data)))/(eps(2)-eps(1)); % do a linear fit for initial guess
            meanData=mean(data);     
            ign = 5; % Ignore first and last ign data points. 
            dataFlt = data-slp(1)*eps; % Look at the skew of the derivative to guess a sign for the atan.
            dataDiff = abs(smooth(diff(dataFlt)));            
            [~, maxDiff] = max(dataDiff(ign:end-ign)); % Find the step
            stepCent=eps(maxDiff+ign);
            % Get sign of slope
            if (mean(dataFlt(1:maxDiff+ign)) < mean(dataFlt(maxDiff+ign:end)))
                stepSlp=-range(dataFlt); %#ok<*CPROPLC>
            else
                stepSlp=range(dataFlt);
            end
            %       offset                           % step width   
            beta0 = [meanData,slp, stepCent, stepSlp, range(eps)/16.0];
            if ~any(isgraphics(tuneData.axes,'axes')), tuneData.rePlot([],'fig'); end             
            axes(tuneData.axes(this.subPlot)); cla; 
            fitPars = fitwrap('plfit noplot woff fine samefig', eps, data, beta0,this.fitFn, [1 1 0 0 0]);
            fitPars = fitwrap('plinit plfit woff fine samefig', eps, data, fitPars,this.fitFn);
            title(sprintf('Line: %g mV',fitPars(5)*1e3));
            a = gca; a.YTickLabelRotation=-30;        
            this.width(runNumber) = fitPars(5);
        end
        
        function makeNewScan(this)
            global tuneData
            scan = defScan;  %#ok<*PROP,*PROPLC>
            scan.loops(2).npoints = this.nrep;
            scan.loops(1).ramptime = -5e-3;
            scan.loops(1).rng = [-this.range this.range];
            scan.loops(1).npoints = 100;
            scan.consts(2).val = awgseqind(tuneData.chrg.pls);
            scan.loops(2).getchan = tuneData.dataChan;
            scan.loops(1).setchan = tuneData.xyChan;
            scan.loops(1).trafofn(1).fn = @(x,y,cntr) x(1)+cntr;
            scan.loops(1).trafofn(2).fn = @(x,y,cntr) -x(1)+cntr;            
            this.scan = scan;
        end
    end    
end