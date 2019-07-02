classdef TwoSen < autotune.Op
    %Zoom < autotune.Op % represents a zoom scan to find Pauli Blockade
    %   see help for individual properties and methods.
    %   e.g. >> help autotune.Zoom.plsGrp
    %   e.g. >> help autotune.Zoom.Run
    
    properties
        scan; % scan for taking data
        subPlot = nan;
        fineIndex = 1;        
    end
    
    properties (Constant)
        figHandle = 2; %fot plotting zoom scan
    end
    
    methods
        function this = TwoSen            
            global tuneData
            
            scan = defScan('chrg',tuneData.activeSetName); %#ok<*PROPLC>
            scan.loops(1).npoints = 200;
            scan.loops(1).ramptime = -7e-3;
            scan.loops(1).rng = [-.35 -.45];                      
            scan.loops(2).npoints = 25;
            scan.loops(2).rng = [-.2 -.5];
            scan.loops(1).settle =0.25; 
            if strcmp(tuneData.activeSetName,'right')
                scan.loops(1).setchan = {'SD4top'};
                scan.loops(2).setchan = {'SD4bot'};
            else
                scan.loops(1).setchan = {'SD1top'};
                scan.loops(2).setchan = {'SD1bot'};
            end            
            this.scan = scan;           
        end
        
        function out = getData(this,runNumber)
        end
        
        function makeNewRun(this,runNumber)
            this.fineIndex = 1; 
        end
        
        function run(this,opts)
            % TwoD Scan.           
            global tuneData;
           if ~exist('opts','var'), opts = ''; end
            file = sprintf('%s/sm_SD%s_%04i_%03i', tuneData.dir,...
                upper(tuneData.activeSetName(1)), tuneData.runNumber, this.fineIndex);
            xvals = scanRng(this.scan,1);      
            d=smrun(tuneData.twoSen.scan,file);
            diffData = diff(d{1},1,2)/(xvals(2)-xvals(1)); % Take first order diff across row
            this.fineIndex = this.fineIndex+1;                   
            xDiff = (xvals(1:end-1)+xvals(2:end))/2;
            yvals = scanRng(this.scan,2);
            
            [maxVal,indX] = max(abs(diffData),[],2); 
            [maxDiff,indY] = max(abs(maxVal));
            indX = indX(indY);
            
            smset(this.scan.loops(1).setchan{1},xDiff(indX));
            smset(this.scan.loops(2).setchan{1},yvals(indY));
            fprintf('Setting gates to point of max sensitivity, %2.2f. %s to %4.4f. %s to %4.4f \n',...
                maxDiff, this.scan.loops(1).setchan{1}, xvals(indX), this.scan.loops(2).setchan{1}, yvals(indY));
        
        if isopt(opts,'fine')
             scan= this.scan; 
             this.scan.loops(1).rng = [xDiff(indX)-0.01,xDiff(indX) + 0.01]; 
             this.scan.loops(1).npoints = floor(this.scan.loops(1).npoints/2);             
             this.scan.loops(1).ramptime = -0.025; 
             this.scan.loops(2).rng = [yvals(indY)-0.007,yvals(indY) + 0.007]; 
             this.scan.loops(2).npoints = floor(this.scan.loops(2).npoints/2);             
             try
                this.run;
             end
             this.scan = scan; 
         end
        end
    end
end