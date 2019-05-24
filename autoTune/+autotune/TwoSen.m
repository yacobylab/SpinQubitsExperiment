classdef TwoSen < autotune.Op
    %Zoom < autotune.Op % represents a zoom scan to find Pauli Blockade
    %   see help for individual properties and methods.
    %   e.g. >> help autotune.Zoom.plsGrp
    %   e.g. >> help autotune.Zoom.Run
    
    properties
        scan; % scan for taking data
        subPlot = nan;
        
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
        end
        
        function run(this)
            % TwoD Scan.           
            global tuneData;
           
            file = sprintf('%s/sm_SD%s_%04i', tuneData.dir, upper(tuneData.activeSetName(1)), tuneData.runNumber);
            d=smrun(tuneData.twoSen.scan,file);
            diffData = diff(d{1},1,2)/(xvals(2)-xvals(1)); % Take first order diff across row
            xvals = scanRng(this.scan,1);                        
            xDiff = (xvals(1:end-1)+xvals(2:end))/2;
            yvals = scanRng(this.scan,2);
            
            [maxVal,indX] = max(abs(diffData),[],2); 
            [maxDiff,indY] = max(abs(maxVal));
            indX = indX(indY);
            
            smset(this.scan.loops(1).setchan{1},xDiff(indY));
            smset(this.scan.loops(2).setchan{1},yvals(indX));
            fprintf('Setting gates to point of max sensitivity, %2.2f. %s to %4.4f. %s to %4.4f',...
                maxDiff, this.scan.loops(1).setchan{1}, xvals(indY), this.scan.loops(2).setchan{1}, yvals(indX));
        end
    end
end