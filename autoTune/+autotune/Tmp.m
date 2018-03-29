classdef Tmp < autotune.Op
    %Tmp < autotune.Op % represents tuning the measPt from Stp & TL scans
    %   see help for individual properties and methods.
    %   e.g. >> help autotune.Tmp.fineIndex
    %   e.g. >> help autotune.Stp.Run
    % To work, needs slp to be defined. 
    
    properties
        quiet = 1; %should we be chatty?
        subPlot = nan; %doesn't really need one. just a place holder    
    end
    
    properties (SetAccess= {?autotune.Data, ?autotune.Op})
       fineIndex; % for multiple stp-tls in same run 
       lastCenter;
    end
    
    methods
        function out = getData(this,~)
            %ignores runNumber. just returns findIndex.
            out = struct('fineIndex',this.fineIndex);
        end
        
        function makeNewRun(this,~)
            this.fineIndex = 1;
        end
        
        function run(this,ctrl) %fixme: do some error checking
            if ~exist('ctrl','var'), ctrl = '';       end            
            global tuneData;
            runNumber = tuneData.runNumber; 
            
            % Determine direction along exch direction using STP                        
            juncSlp = tuneData.chrg.trTriple(runNumber,:) - tuneData.chrg.blTriple(runNumber,:); juncSlp = juncSlp/norm(juncSlp);
            stpSlp = tuneData.stp.slope; stpSlp = stpSlp/norm(stpSlp); % Usually exchange slope. 
            dSTP =sqrt(2)* 1e-6*(tuneData.stp.location(runNumber)-tuneData.stp.target); % in units of volts. 
            dSTPeps = dSTP*stpSlp; % desired change            
            dSTPeps = dSTPeps - dot(dSTPeps,juncSlp)*juncSlp; % We can only determine distance along vector perp to junc, so subtract other component.            
            if ~isopt(ctrl,'quiet'), fprintf('Change for stp was %3.3f,%3.3f mV\n',1e3*dSTPeps(1),1e3*dSTPeps(2)); end                        
            
            %Determine direction along tl scan direction. 
            if tuneData.sepDir(1)>tuneData.sepDir(2) % TL
                leadSlp = [1 tuneData.chrg.xLeadSlope(runNumber)]; % the slope is the x slope                
            else % BR
                leadSlp = [1 tuneData.chrg.yLeadSlope(runNumber)];                
            end
            leadSlp = leadSlp / norm(leadSlp);
            tlSlp = [1 -1./tuneData.tl.slope];  tlSlp=tlSlp/(norm(tlSlp)); % tl.slope is the lead slope.       
            dTL = 1e-6*(tuneData.tl.location(runNumber)-tuneData.tl.target); % in units of volts.
            dTLeps = dTL * tlSlp; % desired change
            dTLeps = dTLeps - dot(leadSlp,dTLeps)*leadSlp; % We can only determine distance along vector perp to lead, so subtract parallel comp. 
            if ~isopt(ctrl,'quiet'), fprintf('Change for tl was %3.3f,%3.3f mV\n',1e3*dTLeps(1),1e3*dTLeps(2)); end                        
            
            leadPerptoJunc = leadSlp - dot(leadSlp,juncSlp)*juncSlp; %
            juncPerptoLead = juncSlp - dot(leadSlp,juncSlp)*leadSlp;
            dEps0 =( leadSlp * dot(dSTPeps,leadPerptoJunc) + juncSlp * dot(dTLeps,juncPerptoLead)) / (1 - dot(leadSlp,juncSlp)^2);
            
            leadPerp = leadSlp/leadSlp(1); leadPerp(2) = -1./leadPerp(2); leadPerp = norm(leadPerp); 
            ang = acos(dot(leadSlp,juncSlp))-pi/2; 
            dSTPalt = norm(dSTPeps); 
            dTLalt = norm(dTLeps); 
            dEpsAlt = ((dTLalt - cos(ang)*dSTPalt)*leadPerp + (dSTPalt - cos(ang)*dTLalt)*leadSlp)/sin(ang)^2; 
            dEpsAlt = dTLalt*leadPerp + (dSTPalt -sin(ang)*dTLalt)*leadSlp/cos(ang); 
            dEpsAlt2 = dTLalt*leadPerp + leadSlp/sin(ang)*(dSTPalt-cos(ang)*dTLalt);
            
            ms = stpSlp(2)./stpSlp(1); mt = tlSlp(2)/tlSlp(1); mj = juncSlp(2)/juncSlp(1); ml = leadSlp(2)/leadSlp(1); 
            dEpsPar = 1/(mt-ms)*(dTLalt/tlSlp(1)*[-ms,1] - dSTPalt/stpSlp(1)*[-mt,1]);
            dEps= 1 / (mj-ml)*(dTL*tlSlp(1)*(mt - ml)*[1,mj]-dSTP*stpSlp(1)*(ms-mj)*[1,ml]);
            %dEps = dEps/2; 
            %dEps = dEpsAlt; 
            big = any(dEps > 2e-3);
            while(big) % big change
                switch(input('Accept change (y/n)? ','s'))
                    case 'y'
                        big=0;
                    case 'n'
                        return;
                end
            end
            tuneData.measPt = tuneData.measPt + dEps;
            fprintf('New measurement point %f,%f\n',1e3*tuneData.measPt(1),1e3*tuneData.measPt(2));
            this.fineIndex = this.fineIndex+1;
            figure(autotune.Chrg.figHandle); hold on;
            plot(tuneData.measPt(1),tuneData.measPt(2),'kv','markersize',10);            
            figure(autotune.Zoom.figHandle); hold on;
            plot(tuneData.measPt(1),tuneData.measPt(2),'kv','markersize',10);
        end
    end    
end