
classdef Chrg < autotune.Op
    %Chrg < autotune.Op % charge scan for autotune.
    %   see help for individual properties and methods.
    %   e.g. >> help autotune.Chrg.plsgrp
    %   e.g. >> help autotune.Chrg.Run
    properties
        scan;% scan for charge scan;
        pls = 'chrg_1_L'; %plsgrp for charge scan (so markers are on)
        defaultOffset = 1e-3*[1,1]; %offset from blTriple to default measPt
        subPlot = nan; %chrg gets its own plot
        imgr;
        imgl;
    end
    
    properties (SetAccess= {?autotune.Data, ?autotune.Op})
        blTriple = [-2e-3, -2e-3]; % a Nrun x 2 double of bottom left triple pt locations
        trTriple = [2e-3,2e-3]; % Nrun x 2 double of top right triple pt locations
        xLeadSlope; % Nrun double of slope of x-lead
        yLeadSlope; % Nrun double of slope of x-lead
        fitTrip;
    end
    
    properties (Constant)
        figHandle = 1; %for plotting the derivative of chrg diagram
    end
    
    methods
        function this = Chrg
            global tuneData
            if strcmp(tuneData.activeSetName,'right')
                this.pls ={'chrg_1_R'};
            end
            this.makeNewScan;
        end
        
        function out = getData(this,runNumber)
            out = struct();
            out.blTriple = this.blTriple(runNumber,:);
            out.trTriple = this.trTriple(runNumber,:);
            out.xLeadSlope = this.xLeadSlope(runNumber);
            out.yLeadSlope = this.yLeadSlope(runNumber);
        end
        
        function makeNewRun(this,runNumber)
            if runNumber >1 && (runNumber ~= size(this.blTriple,1)+1 || runNumber ~= size(this.trTriple,1)+1)
                error('runNumber not consistent with know chrg runs');
            end
            this.blTriple(end+1,:) = [nan,nan];
            this.trTriple(end+1,:) = [nan,nan];
            this.xLeadSlope(end+1) = nan;
            this.yLeadSlope(end+1) = nan;
        end
        
        function run(this,opts,runNumber)
            if ~exist('opts','var')
                opts = '';
            end
            global tuneData;
            if ~exist('runNumber','var') || isempty(runNumber)
                runNumber = tuneData.runNumber + 1;
                tuneData.newRun();
            end
            this.scan.consts(2).val=awgseqind(this.pls);
            %scan = smscanpar(tuneData.chrg.scan, tuneData.cntr)
            file = sprintf('%s/sm_chrg%s_%04i', tuneData.dir, upper(tuneData.activeSetName(1)),runNumber);
            awgcntrl('start on amp');
            scan = this.scan;
            
            if isopt(opts,'wide')
                this.scan.loops(1).rng = [-0.03 0.03];
                this.scan.loops(2).rng = [-0.03 0.03];
            end
            if isopt(opts,'fine')
                this.scan.loops(1).npoints =256;
                this.scan.loops(2).npoints =128;
                this.scan.loops(1).ramptime = 2 * this.scan.loops(1).ramptime;
            end
            %clearMask
            data = smrun(this.scan, file);
            %if any(isnan(data{1}(:))); return; end
            this.ana(opts,data{1});
            this.scan = scan;
        end
        
        function ana(this,opts,data)
            global tuneData
            if ~exist('opts','var') || isempty(opts), opts = '';             end
            if ~exist('data','var') || isempty(data) || ischar(data) || numel(data)==1
                if (~exist('data','var') || isempty(data)) && ~isopt(opts,'last')
                    [data,scan]=loadAna('sm_chrg*');
                elseif exist('data','var') && ~isempty(data) && ischar(data)
                    [data,scan]=loadAna(data);
                else
                    side = upper(tuneData.activeSetName(1));
                    if isopt(opts,'last')
                        num = tuneData.runNumber;
                    else
                        num = data; 
                    end
                    fileName = sprintf('sm_chrg%s_%03.f.mat',side,num);
                    [data,scan]=loadAna(fileName);
                    if isempty(data)
                        fileName = sprintf('sm_chrg%s_%04.f.mat',side,num);
                        [data,scan]=loadAna(fileName);
                    end
                    if isempty(data)
                        return 
                    end
                end
            else
                scan = this.scan;
            end            
            [trip,slp,autofit] = atChargeAna(scan,data,opts); %#ok<*PROPLC>
            this.fitTrip = autofit;
            runNumber = tuneData.runNumber;
            if ~any(isnan(trip))
                this.xLeadSlope(runNumber) = 1./slp(4);
                this.yLeadSlope(runNumber) = 1./slp(2);            
                this.trTriple(runNumber,:) = trip(1:2);
                this.blTriple(runNumber,:) = trip(3:4);
                tuneData.measPt = 1/2*(this.blTriple(runNumber,:)+this.trTriple(runNumber,:)) + this.defaultOffset;
                if tuneData.sepDir(1)>tuneData.sepDir(2) % TL
                    tuneData.tl.slope = 1./slp(4);
                else % BR
                    tuneData.tl.slope = 1./slp(2);
                end
            end
        end
        
        function makeNewScan(this)
            global tuneData
            scan = defScan;  %#ok<*PROP>
            scan.saveloop = [2 3];
            scan.loops(1).rng = [-1e-2 1e-2];
            scan.loops(2).rng = [-1e-2 1e-2];
            scan.loops(1).npoints = 100;
            scan.loops(2).npoints = 64;
            scan.loops(1).ramptime = -4e-3;
            scan.cleanupfn(2).fn = @smaconfigwrap;
            scan.cleanupfn(2).args = {@sleep};
            scan.consts(2).val = awgseqind(this.pls);
            scan.loops(2).getchan = {tuneData.dataChan};
            scan.loops(1).setchan = tuneData.xyChan(1);
            scan.loops(2).setchan = tuneData.xyChan(2);
            scan.loops(1).stream = 1;
            this.scan = scan;
        end      
    end
end

function [trip, slp,autofit] = atChargeAna(scan, data,opts)
% Analyze a charge scan to find triple points and lead slopes
% function [trip, slp] = atChargeana(scan, data,opts)
% Tries to perform cross correlation to find triple points -- if R trip point to left of L, or T/B switched, goes to manual mode, or if program crashes.
% use at_template to set up tuneData.chrg.imgl and imgr

global tuneData
% Clean up data by removing outliers.
dataDiff = diff(data, [], 2);
dataDiff = dataDiff- nanmedian(dataDiff(:));
dataDiff = dataDiff .* sign(nanmean(dataDiff(:)));
m=nanmean(dataDiff(:)); s=nanstd(dataDiff(:));
dataDiff(abs(dataDiff)-m>6*s)=NaN;
plotDeriv([scan.loops(1).rng;scan.loops(2).rng],dataDiff,tuneData.sepDir);

autofit=0;
if ~isopt(opts,'man') % Try automatic triple point identification
    try
        [tripBL,tripTR] = atCorrelate(scan,dataDiff);
        autofit=1;
        if tripTR(1) <= tripBL(1)
            fprintf('Correlator found bottom left triple point right of top right\n');
            autofit=0;
        end
        if tripTR(2) <= tripBL(2)
            fprintf('Correlator found bottom left triple point above top right \n');
            autofit=0;
        end
    catch
        autofit = 0;
    end
    if isopt(opts,'auto') && autofit == 0 % if it fails and it's not required to do auto, end.  At other points, may want to comment out.
        trip = nan(1,4); slp = nan(1,4); 
        return
    end
end % Automatic triple point identification
if autofit == 0 % if triple point needs to be found, you can turn manual mode on or off by pressing m.
    plotDeriv([scan.loops(1).rng;scan.loops(2).rng],dataDiff,tuneData.sepDir);
    fprintf('Click on top right and then bottom left triple point. n to negate data.  m to begin manual fit.\n');
    [xp, yp] = ginput(2);
    tripTR = [xp(1),yp(1)];
    tripBL = [xp(2),yp(2)];
    if(tripTR(1) < tripBL(1)) % User is a doofus and clicked l-r rather than r-l
        tmp = tripBL;
        tripBL = tripTR;
        tripTR = tmp;
        fprintf('You are a doofus\n');
    end
end
if isopt(opts,'mnslp') % Fit leads, and if doing auto, extract more exact triple point from them.
    fprintf('Click on horizontal transition out of L near TL edge of plot\n');
    TLslp=ginput(1);
    fprintf('Click on vertical transition out of L near BL edge of plot\n');
    BLslp=ginput(1);
    fprintf('Click on horizontal transition out of R near BR edge of plot\n');
    BRslp=ginput(1);
    fprintf('Click on vertical transition out of R near TR edge of plot\n');
    TRslp=ginput(1);
    trip = [tripTR,tripBL];
    vertR = (TRslp(1)-tripTR(1))/(TRslp(2)-tripTR(2));
    vertL = (BLslp(1)-tripBL(1))/(BLslp(2)-tripBL(2));
    horzL = (TLslp(1)-tripBL(1))/(TLslp(2)-tripBL(2)); 
    horzR = (BRslp(1)-tripTR(1))/(BRslp(2)-tripTR(2));    
    slp = [vertR(1), vertL(1), horzR(1), horzL(1)];        
    % Plot Leads    
    plot([tripTR(1),TRslp(1)],[tripTR(2),TRslp(2)],'k-',...
        [tripTR(1),BRslp(1)],[tripTR(2),BRslp(2)],'k-',...
        [tripBL(1),BLslp(1)],[tripBL(2),BLslp(2)],'k-',...
        [tripBL(1),TLslp(1)],[tripBL(2),TLslp(2)],'k-');
    plot(trip([1 3]), trip([2 4]), '.r','MarkerSize',15);
else
    cntr = mean([tripTR;tripBL]); % % Center of junction: 1st row x, 2nd row  y
    tripDist=tripTR-tripBL; %distance between trip points
    slpJunc=tripDist(2)/tripDist(1); %junction slope
    slpPerpJunc=-1/slpJunc;
    x = linspace(scan.loops(1).rng(1), scan.loops(1).rng(2), size(dataDiff, 2)+1);
    x = (x(2:end)+x(1:end-1))/2;
    y = linspace(scan.loops(2).rng(1), scan.loops(2).rng(2), size(dataDiff, 1));
    [x2D, y2D] = meshgrid(x, y);
    % This splits the data into 4 quadrants, above/right of center of junction or not. Cartesian axis along junction and perp to it.
    above=(y2D > (x2D-cntr(1))*slpJunc + cntr(2));  % Above 11-02 transition
    right=(y2D > cntr(2)+slpPerpJunc*(x2D-cntr(1)));    % right of 11-02 transition
    aboveR=(y2D > (x2D-tripTR(1))*slpPerpJunc + tripTR(2)); % above or right of tr
    belowL=(y2D < (x2D-tripBL(1))*slpPerpJunc + tripBL(2)); % below or left of bl
    if 0 % Debug image segmentation
        rng = vertcat(scan.loops.rng);
        figure(7); clf;
        subplot(2,2,1)
        scale=nanstd(dataDiff(:));
        imagesc(rng(1, :), rng(2, :), dataDiff+scale*above);
        title('Above');
        axis xy;
        subplot(2,2,2);
        imagesc(rng(1, :), rng(2, :), dataDiff+scale*right);
        title('Right');
        axis xy;
        subplot(2,2,3)
        imagesc(rng(1, :), rng(2, :), dataDiff+scale*aboveR);
        title('Above + Right');
        axis xy;
        subplot(2,2,4);
        imagesc(rng(1, :), rng(2, :), dataDiff+scale*belowL);
        title('Below + Left');
        axis xy;
        figure(8);      clf;
        subplot(2,2,1)
        imagesc(rng(1, :), rng(2, :), dataDiff+scale*(above & aboveR));
        title('Top Right'); axis xy;
        subplot(2,2,2);
        imagesc(rng(1, :), rng(2, :), dataDiff+scale*(~above & right)); %& ~aboveR));
        title('Bottom Right');
        axis xy;
        subplot(2,2,3)
        imagesc(rng(1, :), rng(2, :), dataDiff+scale*(above & ~right));% & ~belowL));
        title('Top Left');
        axis xy;
        subplot(2,2,4);
        imagesc(rng(1, :), rng(2, :), dataDiff+scale*(~above & ~right & belowL));
        title('Bottom Left');
        axis xy;
    end
    %UR sector
    thresh =  2.5*nanstd(dataDiff(:)); %changed from 3 sigma 09/20/10
    ptsVertR=find((dataDiff(:) > thresh) & above(:) & aboveR(:));
    ptsHorzR=find((dataDiff(:) > thresh) & ~above(:) & right(:));
    ptsHorzL=find((dataDiff(:) > thresh) & above(:) & ~right(:));
    ptsVertL=find((dataDiff(:) > thresh) & ~above(:) & ~right(:) & belowL(:));
    
    robust = 1;
    if length(ptsVertR) > 2 && robust % Fit all the leads to lines. fit func has form x = fitfunc(pars,y)
        vertR = fliplr(robustfit(y2D(ptsVertR), x2D(ptsVertR))');
    else
        vertR = polyfit(y2D(ptsVertR),x2D(ptsVertR),1);
    end
    if length(ptsHorzR) > 2 && robust
        horzR = fliplr(robustfit(y2D(ptsHorzR), x2D(ptsHorzR))');
    else
        horzR = polyfit(y2D(ptsHorzR),x2D(ptsHorzR),1);
    end
    if length(ptsHorzL) > 2 && robust
        horzL = fliplr(robustfit(y2D(ptsHorzL), x2D(ptsHorzL))');
    else
        horzL = polyfit(y2D(ptsHorzL),x2D(ptsHorzL),1);
    end
    if length(ptsVertL) > 2 && robust
        vertL = fliplr(robustfit(y2D(ptsVertL), x2D(ptsVertL))');
    else
        vertL = polyfit(y2D(ptsVertL),x2D(ptsVertL),1);
    end
    % find intersections;
    diffR = vertR-horzR;
    diffL = horzL-vertL;
    
    %if autofit % Intersection of lead lines
    yR = -diffR(2)/diffR(1);
    yL = -diffL(2)/diffL(1);
    xR = yR * vertR(1) + vertR(2);
    xL = yL * vertL(1) + vertL(2);
    trip = [xR,yR,xL,yL];
    %else  % trip has format [xR,yR,xL,yL]
    %        trip = [tripTR,tripBL];
    %    end
    slp = [vertR(1), vertL(1), horzR(1), horzL(1)];
    
    % Plot lines on leads and points at trips.
    axis manual;
    plot(cntr(1), cntr(2), 'x');
    plot(polyval(vertR, [y(end), trip(2)]), [y(end), trip(2)], 'k')
    plot(polyval(horzR, [trip(2), y(1)]), [trip(2), y(1)], 'k')
    plot(polyval(horzL, [trip(4), y(end)]), [trip(4), y(end)], 'k')
    plot(polyval(vertL, [y(1), trip(4)]),[y(1), trip(4)], 'k');
    plot(trip([1 3]), trip([2 4]), '.r','MarkerSize',15);
end
end

function plotDeriv(rng,data,sepDir)
global tuneData
figure(tuneData.chrg.figHandle); clf;
imagesc(rng(1, :), rng(2, :), data);
set(gca, 'YDir', 'Normal')
axis image;
box on; hold on;

% Add 0,2 and 1,1 labels
if sepDir(2) < sepDir(1)
    label11=[max(rng(1,:)) min(rng(2,:))];
    label02=[min(rng(1,:)) max(rng(2,:))];
else
    label02=[max(rng(1,:)) min(rng(2,:))];
    label11=[min(rng(1,:)) max(rng(2,:))];
end
lambda=0.9;
label11=lambda*label11+(1-lambda)*label02;
label02=lambda*label02+(1-lambda)*label11;
text(label11(1),label11(2),'(1,1)','HorizontalAlignment','Center','FontWeight','Bold');
text(label02(1),label02(2),'(0,2)','HorizontalAlignment','Center','FontWeight','Bold');
end