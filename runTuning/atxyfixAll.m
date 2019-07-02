function basis= atxyfixAll(side, opts, vStep)
% basis = atxyfixAll(side, opts, vStep)
%   updates the XL, YL, XR, YR part of the basis by by seeing how triple points move when incremented or how stp moves.
%   side: a side name, empty to do both, left, right.
%   vStep: how much to move.
%   Note: there are many functions to create gradients.
%   For compensation matrix, use atcompfix3 (stp only).
%   For other gates, use atgradfix2 (chrg scans), atgradfix3 (stp/tl).
%   opt:
%      amp: let you do it away from measurement point = 0.
%      cent: start offset by vstep/2 in each direction
%      chrg: triple point scan.
global tuneData;
if ~exist('opts','var'), opts = ''; end
figure(10); clf; figure(11); clf;
if ~exist('vStep','var') || isempty(vStep)
    if isopt(opts,'chrg')
        vStep = 2e-3;
    else
        vStep = 1.5e-4;
    end
else
    while vStep > 1e-2
        vStep = 1e-3*vStep;
        warning('Detected huge step. Assuming you meant mV');
    end
end
if exist('side','var') && ~isempty(side)
    sides={side};
else
    sides ={tuneData.alternates.activeSetName};  % fix me
end
chans={}; rfs = {};
for i = 1:length(tuneData.alternates) %save sets    
    if find(strcmp(tuneData.alternates(i).activeSetName, sides))
        chans = [chans tuneData.xyBasis];
        rfs = [rfs tuneData.xyChan];
    end
end
if isopt(opts,'cent')
    for i =1:length(chans)
        tuneData.change(chans{i},-vStep/2)
    end
end
nGates = 19;
gtVals = cell2mat(smget(nGates));
if any(awgcntrl('israw')), awgcntrl('amp'); end
if ~awgcntrl('ison'), awgcntrl('on start wait err'); end
%% Run scans. 
% For chrg scans, click triple points and build gradient. 
if isopt(opts,'chrg')
    for s=1:length(sides)
        autotune.swap(sides{s});
        scan = smscanpar(tuneData.chrg.scan, tuneData.measPt);
        for i = 0:length(chans)
            if i == 0
                d(i+1) = smrun(scan); %#ok<*AGROW> % Run initial scan with no changes
            else
                % Change gate, run scan, change gate back.
                tuneData.change(chans{i},vStep);
                fprintf('atchg '); fprintf(chans{i}); fprintf(' by %.2f mV \n' ,vStep*1e3);
                d(i+1) = smrun(scan);
                tuneData.change(chans{i},-vStep);
                fprintf('atchg '); fprintf(chans{i}); fprintf(' by %.2f mV \n' ,-vStep*1e3);
            end
            if any(isnan(d{i+1}(:))); smset(1:19, gtVals); return; end % Fix me.
        end
    end
    if isopt(opts,'cent')
        for i =1:length(chans)
            tuneData.change(chans{i},vStep/2)
        end
    end
    for i = 1:length(d)
        dataDiff = diff(d{i}, [], 2);
        m=nanmean(dataDiff(:)); s=nanstd(dataDiff(:));
        dataDiff(abs(dataDiff)-m>3*s)=NaN;
        figure(70); clf;
        imagesc(scan.loops(1).rng, scan.loops(2).rng, dataDiff);
        set(gca,'YDir','Normal')
        x=ginput(2);
        trip(i,1:2) = mean(x);
        dtrippoint = trip(i,:)-trip(1,:);
        if i > 1
        fprintf('Triple point moves %.3f mV in x direction, %.3f mV in y direction for a change %.2f mV on gate %s \n',...
            dtrippoint(1)*1e3,dtrippoint(2)*1e3,vStep*1e3,chans{i-1});
        end
    end
    grad= (trip(2:end,:)-trip(1,:))/vStep;
    grad = grad';   
else
    for s=1:length(sides)
        autotune.swap(sides{s});
        otherSide=setdiff([1,2],s); % Grab gates from the other side. 
        dv(s,:) = [vStep, vStep];
        dv(otherSide,:) = [8*vStep, 8*vStep]; % Increment gates on other by larger amount.
        dv = dv'; dv = dv(:);
        smset(rfs,0); 
        [stp0,stpfit]=stpscan([],[],sides{s},1,'amp'); % Find initial STP point
        [tl0,tlfit]=tlscan([],[],sides{s},1,'amp'); % Initial TL point. 
        if ~stpfit || ~tlfit, return; end        
        for i=1:2 % Run stp/tl scans with amt dvCurr added to each RF gate.            
            dvCurr = dv(2*(s-1)+i);
            stpLoc = stpscan(rfs(s,i),dvCurr,sides{s},i+1,'amp');
            G(1,i) = (stpLoc-stp0)/dvCurr; 
            %fprintf('A shift of %3.f uV on %s moves the %s STP by %3.1f uV \n',dv(i)*1e6,rfs{s,i},sides{s},G(1,i)*dvCurr*1e3);
            smset(rfs(:),0); % comp has format [stpX stpY; tlX tlY]
            tlLoc = tlscan(rfs(s,i),dvCurr,sides{s},i+1,'amp'); 
            G(2,i) = (tlLoc-tl0)/dvCurr; 
            fprintf('A shift of %3.f uV on %s moves the %s STP by %3.1f uV, TL by %3.1f uV \n',dv(i)*1e6,rfs{s,i},sides{s},G(1,i)*dvCurr*1e3,G(2,i)*dvCurr*1e3);
            smset(rfs(:),0);
        end
        for i=1:length(chans)
            tuneData.change(chans{i},-dv(i));
            pause(0.2);
            [stppt,stpfit]=stpscan([],[],sides{s},i+3,'amp');
            A(1,1) = -(stppt-stp0)/dv(i);  % Negative sign is because -dv
            %fprintf('A shift of %3.f uV on %s moves the %s STP by %3.1f uV \n',dv(i)*1e6,chans{i},sides{s},A(1)*dv(i)*1e3);
            [tlpt,tlfit]=tlscan([],[],sides{s},i+3,'amp');
            tuneData.change(chans{i},dv(i))
            A(2,1)=-(tlpt-tl0)/dv(i);
            grad(2*s-1:2*s,i)=-pinv(G)*A; % Negative because X/Y do oppsoite of PlsRamp
            fprintf('A shift of %3.f uV on %s moves the %s STP by %3.1f uV, TL by %3.1f uV\n',-dv(i)*1e6,chans{i},sides{s},-A(1)*dv(i)*1e3,-A(2)*dv(i)*1e3);
            if ~stpfit
                fprintf('STP for %s %s didn''t fit \n', chans{i},sides{s})
            end
            if ~tlfit
                fprintf('TL for %s %s didn''t fit \n', chans{i},sides{s})
            end
        end
    end
end
%Now we have all the changes, and need to create a new basis from them.
%trip point has format (run,xy,side)
basis = pinv(grad);
basxy = basislookup(chans);
if length(sides) == 1 && strcmp(sides,'left')
    inds = 1:2;
elseif length(sides) == 1 && strcmp(sides,'right')
    inds = 3:4;
else
    inds = 1:4;
end
basisSave=tuneData.basis(inds,basxy)*basis;
printbasis(basis', chans, chans);
printbasis(basisSave',tuneData.gateChans(inds), chans);
doit = strcmp(input('Accept? (y/[n]', 's'), 'y'); % Make this easier
if doit
    tuneData.basis(inds,basxy)=basisSave;
    tuneData.basisStore{tuneData.runNumber}=tuneData.basis;
end
%close(10); 
end

function printbasis(basis, channames, basenames)
fprintf('       ')
for i=1:length(channames)
    fprintf(['%-4s',repmat('%-5s', 1, 5), '\n'], '',channames{i});
end
fprintf('\n');
fprintf('-----------------------------------------------------\n')
for i = 1:length(basenames)
    fprintf(['%-9s:', repmat('%8.2g', 1, 5), '\n'], [basenames{i} 'n'], basis(i,:));
    fprintf('\n')
end
end