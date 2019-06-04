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
    autotune.swap(tuneData.alternates(i).activeSetName);
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
gtvals = cell2mat(smget(1:18));
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
            if any(isnan(d{i+1}(:))); smset(1:19, gtvals); return; end % Fix me.
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
        otherSide=setdiff([1,2],s);
        dv(s,:) = [vStep, vStep];
        dv(otherSide,:) = [8*vStep, 8*vStep];
        dv = dv';  dv = dv(:);
        [stp0,stpfit]=stpscan([],[],sides{s},1,opts); if isnan(stp0); smset(1:19,gtvals); return; end
        [tl0,tlfit]=tlscan([],[],sides{s},1,opts);  if isnan(tl0); smset(1:19,gtvals); return; end
        if ~stpfit || ~tlfit
            fprintf('STP or TL didn''t fit. Try retuning and try again \n')
            return
        end
        for i=1:2
            G(1,i) = (stpscan(rfs(s,i),dv(2*(s-1)+i),sides{s},i+1,opts)-stp0)/dv(2*(s-1)+i); smset(rfs(:),0); % comp has format [stpX stpY; tlX tlY]
            G(2,i) = (tlscan(rfs(s,i),dv(2*(s-1)+i),sides{s},i+1,opts)-tl0)/dv(2*(s-1)+i); smset(rfs(:),0);
        end
        for j=1:length(chans)
            tuneData.change(chans{j},dv(j));
            [stppt,stpfit]=stpscan([],[],sides{s},j+3,opts);
            [tlpt,tlfit]=tlscan([],[],sides{s},j+3,opts);
            tuneData.change(chans{j},-dv(j))
            A=[stppt-stp0 ; tlpt-tl0]/dv(j);
            grad(2*s-1:2*s,j)=-pinv(G)*A;
            fprintf('A shift of %.2f mV on %s moves the STP and TL on %s by %3.1f uV and %3.1f uV, respectively \n',dv(j)*1e3,chans{j},sides{s},A(1)*dv(j)*1e3,A(2)*dv(j)*1e3);
            if ~stpfit
                fprintf('STP for %s %s didn''t fit \n', chans{j},sides{s})
            elseif ~tlfit
                fprintf('TL for %s %s didn''t fit \n', chans{j},sides{s})
            end
        end
    end
end
%Now we have all the changes, and need to create a new basis from them.
%trip point has format (run,xy,side)
basis=pinv(grad);
basxy = basislookup(chans);
if length(sides)==1 && strcmp(sides,'left')
    inds = 1:2;
elseif length(sides)==1 && strcmp(sides,'right')
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

function b=basislookup(basis)
% Give this either a single gate in string form or a cell.

global tuneData;
if ischar(basis)
    basis={basis};
end
if iscell(basis)
    for j = 1:length(basis)
        b(j)=find(strncmp(basis{j},tuneData.baseNames,length(basis{j})));
    end
else
    fprintf('Please input a cell or a string')
    b=nan;
end
end