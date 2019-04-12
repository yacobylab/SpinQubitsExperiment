function [out, histVoltages, histData, meanvals, fitpars]=procPlsData(filename,config)
% Loads a set of files and returns a struct with information about them:
% function [out, histVoltages, histData, meanvals, fitpars]=procPlsData(filename,config)
% pulsegroups, dBz group, scan, scantime, T1, scaled data, xvals
% Plots the scaled data, and pops up a PPT dialog
% options
%   gatesweep: scan is 1 group repeated while sweeping a gate
%   2d: include a 2d colorscale plot.
%   subleft, subline, subcol: handy background subtraction tricks.
%   noplot, noppt,
%   nodbz: drop dbz reference.
%   colorplot: 2d color plot of all data.
%   noscale: output raw (unrescaled) data
%   samefig
%   offset

if ~exist('filename','var') || isempty(filename), [filename,fpath]=uigetfile('sm*.mat','MultiSelect','on'); end
if ~exist('fpath','var'), fpath =pwd; end
if ischar(filename), filename={filename}; end
global tuneData;
offset=0; figs=[]; sind=1;

if ~exist('config','var') || isempty(config)
    config.opts='';
elseif ischar(config)
    config=struct('opts',config);
elseif iscell(config)
    config = struct(config{:});
end
config=def(config,'opts','samefig hold');
config=def(config,'grps',[]);
config=def(config,'xvals',[-Inf,Inf]);
config=def(config,'legend','prettyname');
config = def(config,'side',{tuneData.activeSetName});

for f=1:length(filename)
    d=load(fullfile(fpath,filename{f}));
    out(f).filename=filename{f}; out(f).scan=d.scan; %#ok<*AGROW>
    out(f).scan.data.prettyname=regexprep(filename{f},'(sm_)|(\.mat)','');
    out(f).scantime=getFileTime(fullfile(fpath,filename{f}));
    if length(d.scan.data.pulsegroups) == 1 && ismatrix(d.data{1})
        for j=1:length(d.data)
            out(f).data{j}=reshape(d.data{j},[size(d.data{j},1),1,size(d.data{j},2)]);
        end
    else
        try
            out(f).data=d.data;
        catch            
            warning('No data in file %s \n',filename{f}); 
            continue
        end
    end
    for s=1:length(config.side)
        out(f).t1(s) = att1(config.side{s},out(f).scantime,'before',d.scan);
    end
    config.dbz=find(cellfun(@(p) ~isempty(p),regexp({out(f).scan.data.pulsegroups.name},'[dD][bB][zZ]')));
    config.nodbz = setdiff(1:length(out(f).scan.data.pulsegroups),config.dbz);
    if isempty(config.grps)
        config.grps = 1:length(out(f).scan.data.pulsegroups);
        if isopt(config.opts,'nodbz'), config.grps=setdiff(config.grps,config.dbz); end
    elseif length(config.grps) == 2
        if isinf(config.grps(2)), config.grps=config.grps(1):length(out(f).scan.data.pulsegroups); end
    end
    out(f).grps=config.grps;
    if ~isopt(config.opts,'noplot')
        if isopt(config.opts,'samefig')
            figure(1); figs=unique([figs 1]);
        else
            figs=unique([figs f]); figure(f);
        end
        if isopt(config.opts,'hold')
            hold on;
        else
            clf;
        end
    end
    sz=size(out(f).data{1});
    if isfield(config,'xval') && ~isempty(config.xval)
        out(f).xv = mat2cell(config.xval);
    else
        for j=1:length(out(f).scan.data.pulsegroups)
            try
                out(f).xv{j} = d.scan.data.pulsegroups(1).varpar(:,1)';
            catch
                out(f).xv{j}=(1:1:size(out(f).data{1},3));
            end
        end
    end
    out(f).tv = guesstv(out(f).scan,config.grps);
    channels=0;
    for i=1:length(out(f).data)
        szs = size(out(f).data{i});
        if all(size(szs) == size(sz)) && all(szs == sz) % This is apparently data.
            channels=channels+1;
        end
    end
    uchan=0;
    if isopt(config.opts,'linescale')
        [out(f).data, ~, meanvals, fitpars, histVoltages, histData]=anaHistScaleLine(out(f).scan,out(f).data,out(f).t1);%vvv and n are the histogram data
    elseif ~isopt(config.opts,'noscale') % Scale data
        [out(f).data, ~, meanvals, fitpars, histVoltages, histData, fidelity]=anaHistScale(out(f).scan,out(f).data,out(f).t1);%vvv and n are the histogram data
    end
    out(f).fidelity = fidelity; 
    out(f).meanvals = meanvals; 
    for i=1:length(out(f).data)
        szs = size(out(f).data{i});
        if all(size(szs) == size(sz)) && all(szs == sz) % This is apparently data.
            if ~isopt(config.opts,'noplot')
                rdata=reshape(permute(out(f).data{i},[1 3 2]),szs(1),szs(2)*szs(3));
                figure(10+i); imagesc(rdata);
            end
            uchan=uchan+1;
            if ~isopt(config.opts,'noplot'), figure(1); subplot(1,channels,uchan); end
            if isopt(config.opts,'gatesweep')
                legs=linspace(out(f).scan.loops(1).rng(1),out(f).scan.loops(1).rng(2),out(f).scan.loops(1).npoints);
                out(f).legs=legs;
            else
                legs=out.tv;
                if isempty(legs)
                    legs=1:size(out(f).data{i},2);
                end
            end
            if isempty(config.grps), config.grps=1:length(out(f).scan.data.pulsegroups); end
            for k=config.nodbz
                if isfield(config,'frames') && ~isempty(config.frames)
                    out(f).d{i} = squeeze(nanmean(out(f).data{i}(config.frames,k,:),1));
                else
                    out(f).d{i} = squeeze(nanmean(out(f).data{i}(:,k,:),1));
                end
                if isopt(config.opts,'gatesweep')
                    leg=legs(k);
                else
                    leg=out(f).scan.data.prettyname;
                end
                if isnumeric(leg), leg=sprintf('%g',leg); end
                if isopt(config.opts,'center')
                    fo = mean(out(f).d{i});
                else
                    fo=0;
                end
                if ~isopt(config.opts,'noplot')
                    hold on; plot(out(f).xv{min(k,end)},out(f).d{i}+offset-fo,'DisplayName',leg);
                end
                sind=sind+1;
                if isopt(config.opts,'offset'), offset=offset + mean([std(out(f).d{i}),range(out(f).d{i})]); end
            end
            if ~isopt(config.opts,'noplot'), legend show; end
            if ~isopt(config.opts,'nodbz') % Plot dBz data
                for k=config.dbz
                    figure(2); hold on;
                    title('dBz groups');
                    if isfield(config,'frames') && ~isempty(config.frames)
                        out(f).d{i} = squeeze(nanmean(out(f).data{i}(config.frames,k,:),1));
                    else
                        out(f).d{i} = squeeze(nanmean(out(f).data{i}(:,k,:),1));
                    end
                    if isfield(out(f).scan.data,config.legend)
                        leg=out(f).scan.data.(config.legend);
                    elseif isopt(config.opts,'gatesweep')
                        leg=legs(k);
                    else
                        leg=out(f).scan.data.prettyname;
                    end
                    if isnumeric(leg), leg=sprintf('%g',leg); end
                    if isopt(config.opts,'center')
                        fo = mean(out(f).d{i});
                    else
                        fo=0;
                    end
                    if ~isopt(config.opts,'noplot')
                        hold on; plot(out(f).xv{min(k,end)},out(f).d{i}+offset-fo,'DisplayName',leg);
                    end
                    sind=sind+1;
                    if isopt(config.opts,'offset'), offset=offset + mean([std(out(f).d{i}),range(out(f).d{i})]); end
                end
            end
            if isopt(config.opts,'2d')
                figure(99); clf;
                z=squeeze(nanmean(out(f).data{i},1));
                z=z(config.grps,:);
                if isopt(config.opts,'subline'), z=z-repmat(mean(z,2),1,size(z,2)); end
                if isopt(config.opts,'subleft'), z=z-repmat(z(:,1),1,size(z,2)); end
                if isopt(config.opts,'subcol'), z=z-repmat(mean(z,1),size(z,1),1); end
                if isopt(config.opts,'plane')
                    coeff=fit_plane(z);
                    [mx,my]=meshgrid(1:size(z,2),1:size(z,1));
                    z=z-mx*coeff(1)-my*coeff(2)-coeff(3);
                end
                if isfield(config,'smooth'), z=filter(z,config.smooth); end
                imagesc(out(f).xv{end},legs,z);
                out(f).z=z;
            end
        end
    end
    if ~isopt(config.opts,'noplot') % Show Legend
        legend('off'); out(f).l=legend('show');
        set(out(f).l,'Interpreter','none');
    end
end
if 0%~isopt(config.opts,'noppt') % Pop up PPT dialogue
    ppt=guidata(pptplot);
    set(ppt.e_file,'String',filename{1});
    set(ppt.e_figures,'String',['[',sprintf('%d ',figs),']']);
    set(ppt.e_title,'String',out(1).scan.data.prettyname);
    set(ppt.e_body,'String','');
    set(ppt.exported,'Value',0);
end
end

function tv=guesstv(scan,grps) % this is variation between groups.
if ~isfield(scan.data,'pulsegroups') || length(scan.data.pulsegroups)==1
    tv = 0 ;
else
    if isfield(scan.data.pulsegroups,'xval') && ~isempty(scan.data.pulsegroups(1).xval)
        for j=1:length(grps)
            tv(j) = scan.data.pulsegroups(j).xval;
        end
    else
        for j=1:length(grps)
            params(:,j) = scan.data.pulsegroups(j).params;
        end
        dxv=sum(diff(params,[],1) ~= 0,1);
        [dm,di]=max(dxv);
        if dm == 0
            di=1;
            fprintf('Warning: no xval variation\n');
        end
        tv=params(:,di);
    end
end
end

function se=cull(data)
m = median(data(:));
s = median(abs(data(:)-m));
se = find(abs(data(:)-m) < 2*s);
m = mean(data(se));
se = find(abs(data(:)-m) < 2*s);
end
function coeff=fit_plane(data)
data=data(~any(isnan(data),2),:);
[gx,gy] = gradient(data);
for l=1:size(gx,1)
    gx(l,:)=smooth(gx(l,:),3);
end
for l=1:size(gy,2)
    gy(:,l)=smooth(gy(:,l),3);
end
coeff(1)=median(gx(cull(gx)));
coeff(2)=median(gy(cull(gy)));
coeff(3)=mean(mean(data));
end
function out=filter(data, sigma)
if (~exist('sigma','var'))
    sigma=3;
end
wid=sigma;
if length(wid) == 1
    wid=[wid wid];
end
ks=max(5,floor(max(sigma)*3/2)*2+1);
kw=floor(ks/2);
[x,y]=meshgrid(-kw:kw,-kw:kw);
kernel=exp(-(x.*x)./(2*wid(1)*wid(1)) - (y.*y)/(2*wid(2)*wid(2)));
kernel=kernel / sum(kernel(:));
out=filter2(kernel,data);
remove=3;
out([(1:remove) (end-remove:end)],:)=nan;
out(:,[(1:remove) (end-remove:end)])=nan;
end