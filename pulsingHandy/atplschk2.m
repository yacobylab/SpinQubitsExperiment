function atplschk2(grp,opts,config,side,ind)
% plots pulse on top of relevant charge diagrams.
%function atplschk2(grp,opts,config,side,ind)
% pg can be the name of the group or the number in awggroups, if left empty you can load a file. 
%   opts can include noflat, diff, meas
% config can be a config struct with fields:
%   pulses is a list of pulses to plot.  defaults to 1.
%   run is a tunedata run.  Defaults to last.
%   awg is which awg to render.  Defaults to 1
%   offset is an offset to apply to the CHARGE DIAGRAM (ie to check zoom pulse) defaults to zero.  1x2 or 1x2 vector, xl,yl,xr,yr
%   title: title for plot. 
%   subplot 
%   axis
% For now only works for one qubit. 

global tuneData; global awgdata; global plsdata; 
purple = [0.4940 0.1840 0.5560];
orange = [0.8500 0.3250 0.0980]; 
red = [0.6350 0.0780 0.1840];
colors{2} = purple; colors{1}=orange; colors{3}=red; 
if ~exist('side','var') || isempty(side),     side = tuneData.activeSetName; end 
runname = upper(side(1)); %name of current set
if ~exist('opts','var'),     opts = ''; end
if ~exist('config','var'),     config=struct; end
if ~exist('ind','var') || isempty(ind),     ind=1; end

config=def(config,'offset',[0 0]);
config=def(config,'pulses',1);
config=def(config,'awg',1); awg=config.awg;
%% Grab pulse and scan info. 
if ~exist('grp','var') || isempty(grp) % Get pulse information 
    [f,fp] = uigetfile(''); 
    d=load(fullfile(fp,f)); 
    plsPlotInfo = d.scan.data.pulsegroups(ind);
    grp = '';
    scan = d.scan; 
elseif ischar(grp) 
    d = load([plsdata.grpdir 'pg_' grp]); 
    plsPlotInfo = d.grpdef; 
else
    plsPlotInfo=awgdata(1).pulsegroups(grp);
end
if isempty(grp) %Find measurement point at time of scan. 
    for i=1:2
        for j = 1:length(scan.consts)
            if strcmp(tuneData.xyChan{i},scan.consts(j).setchan)
                measPt(i) = d.scan.consts(j).val;
            end
        end
    end
else
    measPt = tuneData.measPt;
end
if isopt(opts,'meas'),     config.offset = -measPt; end

f = figure(11); 
if isopt(opts,'clf'), clf; end
f.Name = 'Pulse Check 2';
if isfield(config,'subplot') % Set up figure for pulses 
    subplot(config.subplot(1),2*config.subplot(2),2*config.subplot(3)-1); a2 = gca; hold on;
    subplot(config.subplot(1),2*config.subplot(2),2*config.subplot(3)); a3 =gca; hold on;
elseif isfield(config,'axis')
    axes(config.axis(2)); a2 = gca; 
    axes(config.axis(3)); a3 = gca; 
else 
    subplot(3,1,1); a2 =gca; 
    subplot(3,1,2); a3=gca; 
    subplot(3,1,3); a4=gca; 
end
f = figure(10); 
if isopt(opts,'clf'), clf; end
f.Name = 'Pulse Check'; 
if isfield(config,'subplot') % Set up figure for charge scan
    subplot(config.subplot(1),config.subplot(2),config.subplot(3)); 
elseif isfield(config,'axis') 
    axes(config.axis(1));
end
a1 = gca; 
%% Load and plot charge scans
if ~isfield(config,'run') || isempty(config.run) % Grab most recent charge scan if no num given. 
    tuneFiles = dir(tuneData.dir); tuneFiles = {tuneFiles.name};
    try
        fpat = sprintf('sm_chrg%s_(\\d{4}).mat', runname);
    catch
        fpat = sprintf('sm_chrg%s_(\\d{3}).mat', runname);
    end
    fnumsCell = regexp(tuneFiles,fpat,'tokens');
    fnumsCell = [fnumsCell{:}];
    fnums = str2double([fnumsCell{:}]);
    num = max(fnums);
else
    num = config.run;
end
file = sprintf('sm_chrg%s_%04d.mat',runname,num);
fileName = [tuneData.dir '/' file];
chrgScan = load(fileName);
rngX=chrgScan.scan.loops(1).rng;
rngY=chrgScan.scan.loops(2).rng;
data=chrgScan.data{1};
if isopt(opts,'diff')
    data=diff(data);
end
if ~isopt(opts,'noflat')
    coeff=fitPlane(data);
    [mx,my]=meshgrid(1:size(data,2),1:size(data,1));
    data=data-mx*coeff(1)-my*coeff(2)-coeff(3);
end
imagesc(1e3*(rngX+config.offset(1)),1e3*(rngY+config.offset(2)),data);
axis image; set(gca,'YDir','Normal'); hold on;
%% Set up pulsegroups
plsPlotInfo=plsmakegrp(plsPlotInfo,'',config.pulses);
for i = 1:length(config.pulses)
    if(isfield(plsPlotInfo,'pulseind'))
        wfInfo=plstowf(plsPlotInfo.pulses(plsPlotInfo.pulseind(i)),plsPlotInfo.dict);
    else
        wfInfo=plstowf(plsPlotInfo.pulses(i),plsPlotInfo.dict);
    end
    outchans=cell(2,1); outmark=cell(2,1);
    for j=1:2
        outchans{j}=wfInfo.data(awg).wf(j,:);
       
    end

    if  ~isempty(outchans{1}) && ~isempty(outchans{2})        
        plot(a1,outchans{1},outchans{2},'Color',colors{i},'LineWidth',2);
        if isfield(config,'title') && ~isempty(config.title), title(config.title); end
    end    
    plot(a2,outchans{1},'DisplayName',sprintf('X%d',config.pulses(i))); hold(a2,'on');
    if isfield(config,'title') && ~isempty(config.title), title(config.title); end    
    plot(a3,outchans{2},'DisplayName',sprintf('Y%d',config.pulses(i))); hold(a3,'on'); 
    if isfield(config,'title') && ~isempty(config.title), title(config.title); end
    for j = 1:size(wfInfo.data(awg).marker,1) 
         outmark{j}=wfInfo.data(awg).marker(j,:); 
         plot(a4,outmark{j},'DisplayName',sprintf('M%d',j)); hold(a4,'on'); 
    end
end
legend(a2,'show'); legend(a3,'show'); 
end