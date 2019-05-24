function atplschk(grp,opts,config,side)
% Plots pulse on top of relevant charge diagram, as well as line plot of voltage vs. time
% for channels and markers.
% function atplschk(grp,opts,config,side,ind)
% grp: name of the group, number in awggroups, struct containting pulsegroup, 
%   or if left empty you can load scan that has pulsegroup in scan.data. 
% opts can be noflat, diff, meas
%   noflat: Don't do plane subtraction of charge scan.
%   diff: Plot differentiated charge scan.
%   meas: (CHECK) Set measPt to 0,0.
% config can be a config struct with fields:
%   pulses is a list of pulses to plot.  Defaults to 1.
%   run is a tunedata run. Defaults to last.
%   awg is which awg to render.  Defaults to 1
%   offset is an offset to apply to the CHARGE DIAGRAM (ie to check zoom pulse). Defaults to zero.  1x2 vector, xl,yl,xr,yr
%   title: title for plot.
%   subplot: Provide number of subplots
%   axis: use this to configure plotting externally, axis handles
%   timeOffset: offset the time of different pulses; useful with fill. 
% side: which qubit
% ind: 
% for multiple pulses, 1st purple, then orange, then red. Add ability to do > 3 
% For now only works for one qubit.

global tuneData; global awgdata; global plsdata;
purple = [0.4940 0.1840 0.5560];
orange = [0.8500 0.3250 0.0980];
red = [0.6350 0.0780 0.1840];
colors{1}=orange; colors{2} = purple; colors{3}=red;
if ~exist('side','var') || isempty(side), side = tuneData.activeSetName; end
runname = upper(side(1)); %name of current set
if ~exist('opts','var'), opts = ''; end
if ~exist('config','var'), config=struct; end
config=def(config,'offset',[0 0]);
config=def(config,'ind',1); ind = config.ind; 
config=def(config,'pulses',1);
config=def(config,'awg',1); awg=config.awg;
config=def(config,'timeOffset',zeros(1,length(config.pulses))); 
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
elseif isstruct(grp)
    plsPlotInfo = grp(ind);
elseif isnumeric % If it's a number, must currently be on AWG
    plsPlotInfo=awgdata(1).pulsegroups(grp);
else
    error('Incorrect format for pulse, please review help page for atplschk');
end

% Find measurement point at time of scan. If on AWG now, use current measPt.
% FIXME: when is this plotted?j
if (isempty(grp) || ischar(grp)) && exist('scan','var')
    for i=1:2
        for j = 1:length(scan.consts)
            if strcmp(tuneData.xyChan{i},scan.consts(j).setchan)
                measPt(i) = d.scan.consts(j).val;
            end
        end
    end
elseif ~isempty(grp) && ~ischar(grp)
    measPt = tuneData.measPt;
end
if isopt(opts,'meas'),  config.offset = -measPt; end
if ~isfield(config,'axis')
    f = figure(11);
    if isopt(opts,'clf'), clf; end
    f.Name = 'Pulse Check 2';
end
if isfield(config,'subplot') % Set up figure for pulses
    subplot(config.subplot(1),2*config.subplot(2),2*config.subplot(3)-1); a2 = gca; hold on;
    subplot(config.subplot(1),2*config.subplot(2),2*config.subplot(3)); a3 =gca; hold on;
elseif isfield(config,'axis') % Provide axes (for use in bulk plotting)
    axes(config.axis(2)); a2 = gca;
    axes(config.axis(3)); a3 = gca;
else
    subplot(3,1,1); a2 = gca;
    subplot(3,1,2); a3 = gca;
    subplot(3,1,3); a4 = gca;
end
if ~isfield(config,'axis')
    f = figure(10);
    if isopt(opts,'clf'), clf; end
    f.Name = 'Pulse Check';
end
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
if ~isopt(opts,'noflat') && ~isopt(opts,'diff')
    coeff=fitPlane(data);
    [mx,my]=meshgrid(1:size(data,2),1:size(data,1));
    data=data-mx*coeff(1)-my*coeff(2)-coeff(3);
end
imagesc(1e3*(rngX+config.offset(1)),1e3*(rngY+config.offset(2)),data);
axis image; set(gca,'YDir','Normal'); hold on;
%% Set up pulsegroups and plot voltage vs. time 
plsPlotInfo=plsmakegrp(plsPlotInfo,'',config.pulses);
for i = 1:length(config.pulses)
    if(isfield(plsPlotInfo,'pulseind')) % Make given pulse.
        wfInfo=plstowf(plsPlotInfo.pulses(plsPlotInfo.pulseind(i)),plsPlotInfo.dict);
    else
        wfInfo=plstowf(plsPlotInfo.pulses(i),plsPlotInfo.dict);
    end
    outchans=cell(2,1); outmark=cell(2,1);
    for j=1:2 %  2 channels of data / qubit.
        outchans{j}=wfInfo.data(awg).wf(j,:);
    end
    if  ~isempty(outchans{1}) && ~isempty(outchans{2}) % Plot on charge scan
        plot(a1,outchans{1},outchans{2},'Color',colors{i},'LineWidth',2);
        if isfield(config,'title') && ~isempty(config.title), title(a1,config.title); end
    end
    inds = (1:length(outchans{1}))+config.timeOffset(i);
    plot(a2,inds,outchans{1},'DisplayName',sprintf('X%d',config.pulses(i))); hold(a2,'on');
    if isfield(config,'title') && ~isempty(config.title), title(a2,config.title); end
    plot(a3,inds,outchans{2},'DisplayName',sprintf('Y%d',config.pulses(i))); hold(a3,'on');
    if isfield(config,'title') && ~isempty(config.title), title(a3,config.title); end
    if ~isfield(config,'axis')
        for j = 1:size(wfInfo.data(awg).marker,1)
            outmark{j}=wfInfo.data(awg).marker(j,:);
            plot(a4,inds,outmark{j},'DisplayName',sprintf('M%d',j)); hold(a4,'on');
        end
    end
end
legend(a2,'show','location','best'); legend(a3,'show','location','best');
end