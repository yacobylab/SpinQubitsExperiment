%function anaRamseyRF
%%  Fit frequencies 
%num = [3402:3416];
%f=fileNumbers(num);

f = getFiles;

[jFit,freq,data,a]=fitFreqs(f,'filter'); 
%% Replot data, using only certain inds 
figure(1202); clf;
extPlot=0;
opts = 'filter'; 
if isopt(opts,'filter'), extPlot = extPlot + 1; end
if isopt(opts,'amp'), extPlot = extPlot+1; end
if isopt(opts,'t2s'), extPlot = extPlot + 1; end
ax = tightSubplot(2+extPlot,'vert nox');
for i =1:length(ax)
    hold(ax(i),'on');
    if i < 2+extPlot, ax(i).XAxis.Visible='off'; end
end

inds = [1:6];
for i = inds    
    ind = 1;
    if isopt(opts,'filter')
        ylabel(ax(ind),'J (MHz)');
        plot(ax(ind),freq{i}*1e-9,a(i).jCorr,'.-'); ind = ind+1;
    end
    if isopt(opts,'amp')
        ylabel(ax(ind),'Amplitude');
        plot(ax(ind),freq{i}*1e-9,a(i).params(:,2),'.-');
        plot(ax(ind),freq{i}*1e-9,a(i).params(:,1),'.-');
        ind = ind+1;
    end
    if isopt(opts,'t2s')
        ylabel(ax(ind),'T_2^*');
        plot(ax(ind),freq{i}*1e-9,1./a(i).params(:,6),'.-'); ind = ind+1;
    end
    
    plot(ax(ind),freq{i}*1e-9,jFit{i},'.-','DisplayName',num2str(a(i).pow));
    legInd = ind;
    ylabel(ax(ind),'J (MHz)'); ind = ind+1;
    
    plot(ax(ind),freq{i}*1e-9,a(i).jErr,'.-');
    ylabel(ax(ind),'J Error (MHz)');
end
legend(ax(legInd),'show','location','best');
fitAxis(ax(1)); fitAxis(ax(2)); fitAxis(ax(3)); 
xlabel(ax(end),'Frequency (GHz)');
%% Plot data at a given frequency
figure(1100); clf; hold on;
for i = 1:1:length(a)
    ind = closePt(freq{i},1.446e9);    
    pow(i) = a(i).configvals(strcmp(a(i).configch,'RFpow3'));
    plot(data{i}(ind,:),'DisplayName',num2str(pow(i)));
    jRes(i) = jFit{i}(ind); 
end
legend(gca,'show','location','best');
%% Patch together scans at different frequencies; 
ind = 1;
figure(1103); clf; hold on; ax = gca; 
for i = 1:4 
    for j = 1:4
        freq{ind}=scanRng(a(ind).scan,1);    
        data = squeeze(a(ind).data{1}); 
        plot(freq{ind},b{ind},'.-'); 
        ind = ind+1;         
    end
    ax.ColorOrderIndex = 1; 
end
%% Just plot, don't analyze RFs
%num = 3375:3380; % 8.5
%num = ([3418:3419,3426:3430]); % 8.5
%f=fileNumbers(num);
f = getFiles; 

a=procPlsData(f,'noplot');

figure(1203); clf;
ha = tightSubplot(length(a),'smart');

for i = 1:1:length(a)
    if ndims(a(i).data{1})>2
        data{i} = squeeze(nanmean(a(i).data{1})); 
    else
        data{i} = squeeze(a(i).data{1});
    end    
    freq{i} = scanRng(a(i).scan,1);    
    pow(i) = a(i).configvals(strcmp(a(i).configch,'RFpow3'));    
    imagesc(ha(i),1:size(data{i},2),freq{i}*1e-9,data{i}); title(ha(i),num2str(pow(i)));
    ha(i).YDir = 'normal'; 
end
%% Numbers of scans (old) 
num = 3052:3056; % 7.5 - 8.5
num = 3046:3050; % 6.5 - 7.5
num = 3041:3045; % 5.5-6.5
num = 3036:3040; % 4.5
num = 3031:3035; % 3.5
num = 3025:3029; % 2.5
num = 3008:3022;
num = 3018:3022;